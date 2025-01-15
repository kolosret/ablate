#include "szBurn.hpp"
#include "utilities/constants.hpp"
#include "particles/processes/burningModel/liquidFuels/waxFuel.hpp"
#include "particles/processes/burningModel/dropletFlame/waxFlame.hpp"
//#include "particles/processes/burningModel/liquidFuels/liquidFuel.hpp"
#include <math.h>

#include <utility>




ablate::particles::processes::burningModel::SZBurn::SZBurn(PetscReal convectionCoeff,
       PetscReal ignitionTemperature, PetscReal burnRate, PetscReal nuOx, PetscReal Lv, PetscReal heatOfCombustion,
       const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsProducts,
       PetscReal extinguishmentOxygenMassFraction, std::shared_ptr<eos::EOS> eosIn,std::string  fuelin) :
       BurningProcess(std::move(eosIn), {}, massFractionsProducts, ignitionTemperature, nuOx, Lv, heatOfCombustion, extinguishmentOxygenMassFraction),
       K(burnRate), convectionCoeff(convectionCoeff),speciesNames{fuelin, "O2", "N2", "CO2", "H2O"},speciesOffSet(5, -1)

       {

    //one can easily add more fuels here
    //TODO CHANGE THIS to C32H66 after done with debugging!!!!!!!!!
    if (fuelin == "CH4") {
        fuelType="wax";
        liquidFuel = std::make_shared<ablate::particles::processes::burningModel::waxFuel>();
        flame = std::make_shared<ablate::particles::processes::burningModel::waxFlame>();
    } else {
        throw std::invalid_argument("No liquid fuel, flamelet of the fuel with formula: " + fuelin);
    }

    //initialize farfield
    farField.Temperature = 350.0;
    farField.Pressure = 101325.0;
    farField.Yox = 1.0;


    auto speciesList = eos->GetSpeciesVariables();
    numberSpecies = speciesList.size();
    for (int idxsp = 0; idxsp < 5;idxsp++) {
        for (PetscInt idx = 0; idx < numberSpecies; idx++)
            if (speciesList.at(idx) == speciesNames[idxsp]) {
                speciesOffSet[idxsp] = idx;
                continue;
            }

    }

    //Check if all species are found in the mechanism;
    int cnt = count(speciesOffSet.begin(), speciesOffSet.end(), -1);
    if (cnt > 0)
        throw std::runtime_error("The Shvab-Zeldovich model requires fuelFormula,O2,N2,CO2,H2O to work.");

    //Initilaize droplet to 300 and 0.1 mm
    Td=300;
    Dp=100E-6;

}


void ablate::particles::processes::burningModel::SZBurn::ComputeRHS(PetscReal time, accessors::SwarmAccessor &swarmAccessor, accessors::RhsAccessor &rhsAccessor, accessors::EulerianAccessor &eulerianAccessor)
{
    //Number of dimensions
    const auto dim = eulerianAccessor.GetDimensions();
    
    //Grab wanted fields from the eulerian field accessor
    //For the simple burning droplet model, the burnrate is only a function of Tfar, Pfar, YO2far
    auto farFieldTemperature = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::TEMPERATURE_FIELD];
    auto farFieldPressure = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::PRESSURE_FIELD];
    auto farFieldSpecies = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::YI_FIELD];
    //Get consevred fields for density for the Reynolds number calculation
    auto farFieldConserved= eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD];
    auto farFieldVelocity = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::VELOCITY_FIELD]    ;

    //Get the fields(SOL) for which we are creating sources
    //PS for now it is fine that the inertial is adding source terms for the vel and position
    //I could do the drag model here to reduce overhead....
    //PS2 When droplets are analytically integrated forward make sure you don't use inertial
    auto tempRHS = rhsAccessor[ablate::particles::ParticleSolver::ParticleTemperature];
    auto massRHS = rhsAccessor[ablate::particles::ParticleSolver::ParticleMass];


    //Grab the Particle Fields
//    auto partDensity = swarmAccessor[ablate::particles::ParticleSolver::ParticleDensity];
    auto partVelocity = swarmAccessor[ablate::particles::ParticleSolver::ParticleVelocity];
    auto particlesPerParcel = swarmAccessor[ablate::particles::ParticleSolver::ParticleNPP];
    auto partDiameterAvg = swarmAccessor[ablate::particles::ParticleSolver::ParticleDiameter];
    auto partBurning = swarmAccessor[ablate::particles::ParticleSolver::ParticleBurning];
    auto partTemperature = swarmAccessor[ablate::particles::ParticleSolver::ParticleTemperature];
    auto partCp = swarmAccessor[ablate::particles::ParticleSolver::ParticleCP];
    auto parcelMass = swarmAccessor[ablate::particles::ParticleSolver::ParticleMass];

    auto numParticles = swarmAccessor.GetNumberParticles();
    PetscReal SAtot;
    //Calculate Each Particles RHS's
    for( auto np = 0; np < numParticles; np++) {
        //TODO see if the burning toggle is set elswhere if not set here
//        if (farFieldSpecies(np,oxygenOffset))
//            partBurning(np)=1;
                
        double YiO2=1;
        UpdateFarfield(farFieldTemperature(np),farFieldPressure(np),farFieldSpecies(np,oxygenOffset));
        
        if (partBurning(np)) {
            //Calculate the droplet burn rate
            CalcBurnRate();
            
            //Convective effects on the droplet
            PetscReal slipvelmag = 0;
            for (PetscInt n = 0; n < dim; n++) {
                slipvelmag += PetscSqrtReal(farFieldVelocity(np, n) - partVelocity(np, n));
            }
            PetscReal reP= farFieldConserved(np,0) * slipvelmag * partDiameterAvg(np) / gasPhase.mug;
            PetscReal Shg = 2*(1+PetscSqrtReal(reP)*TransportConst.cubeRootSc / 3)* log(1+B)/B;

//            //Update burn rate
//            K = 8 * gasPhase.kg / (gasPhase.Cpg  * liquidFuel->fuelProperties.rhol) * log(1 + B);

            //For a single droplet
            double mdot = PETSC_PI * partDiameterAvg(np) * gasPhase.kg / gasPhase.Cpg * Shg;
            //TODO make sure this is just the derivative
            massRHS(np) -= particlesPerParcel(np) * mdot;

            //For the temperature: (pi*Dd*(mug/Pr)*Nu*(Cp(Tinf-Ts)+hc/rox)-mdotg*hfg)/Cvl/(md)
            // (mug*Cpg/Pr)=kg
            // For unity Le, Sc=Pr => Sh=Nu
            double Nug = Shg;
            double mD = parcelMass(np)/particlesPerParcel(np);
            double deltaT = farFieldTemperature(np)-Ts+(flame->Hc/gasPhase.rox)/gasPhase.Cpg;
            tempRHS(np) += (PETSC_PI*partDiameterAvg(np)*(gasPhase.kg)*Nug*deltaT-mdot*liquidFuel->fuelProperties.Hvap)/liquidFuel->fuelProperties.Cvl/(mD);
        }
        else {
            SAtot = particlesPerParcel(np)*PETSC_PI*std::pow(partDiameterAvg(np),2);
            tempRHS(np) += convectionCoeff*(farFieldTemperature(np)-partTemperature(np))*SAtot/(parcelMass(np)*partCp(np));
            }

    }
}

void ablate::particles::processes::burningModel::SZBurn::ComputeEulerianSource(PetscReal startTime, PetscReal endTime, PetscInt ndims, accessors::SwarmAccessor &swarmAccessorPreStep, accessors::SwarmAccessor &swarmAccessorPostStep, accessors::EulerianSourceAccessor &eulerianSourceAccessor)
{

    //We need to communicate the energy that was transfer to/from the eulerian field to the particle back to the eulerian field rhoEnergy
    // Get sizes from the accessors
    const auto np = swarmAccessorPreStep.GetNumberParticles();
    //The mass and energy  are in the eulerian field (rhoE is the second index (1), rho is the first index)
    auto sourceEuler = eulerianSourceAccessor[ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD];
//    auto sourceSpecies = eulerianSourceAccessor[ablate::finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD];

    //Grab Next and Old Swarm Temperatures (Assume Cp and mass does not change ever)
    auto parcelTempNew = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleTemperature];
    auto parcelTempOld = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleTemperature];
    auto parcelMassNew = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleMass];
    auto parcelMassOld = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleMass];
    //Determine if the particle was burning or not
    auto parcelBurning = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleBurning];
    auto parcelCp = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleCP];


    //Product Species will be the same for all particles so calculate them here
    //The below should work since the function should be a constant function and not dependent on space or time
    double mdot;
    //Compute the source terms due to burning/ heat transfer
    for (PetscInt p = 0; p < np; ++p) {


        if (parcelBurning(p)) {
            //Calculate the droplet burn rate
            CalcBurnRate();

            // Start by adding in the mass source terms
            mdot = parcelMassOld(p)-parcelMassNew(p); //mass going into farfield
            sourceEuler(p,ablate::finiteVolume::CompressibleFlowFields::RHO) += mdot;
            // now onto the species mass
            //Now to take care of species terms
//            for (PetscInt ns = 0; ns < numberSpecies; ns ++)
//                sourceSpecies(p,ns) += mdot*(1+properties.nuOx)*properties.massFractionsProducts[ns]; //add the products
//            sourceSpecies(p,oxygenOffset) -= properties.nuOx*mdot; //remove the oxygen
            //add in the burning energy
            sourceEuler(p,1) -= mdot*properties.HC;
        }
        else {
            //Since this is assumed no mass loss unless its burning, add in the particle heating/cooling coupling
            sourceEuler(p,1) += parcelMassNew(p)*parcelCp(p)*(parcelTempNew(p)-parcelTempOld(p));
        }
    }
}

void ablate::particles::processes::burningModel::SZBurn::UpdateFarfield(double Tfar, double Pfar, double YiO2far) {
    //Update farfield properties
    farField.Temperature=Tfar;
    farField.Pressure=Pfar;
    farField.Yox=YiO2far;

    //Sutherland gas transport, Le=1,
    gasPhase.mug = TransportConst.muo * PetscSqrtReal(Tfar / TransportConst.to) * (Tfar / TransportConst.to) * (TransportConst.to + TransportConst.so) / (Tfar + TransportConst.so);

    //For now its okay I guess...
    gasPhase.Cpg = liquidFuel->fuelProperties.Cpg;
    //Alternativelly, need to change EOS to zerork, 
    //    gasPhase.Cpg = eos->mech->getMassCpFromTY();

    //Conductivity of the gas
    gasPhase.kg = gasPhase.mug * gasPhase.Cpg  / TransportConst.pr;

    
    //For if droplet burning update flame properties
    //TODO pass burning variable in here, to stay consistent with burning criterion
    if (YiO2far>0.1) {
        double nN2dnO2 = ((1 - YiO2far) * 28) / (YiO2far * 32);
//        printf("nN2dnO2 is: %f\n", nN2dnO2);

        //For air: rox = (x+y/4)*(MWO2+3.76*MWN2)/MWF
        gasPhase.rox = (flame->x + flame->y / 4) * (flame->MWO2 + nN2dnO2 * flame->MWN2) / flame->MWFuel;

        //Flame temperature (For now its constant set at 2700, one can recalculate it for different YO2far
        gasPhase.Tflame=flame->Tad;
        gasPhase.nuOx=flame->nuO2;
        
    }



};

void ablate::particles::processes::burningModel::SZBurn::CalcBurnRate() {

    //Solve the droplet itteration
    double root = BisectionMethodSolve([this](double Y) { return SolveSZBurn(Y); },0.999999,1E-10,1e-10,100);

    //Set the transfer number
    B = result.BoxT;
    Ts = result.Ts;
}

double ablate::particles::processes::burningModel::SZBurn::SolveSZBurn(double YAsOld) {
    double MWs = 1 / (YAsOld / liquidFuel->fuelProperties.MW + (1 - YAsOld) / MWair);
    double Pvap = YAsOld * farField.Pressure * MWs / liquidFuel->fuelProperties.MW;

    // Bound the vapor pressure
    Pvap = std::min(Pvap, farField.Pressure);
    Pvap = std::max(Pvap, 1E-20);

    // Calculate the surface temperature
    double Tsnew;
    liquidFuel->Tvap(&Tsnew, &Pvap);

//    double MdotF_D_Mdot1=1;
//    double MdotOX_D_Mdot2 = 0; -nuO2
//    double gamma1 = 1; //Temperature dependent gamma = mu(T) * Sc
//    double gamma2 = 1;
//
//    double LHS = log((MdotF_D_Mdot1 - YAsOld) / MdotF_D_Mdot1) / log(MdotOX_D_Mdot2 / (MdotOX_D_Mdot2 - farField.Yox));
//    double rf_rs = 1 + LHS * gamma1 / gamma2;
//    double mdot1 = 2 * PETSC_PI * Dp * gamma1 / (1. - 1. /(rf_rs+1E-20)) * log(MdotF_D_Mdot1 / (MdotF_D_Mdot1 - YAsOld));

    double Bm = (farField.Yox/gasPhase.nuOx + YAsOld)/(1-YAsOld+1E10);
    double mdot = 2.*PETSC_PI*Dp*(gasPhase.kg/gasPhase.Cpg)*log(1+Bm)/(Bm);
    double Ad = Dp*Dp*PETSC_PI;
    double mFlux = mdot/Ad;

    double ql = -liquidFuel->fuelProperties.kl * (Tsnew - Td) / (max(0.5*Dp, 1E-10));
    double qlLimiter = abs(dropletConstants.qlimfac * mFlux * liquidFuel->fuelProperties.Hvap) / abs(ql + 1E-10);

    qlLimiter = std::min(qlLimiter,1.e+0);
    ql = ql* qlLimiter;
//    ql = 0;

    double BoxT = (farField.Yox * flame->Hc / gasPhase.rox + gasPhase.Cpg  * (farField.Temperature - Tsnew)) / (liquidFuel->fuelProperties.Hvap - ql/mFlux);

    double YFs_new = (BoxT - farField.Yox/gasPhase.rox) / (1 + BoxT);

    YFs_new = std::min(YFs_new, 1.0);
    YFs_new = std::max(YFs_new, 0.0);
    result.Yfs_new = YFs_new;

    //Update transfer number
    result.BoxT=BoxT;

    //Save surface temperature for Temperature source terms
    result.Ts=Tsnew;

    return YFs_new - YAsOld;
}

void ablate::particles::processes::burningModel::SZBurn::CalcEvapRate() {

    //Solve the droplet itteration
    double root = BisectionMethodSolve([this](double Y) { return SolveSZEvap(Y); },0.999999,1E-10,1e-10,100);

    //Set evaporation rate
    K = result.K;
}

double ablate::particles::processes::burningModel::SZBurn::SolveSZEvap(double YAsOld) {
    double MWs = 1 / (YAsOld / liquidFuel->fuelProperties.MW + (1 - YAsOld) / MWair);
    double Pvap = YAsOld * farField.Pressure * MWs / liquidFuel->fuelProperties.MW;

    // Bound the vapor pressure
    Pvap = std::min(Pvap, farField.Pressure);
    Pvap = std::max(Pvap, 1E-20);

    // Calculate the surface temperature
    double Tsnew;
    liquidFuel->Tvap(&Tsnew, &Pvap);

    double Bm = YAsOld/(1-YAsOld);
    double mdot = 2.*PETSC_PI*Dp*(gasPhase.kg/gasPhase.Cpg)*log(1+Bm)/Bm;
    double Ad = Dp*Dp*PETSC_PI;
    double mFlux = mdot/Ad;

    double ql = -liquidFuel->fuelProperties.kl * (Tsnew - Td) / (max(0.5*Dp, 1E-10));
    double qlLimiter = abs(dropletConstants.qlimfac * mFlux * liquidFuel->fuelProperties.Hvap) / abs(ql + 1E-10);

    qlLimiter = std::min(qlLimiter,1.e+0);
    ql = ql* qlLimiter;
    //    ql = 0;

    double BoxT = (gasPhase.Cpg  * (farField.Temperature - Tsnew)) / (liquidFuel->fuelProperties.Hvap - ql/mFlux);

    double YFs_new = BoxT / (1 + BoxT);

    YFs_new = std::min(YFs_new, 1.0);
    YFs_new = std::max(YFs_new, 0.0);
    result.Yfs_new = YFs_new;

    result.K = 8 * gasPhase.kg / (gasPhase.Cpg  * liquidFuel->fuelProperties.rhol) * log(1 + BoxT);

    return YFs_new - YAsOld;
}


double ablate::particles::processes::burningModel::SZBurn::BisectionMethodSolve(std::function<double(double)> func, double a, double b, double tol, int maxIter) {
    if (func(a) * func(b) >= 0) {
        std::cerr << "Error: The function must have opposite signs at the endpoints a and b." << std::endl;
        return NAN;
    }

    double c = a;
    for (int i = 0; i < maxIter; ++i) {
        c = (a + b) / 2; // Midpoint
        if (func(c) == 0.0 || abs(b - a) / 2 < tol) {
            return c; // Root found or tolerance met
//            continue;
        }
        if (func(c) * func(a) < 0) {
            b = c; // Root is in left
        } else {
            a = c; // Root is in right
        }
    }

    std::cerr << "Warning: Maximum number of iterations reached. Approximate root: " << c << std::endl;
    return c;
}

#include "registrar.hpp"

REGISTER(ablate::particles::processes::Process, ablate::particles::processes::burningModel::SZBurn, "Example of an no evaporation/Burning Model",
         OPT(PetscReal, "convectionCoefficient", "The convection Coefficient for the simple heating mode"),
         OPT(PetscReal, "IgnitionTemperature", "The temperature of the far field where it is said to *ignite*"),
         OPT(PetscReal, "burnRate", "The constant burning rate K for the d^2 law (defaults to 0)"),
         OPT(PetscReal, "nuOx", "stoichiometric oxidizer mass coefficient (defaults to 0)"),
         OPT(PetscReal, "Lv", "LatentHeat of Vaporization (defaults to 1e5)"),
         OPT(PetscReal, "Hc", "heat of combustion per unit mass of fuel (defaults to 1e6)"),
         OPT(ablate::mathFunctions::FieldFunction, "productsMassFractions", "The product species functions"),
         OPT(PetscReal, "extinguishmentOxygenMassFraction", "The minimum oxygen mass fraction needed for burning (defaults to 0.1)"),
         ARG(ablate::eos::EOS, "eos", "the eos used to compute various gas properties"),
         ARG(std::string, "fuelFormula", "Fomula of the fuel just as it is in the mechanism, curently implmented: 'C32H66'(wax)") );
