#include "twoZoneBurn.hpp"
#include "utilities/constants.hpp"
#include "particles/processes/burningModel/liquidFuels/waxFuel.hpp"
#include "particles/processes/burningModel/liquidFuels/liquidFuel.hpp"
#include <math.h>

ablate::particles::processes::burningModel::TwoZoneBurn::TwoZoneBurn(PetscReal convectionCoeff,
       PetscReal ignitionTemperature, PetscReal burnRate, PetscReal nuOx, PetscReal Lv, PetscReal heatOfCombustion,
       const std::shared_ptr<ablate::mathFunctions::FieldFunction> &massFractionsProducts,
       PetscReal extinguishmentOxygenMassFraction, std::shared_ptr<eos::zerorkEOS> eosIn,
       const std::string& fuelType) :
       BurningProcess(std::move(eosIn), {}, massFractionsProducts, ignitionTemperature, nuOx, Lv, heatOfCombustion, extinguishmentOxygenMassFraction),
       burnRate(burnRate), convectionCoeff(convectionCoeff)
       {

    if (!std::dynamic_pointer_cast<eos::zerorkEOS>(eosIn)) {
        throw std::invalid_argument("the twozone model only accepts eos::zerorkEOS as the input sorry in advance, Kolos<3");
    }
//
//    //Cast eosIn into eos
//    eos = dynamic_pointer_cast<eos::zerorkEOS>(eosIn);

    //one can easily add more fuels here
    if (fuelType == "wax") {
        liquidFuel = std::make_shared<ablate::particles::processes::burningModel::waxFuel>();
    } else {
        throw std::invalid_argument("Unknown fuel type: " + fuelType);
    }


    nSpc=eosIn->mech->getNumSpecies();


    //initialize farfield //TODO make these as optional inputs
    farField.Temperature = 350.0;
    farField.Pressure = 101325.0;
    farField.Yox = 1.0;

    // TODO evaluate the properties

}


void ablate::particles::processes::burningModel::TwoZoneBurn::ComputeRHS(PetscReal time, accessors::SwarmAccessor &swarmAccessor, accessors::RhsAccessor &rhsAccessor, accessors::EulerianAccessor &eulerianAccessor)
{
    //Grab wanted fields from the eulerian field accessor
    //In this case we only need the Oxygen and temperature fields
    auto farFieldTemperature = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::TEMPERATURE_FIELD];
//    auto farFieldSpecies = eulerianAccessor[ablate::finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD];

    //Grab The fields of the particle that we will be changing
    auto massRHS = rhsAccessor[ablate::particles::ParticleSolver::ParticleMass];
    //Burning particles
    auto tempRHS = rhsAccessor[ablate::particles::ParticleSolver::ParticleTemperature];

    //Grab the Particle Fields we need
    auto partDensity = swarmAccessor[ablate::particles::ParticleSolver::ParticleDensity];
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
        //Ideally this is a farField Value, but for this test case I'm not letting the particles burn until the particle itself reaches
        // the ignition temperature, Particle may be burning now then when it was updated in the decode call, this will be called twice
        CalcBurnRate();
        if (partBurning(np))
            massRHS(np) -= particlesPerParcel(np)*(partDensity(np)*PETSC_PI/4*partDiameterAvg(np)*burnRate);
        //If it's not burning, we are just heating with our shitty convection model
        else {
            SAtot = particlesPerParcel(np)*PETSC_PI*std::pow(partDiameterAvg(np),2);
            tempRHS(np) += convectionCoeff*(farFieldTemperature(np)-partTemperature(np))*SAtot/(parcelMass(np)*partCp(np));
            }
    }
}

void ablate::particles::processes::burningModel::TwoZoneBurn::ComputeEulerianSource(PetscReal startTime, PetscReal endTime, PetscInt ndims, accessors::SwarmAccessor &swarmAccessorPreStep, accessors::SwarmAccessor &swarmAccessorPostStep, accessors::EulerianSourceAccessor &eulerianSourceAccessor)
{

    //We need to communicate the energy that was transfer to/from the eulerian field to the particle back to the eulerian field rhoEnergy
    // Get sizes from the accessors
    const auto np = swarmAccessorPreStep.GetNumberParticles();
    //The mass and energy  are in the eulerian field (rhoE is the second index (1), rho is the first index)
    auto sourceEuler = eulerianSourceAccessor[ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD];
    auto sourceSpecies = eulerianSourceAccessor[ablate::finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD];

    //Grab Next and Old Swarm Temperatures (Assume Cp and mass does not change ever)
    auto parcelTempNew = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleTemperature];
    auto parcelTempOld = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleTemperature];
    auto parcelMassNew = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleMass];
    auto parcelMassOld = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleMass];
    //Determine if the particle was burning or not
    auto parcelBurning = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleBurning];

    auto parcelCp = swarmAccessorPostStep[ablate::particles::ParticleSolver::ParticleCP];
    PetscReal mdot;
    //Product Species will be the same for all particles so calculate them here
    //The below should work since the function should be a constant function and not dependent on space or time

    //Compute the source terms due to burning/ heat transfer
    for (PetscInt p = 0; p < np; ++p) {

        //TODO talk to kenny to aoid calculating the temperature again/ how to get aux stuff from the accessors
        double Tfar = 1;
        double Pfar = 1;
        double YiO2 = 1;

        //Update the farfield stat, T, Yox,
        UpdateFarfield(&farField, Tfar,Pfar,YiO2);


        if (parcelBurning(p)) {
            // Start by adding in the mass source terms
            mdot = parcelMassOld(p)-parcelMassNew(p); //mass going into farfield
            sourceEuler(p,ablate::finiteVolume::CompressibleFlowFields::RHO) += mdot;
            // now onto the species mass
            //Now to take care of species terms
            for (PetscInt ns = 0; ns < numberSpecies; ns ++)
                sourceSpecies(p,ns) += mdot*(1+properties.nuOx)*properties.massFractionsProducts[ns]; //add the products
            sourceSpecies(p,oxygenOffset) -= properties.nuOx*mdot; //remove the oxygen
            //add in the burning energy
            sourceEuler(p,1) -= mdot*properties.HC;
        }
        else {
            //Since this is assumed no mass loss unless its burning, add in the particle heating/cooling coupling
            sourceEuler(p,1) += parcelMassNew(p)*parcelCp(p)*(parcelTempNew(p)-parcelTempOld(p));
        }
    }
}

    void ablate::particles::processes::burningModel::TwoZoneBurn::CalcBurnRate() {

    //Calculate the combustion, either implement equlibrium or just simple flamelet
//    TwozoneEqulibrium();

    //TODO get these from the function
    double Tfar=0;
    double Pfar = 0;
    double YiO2 = 1;



    //Update the zone properties
    UpdateZoneProperties();

    //TODO implement itteration here for now just YFs=1, boiling
    SolveTwoZone(1.0, farField);

    //Calculate burnrate, Energy source,
    burnRate = 1E-7;
}

void ablate::particles::processes::burningModel::TwoZoneBurn::UpdateFarfield(ablate::particles::processes::burningModel::TwoZoneBurn::farFieldProp* farfield,double Tfar, double Pfar, double YiO2far) {

    farfield->Temperature=Tfar;
    farfield->Pressure=Pfar;
    farfield->Yox=YiO2far;

};


void ablate::particles::processes::burningModel::TwoZoneBurn::UpdateZoneProperties(ablate::particles::processes::burningModel::TwoZoneBurn::farFieldProp *farfield,
                                                                          ablate::particles::processes::burningModel::TwoZoneBurn::innerZone *innerzone,
                                                                          ablate::particles::processes::burningModel::TwoZoneBurn::outerZone *outerzone) {

    innerzone->T=(Ts+flame->Tad)/2;
    innerzone->Le=1;
    innerzone->MdotF_D_Mdot1=1; //all other species are 0 except fuel

    // First we need to get the mass fractions at the flame
    // I am assuming per simit that the mass fraction of the fuel at the surface is 1




    double Cp1;
    double k1;
    double D1;


    innerzone->Cp


}


void ablate::particles::processes::burningModel::TwoZoneBurn::TwozoneFlame() {


}


void ablate::particles::processes::burningModel::TwoZoneBurn::TwoZoneTransport(double* temperature,double* density,double* Yi[],double* Cp,double* k,double* diffusivity){

    *density=eos->mech->getDensityFromTPY(&temperature, &farField.Pressure,&y[]) const;
    double mu = TransportConst.muo * PetscSqrtReal(*temperature / TransportConst.to) * (*temperature / TransportConst.to) * (TransportConst.to + TransportConst.so) / (*temperature + TransportConst.so);
    *diffusivity = mu / *density / TransportConst.sc;

    *Cp = eos->mech->getMassCpFromTY();
    *k = mu * *Cp / TransportConst.pr;


};


void ablate::particles::processes::burningModel::TwoZoneBurn::SolveTwoZone(double YFsguess,ablate::particles::processes::burningModel::TwoZoneBurn::farFieldProp farfield,
                  ablate::particles::processes::burningModel::TwoZoneBurn::innerZone innerzone,
                  ablate::particles::processes::burningModel::TwoZoneBurn::outerZone outerzone){

    double YFs_old = YFsguess;
    double MWs = 1/(YFsguess/liquidFuel->fuelProperties.MW + (1-YFsguess)/MWair);
    double Pvap = YFsguess * farfield.Pressure * MWs / liquidFuel->fuelProperties.MW;


    //Bound the vapor pressure
    Pvap = std::min(Pvap,farfield.Pressure);
    Pvap = std::max(Pvap, 1E-20);

    //calculate the surface temperature
    double Tsnew = 290;
    liquidFuel->Tvap(&Tsnew,&Pvap);

    double LHS = log((innerzone.MdotF_D_Mdot1 - YFs) / innerzone.MdotF_D_Mdot1) / log(outerzone.MdotOX_D_Mdot2 / (outerzone.MdotOX_D_Mdot2 - farfield.Yox));
    double rf_rs = 1 + LHS * innerzone.gamma / outerzone.gamma;
    double mdot1 = 2 * Pivalue * Dp * innerzone.gamma / (1. - 1. /(rf_rs+1E-20)) * log(innerzone.MdotF_D_Mdot1 / (innerzone.MdotF_D_Mdot1 - YFs));
    mdot1 = std::max(mdot1, 1E-20); // to prevent non-physical solution associated with negative mass flow.

    //TODO make sure Apl is the area
    double Qdot_condl = -liquidFuel->fuelProperties.kl * (Dp* Dp * Pivalue) * (Ts - Tp) / (max(0.5*Dp, 1E-20));
    double Qdot_limiter = abs(qlimfac * mdot1 * liquidFuel->fuelProperties.Hvap) / abs(Qdot_condl + 1E-20);
    Qdot_limiter = std::min(Qdot_limiter,1.e+0);
    Qdot_condl = Qdot_condl * Qdot_limiter;

    double QdotF_D_Mdot = innerzone.Cp * Ts - liquidFuel->fuelProperties.Hvap + Qdot_condl / (mdot1 + 1E-20);
    double QdotO_D_Mdot = QdotF_D_Mdot + innerzone.MdotF_D_Mdot1 * liquidFuel->fuelProperties.Hc;
    double Tf = (outerzone.Cp * farfield.Temperature - QdotO_D_Mdot) * pow(outerzone.MdotOX_D_Mdot2 / (outerzone.MdotOX_D_Mdot2 - farfield.Yox), 1 / outerzone.Le);
    Tf = (QdotO_D_Mdot + Tf) /outerzone.Cp;
    Ts = (QdotF_D_Mdot + (innerzone.Cp * Tf - QdotF_D_Mdot) * exp(-mdot1 * innerzone.Cp * (1. - 1. / rf_rs) / (2. * Pivalue * Dp * innerzone.k))) / innerzone.Cp;

    // Recomputing Yfs using Tf and stepping from flame to surface using YFs-T relation.
    QdotF_D_Mdot = innerzone.Cp * Ts - liquidFuel->fuelProperties.Hvap + Qdot_condl / (mdot1 + 1E-20);
    double YFs_new = innerzone.MdotF_D_Mdot1*(1.- pow((innerzone.Cp * Ts - QdotF_D_Mdot) / (innerzone.Cp * Tf - QdotF_D_Mdot),innerzone.Le));

    YFs_new = std::min(YFs_new, 1.);
    YFs_new = std::max(YFs_new, 0.);


}




#include "registrar.hpp"
REGISTER(ablate::particles::processes::Process, ablate::particles::processes::burningModel::TwoZoneBurn, "Example of an no evaporation/Burning Model",
         ARG(PetscReal, "convectionCoefficient", "The convection Coefficient for the simple heating mode"),
         OPT(PetscReal, "IgnitionTemperature", "The temperature of the far field where it is said to *ignite*"),
         OPT(PetscReal, "burnRate", "The constant burning rate K for the d^2 law (defaults to 0)"),
         OPT(PetscReal, "nuOx", "stoichiometric oxidizer mass coefficient (defaults to 0)"),
         OPT(PetscReal, "Lv", "LatentHeat of Vaporization (defaults to 1e5)"),
         OPT(PetscReal, "Hc", "heat of combustion per unit mass of fuel (defaults to 1e6"),
         ARG(ablate::mathFunctions::FieldFunction, "productsMassFractions", "The product species functions"),
         ARG(std::string, "fuelType", "The type of liquid fuel, curently implmented: 'wax'"),
         OPT(PetscReal, "extinguishmentOxygenMassFraction", "The minimum oxygen mass fraction needed for burning (defaults to 0.1)"),
         ARG(ablate::eos::EOS, "eos", "the eos used to compute various gas properties") );