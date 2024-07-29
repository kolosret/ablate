#include "twoZoneBoil.hpp"
#include <math.h>

ablate::particles::processes::burnmodel::TwoZoneBoil::TwoZoneBoil(const std::shared_ptr<ablate::parameters::Parameters> &options){};

void ablate::particles::processes::burnmodel::TwoZoneBoil::ComputeBurnRate(const std::shared_ptr<std::vector<double>> *farfield, PetscReal *burnRate, PetscReal *energySource){

    /** \brief Two Zone burn model assuming boiling droplets
     * required inputs in SI units:
     * hc, heat of combustion
     * hfg, latent heat
     * kfuel, fuel conductivity
     * Cpg, Cp of the gas at 1000K
     * rhol, liquid density
     * rox, mass of the oxidizer consumed for unit mass of fuel assuming F +rox Ox -> (1+rox)P

     */


    double Tfar=298;
    //Set properties form input file
    double kg=(0.4*Fuelproperties.kfuel) + (0.6*Fuelproperties.kair);

    double BoxT=((Fuelproperties.hc/Fuelproperties.rox) + Fuelproperties.Cpg*(Tfar - Fuelproperties.Tboil))/Fuelproperties.hfg;

    *burnRate = (8*kg*log(1+BoxT))/(Fuelproperties.rhol*Fuelproperties.Cpg);

}


void ablate::particles::processes::burnmodel::TwoZoneBoil::Fuelproperties::Set(const std::shared_ptr<ablate::parameters::Parameters> &options) {
    if (options) {
        hc = options->Get("hc", hc);
        hfg = options->Get("hfg", hfg);
        Tboil = options->Get("Tboil", Tboil);
        kfuel = options->Get("kfuel", kfuel);
        rhol = options->Get("rhol", rhol);
        rox = options->Get("rox", rox);
        Cpg = options->Get("Cpg", Cpg);
    }

    //TODO add check to make sure everything is set
}



#include "registrar.hpp"
REGISTER(ablate::particles::processes::burnmodel::BurnModel, ablate::particles::processes::burnmodel::TwoZoneBoil, "Twozone Burnmodel",
        ARG(ablate::parameters::Parameters, "options","Required arguments:..."));