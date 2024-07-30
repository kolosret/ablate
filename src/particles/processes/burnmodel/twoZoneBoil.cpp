#include "twoZoneBoil.hpp"
#include <math.h>

ablate::particles::processes::burnmodel::TwoZoneBoil::TwoZoneBoil(const std::shared_ptr<ablate::parameters::Parameters> &options){};

void ablate::particles::processes::burnmodel::TwoZoneBoil::ComputeBurnRate(const std::shared_ptr<std::vector<double>> *farfield, PetscReal *burnRate, PetscReal *energySource){

    /** \brief Two Zone burn model assuming boiling droplets
     * required inputs in SI units:
     * Hcomb, heat of combustion
     * Hvap, latent heat
     * kliq, fuel conductivity
     * Cpg, Cp of the gas at 1000K
     * rhol, liquid density
     * rox, mass of the oxidizer consumed for unit mass of fuel assuming F +rox Ox -> (1+rox)P

     */


    double Tfar=298;
    //Set properties form input file
    double kg=(0.4*Fuelproperties.kliq) + (0.6*Fuelproperties.kair);

    double BoxT=((Fuelproperties.Hcomb/Fuelproperties.rox) + Fuelproperties.Cpg*(Tfar - Fuelproperties.Tboil))/Fuelproperties.Hvap;

    *burnRate = (8*kg*log(1+BoxT))/(Fuelproperties.rhol*Fuelproperties.Cpg);

}


void ablate::particles::processes::burnmodel::TwoZoneBoil::Fuelproperties::Set(const std::shared_ptr<ablate::parameters::Parameters> &options) {
    if (options) {
        Hcomb = options->Get("Hcomb", Hcomb);
        Hvap = options->Get("Hvap", Hvap);
        Tboil = options->Get("Tboil", Tboil);
        kliq = options->Get("kliq", kliq);
        rhol = options->Get("rhol", rhol);
        rox = options->Get("rox", rox);
        Cpg = options->Get("Cpg", Cpg);
    }

    //TODO add check to make sure everything is set
}



#include "registrar.hpp"
REGISTER(ablate::particles::processes::burnmodel::BurnModel, ablate::particles::processes::burnmodel::TwoZoneBoil, "Twozone Burnmodel",
        ARG(ablate::parameters::Parameters, "options","Required arguments:..."));