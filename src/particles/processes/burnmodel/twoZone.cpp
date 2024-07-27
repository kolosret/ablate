#include "twoZone.hpp"

void ablate::particles::processes::burnmodel::TwoZone::ComputeBurnRate(const PetscReal *partVel, PetscReal *burnRate, PetscReal *energySource){

    /** \brief Two Zone burn model
        ref: Combbustion notes...
     */

    *burnRate=1e-6;

}

#include "registrar.hpp"
REGISTER_WITHOUT_ARGUMENTS(ablate::particles::processes::burnmodel::BurnModel, ablate::particles::processes::burnmodel::TwoZone, "Twozone Burnmodel");