#include "constantK.hpp"


ablate::particles::processes::burnmodel::constantK::constantK(const double burnrate){
    K=burnrate;
}

void ablate::particles::processes::burnmodel::constantK::ComputeBurnRate(const std::shared_ptr<std::vector<double>> *farfield, PetscReal *burnRate, PetscReal *energySource){

    /** \brief This method retorns a constant burn rate
     */

    *burnRate=K;

}

#include "registrar.hpp"
REGISTER(ablate::particles::processes::burnmodel::BurnModel, ablate::particles::processes::burnmodel::constantK, "Twozone Burnmodel",
                           ARG(double, "burnrate", "burn rate in Si units"));