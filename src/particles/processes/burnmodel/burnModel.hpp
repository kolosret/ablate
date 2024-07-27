#ifndef ABLATELIBRARY_BURNMODEL_HPP
#define ABLATELIBRARY_BURNMODEL_HPP

#include <petsc.h>
#include <filesystem>
#include <map>
#include <memory>
#include "eos/zerork.hpp"

namespace ablate::particles::processes::burnmodel {

class BurnModel {
   public:
    virtual void ComputeBurnRate(const PetscReal *partVel, PetscReal *burnRate, PetscReal *energySource) = 0;

    virtual ~BurnModel() = default;
};

}

#endif