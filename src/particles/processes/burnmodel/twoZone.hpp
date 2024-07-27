#ifndef ABLATELIBRARY_TWOZONE_HPP
#define ABLATELIBRARY_TWOZONE_HPP

#include "burnModel.hpp"

namespace ablate::particles::processes::burnmodel {

class TwoZone : public BurnModel {
   public:
    void ComputeBurnRate(const PetscReal *partVel, PetscReal *burnRate, PetscReal *energySource) override;

   private:
    double Hvap;
    double Pvap;
    double Tboil;
    double kliq;
    double Cl;
    double rhol;

};

}

#endif  // ABLATELIBRARY_TWOZONE_HPP
