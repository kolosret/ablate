#ifndef ABLATELIBRARY_CONSTANTK_HPP
#define ABLATELIBRARY_CONSTANTK_HPP

#include "burnModel.hpp"

namespace ablate::particles::processes::burnmodel {



class constantK : public BurnModel {
   public:
    void ComputeBurnRate(const std::shared_ptr<std::vector<double>> *farfield, PetscReal *burnRate, PetscReal *energySource) override;

    explicit constantK(const double burnrate);

    private:
    double Hvap;
    double Pvap;
    double Tboil;
    double kliq;
    double Cl;
    double rhol;
    double K;

};


}

#endif
