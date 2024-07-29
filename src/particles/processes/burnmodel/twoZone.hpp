#ifndef ABLATELIBRARY_TWOZONE_HPP
#define ABLATELIBRARY_TWOZONE_HPP

#include <filesystem>
#include <map>
#include <memory>
#include "burnModel.hpp"
#include "eos/eos.hpp"


namespace ablate::particles::processes::burnmodel {

class TwoZone : public BurnModel {
   public:
    void ComputeBurnRate(const std::shared_ptr<std::vector<double>> *farfield, PetscReal *burnRate, PetscReal *energySource) override;



    explicit TwoZone(std::string fuel, const std::shared_ptr<eos::EOS> &eos);


    struct fuelprops{
        double Hvap=0.5;
        double Tboil=0.35;
        double kliq=0.5;
        double Cl=0.5;
        double rhol=0.8;
    };
   private:
    ablate::particles::processes::burnmodel::TwoZone::fuelprops fuelprops;

    ablate::eos::ChemistryModel::ThermodynamicTemperatureMassFractionFunction specificHeatConstantVolumeFunction;

    eos::ChemistryModel::ThermodynamicTemperatureMassFractionFunction speciesSensibleEnthalpyFunction;

   double constSmall =1e-10;
    double K;
    std::string fuelstr;

    void solveYFs(double *YFs);

    void vaporpressure(double *Temp, double *Pvap);
};


}

#endif  // ABLATELIBRARY_TWOZONE_HPP
