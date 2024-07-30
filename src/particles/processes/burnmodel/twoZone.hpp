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



   private:

    struct Fuelproperties{

        //default values set for nheptane for now...
        //TODO make sure user defines these
        std::string fuelname = "N7C16";
        double Hvap = 316000; //latent heat
        double kliq = 0.0971;
        double rholiq = 684;
        double Tboil = 371;

        //Optional input arguments:
        double qlimfac=0.1;

        //TODO calculate the following values
        double Dp = 1e-3;
        double Hcomb = 4e7;
        double MW = 100;

        void Set(const std::shared_ptr<ablate::parameters::Parameters>&);
    };

    struct Zoneproperties{

        //inner zone between fuel and flame
        // TODO make sure values are set
        double T = 300;
        std::vector<double> Y;
        double rhoD = 0; //latent heat
        double k = 0;
        double Cp = 0;
        double Le = 1;
    };



    ablate::particles::processes::burnmodel::TwoZone::Fuelproperties fuelprops;

    ablate::particles::processes::burnmodel::TwoZone::Zoneproperties zone1;

    ablate::particles::processes::burnmodel::TwoZone::Zoneproperties zone2;

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
