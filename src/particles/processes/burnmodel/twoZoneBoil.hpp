#ifndef ABLATELIBRARY_TWOZONEBOIL_HPP
#define ABLATELIBRARY_TWOZONEBOIL_HPP

#include <filesystem>
#include "burnModel.hpp"
#include "eos/eos.hpp"


namespace ablate::particles::processes::burnmodel {

class TwoZoneBoil: public BurnModel {
   public:
    void ComputeBurnRate(const std::shared_ptr<std::vector<double>> *farfield, PetscReal *burnRate, PetscReal *energySource) override;

    explicit TwoZoneBoil(const std::shared_ptr<ablate::parameters::Parameters> &options);

   private:

    struct Fuelproperties{

        //default values set for nheptane
        // TODO dont set default values or let user know default heptane is used...
        double hc = 4.49260e7; //heat of combustion
        double hfg = 316000; //latent heat
        double Tboil = 371;
        double kfuel = 0.0971;
        double Cpg = 3814;
        double rhol = 684;
        double rox = 11;

        double kair=0.066;

        void Set(const std::shared_ptr<ablate::parameters::Parameters>&);
    };

    ablate::particles::processes::burnmodel::TwoZoneBoil::Fuelproperties Fuelproperties;


};


}

#endif  // ABLATELIBRARY_TWOZONEBOIL_HPP
