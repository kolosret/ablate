#ifndef ABLATELIBRARY_WAXFUEL_H
#define ABLATELIBRARY_WAXFUEL_H

#include "liquidFuel.hpp"
#include <cmath>

namespace ablate::particles::processes::burningModel {

class waxFuel : public ablate::particles::processes::burningModel::LiquidFuel {
   public:
    // Constructor with initializer list
    waxFuel()
        : LiquidFuel("wax", {0.33E6, 41.1E6, 450, 800, 0.12, 2920,1356}, {7.1356, 2276.9, 75.9})  // Initialize base class members
    {}

    //Liquid phase CP from: https://webbook.nist.gov/cgi/cbook.cgi?ID=C544854&Units=SI&Mask=2#Thermo-Condensed


    // Override Tvap method
    void Tvap(double* Tvap, double* Pvap) override {
        double P_hgmm = *Pvap / 133.322;  // converting Pin (Pa) into Hgmm
        double Tresult = (antoineConstants.B / (antoineConstants.A - log10(P_hgmm))) - antoineConstants.C;  // input pressures in Hgmm and returns Temperature in C
        *Tvap = Tresult + 273.15;  // converting the result (Celsius) into Kelvin
    }

    // Override Pvap method
    void Pvap(double* Tvap, double* Pvap) override {
        double T_celsius = *Tvap - 273.15;  // converting Tin (K) into C
        double Presult = (pow(10, antoineConstants.A - antoineConstants.B / (T_celsius + antoineConstants.C)));  // input Temperature in C and returns Pressure in Hgmm
        *Pvap = Presult * 133.322;  // converting the result (Hgmm) into Pa
    }
};

}  // namespace ablate::particles::processes::burningModel

#endif  // ABLATELIBRARY_WAXFUEL_H