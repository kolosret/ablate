#include "twoZone.hpp"
#include <math.h>

ablate::particles::processes::burnmodel::TwoZone::TwoZone(const std::string fuel, const std::shared_ptr<eos::EOS> &eos){
    fuelstr=fuel;
}

void ablate::particles::processes::burnmodel::TwoZone::ComputeBurnRate(const std::shared_ptr<std::vector<double>> *farfield, PetscReal *burnRate, PetscReal *energySource){

    /** \brief Two Zone burn model

     */

    *burnRate=1e-6;




}

void ablate::particles::processes::burnmodel::TwoZone::solveYFs(double *YFs){


    // vector to contain important properties of the droplet
    std::vector<double> dropletvector(7);

    //TODO figure out the molar weights
    double MWfarfield = 28.96;
    double MwF = 100; //Molar weight of the fuel
    double Mws = 1/(*YFs/MwF + (1-*YFs)/MWfarfield); // Molar weight at the surface
    double Pg = 101325; //This should be the farfield pressure
    double Pf = Pg * (*YFs) * Mws / MwF; //Partial pressure of the fuel at the surface
    double Pvap=min(Pf,Pg);

    double Tsurf;
    //Vapor pressure curve
    ablate::particles::processes::burnmodel::TwoZone::vaporpressure(&Tsurf,&Pvap);

    //Stuff to define
    // just keep it here:
    double MdotF_D_Mdot1=1;
    double Pi=3.141592653;


    // TODO needs to change
    double MdotOX_D_Mdot2;
    double Yoxinf;
    double Dp=0.001;
    double Apl = Dp*Dp*Pi;
    double Ts;
    double Tp;
    double Tinf;




    double LHS = log((MdotF_D_Mdot1-*YFs) / MdotF_D_Mdot1) / log(MdotOX_D_Mdot2 / (MdotOX_D_Mdot2 - Yoxinf));
    double rf_rs = 1 + LHS * zone1.rhoD / zone2.rhoD;
    double mdot1 = 2 * Pi * Dp * zone1.rhoD / (1. - 1. /(rf_rs + constSmall)) * log(MdotF_D_Mdot1 / (MdotF_D_Mdot1 - *YFs));
    mdot1 = max(mdot1, 1e-22); // to prevent non-physical solution associated with negative mass flow.

    double Qdot_condl = -fuelprops.kliq * Apl * (Ts - Tp) / (max(0.5*Dp, constSmall));
    double Qdot_limiter = abs(fuelprops.qlimfac * mdot1 * fuelprops.Hvap / abs(Qdot_condl + 1e-22));
    Qdot_limiter = min(Qdot_limiter,1.e+0);
    Qdot_condl = Qdot_condl * Qdot_limiter;

    double QdotF_D_Mdot = zone1.Cp * Ts - fuelprops.Hvap + Qdot_condl / (mdot1 + 1e-22);		//no radiation
    double QdotO_D_Mdot = QdotF_D_Mdot + MdotF_D_Mdot1 * fuelprops.Hcomb;
    double Tf = (zone2.Cp * Tinf - QdotO_D_Mdot) * pow(MdotOX_D_Mdot2 / (MdotOX_D_Mdot2 - Yoxinf), 1 / zone2.Le);
    Tf = (QdotO_D_Mdot + Tf) /zone2.Cp;
    Ts = (QdotF_D_Mdot + (zone1.Cp * Tf - QdotF_D_Mdot) * exp(-mdot1 * zone1.Cp * (1. - 1. / rf_rs) / (2. * Pi * Dp * zone1.k))) / zone1.Cp;

    // Recomputing Yfs using Tf and stepping from flame to surface using YFs-T relation.
    QdotF_D_Mdot = zone1.Cp * Ts - fuelprops.Hvap + Qdot_condl / (mdot1 + 1e-22);		//no radiation
    double YFs_new = MdotF_D_Mdot1*(1.-pow((zone1.Cp * Ts - QdotF_D_Mdot) / (zone1.Cp * Tf - QdotF_D_Mdot),zone1.Le));




}


void ablate::particles::processes::burnmodel::TwoZone::vaporpressure(double *Temp, double *Pvap){
     //TODO Change constants for Clausius-Clapeyron

     double Pref=101325;
     double Tref=298;
     double RF=8314/fuelprops.MW;

     *Temp=1/(1/Tref + RF/fuelprops.Hvap*log(Pref/(*Pvap)));
}

void ablate::particles::processes::burnmodel::TwoZone::Fuelproperties::Set(const std::shared_ptr<ablate::parameters::Parameters> &options) {
    if (options) {
        //User provided values
        fuelname = options->Get("fuelname", fuelname);
        Hvap = options->Get("Hvap", Hvap);
        Tboil = options->Get("Tboil", Tboil);
        kliq = options->Get("kliq", kliq);
        rholiq = options->Get("rhol", rholiq);
        qlimfac = options->Get("qlimfac", qlimfac);
    }
}


#include "registrar.hpp"
REGISTER(ablate::particles::processes::burnmodel::BurnModel, ablate::particles::processes::burnmodel::TwoZone, "Twozone Burnmodel",
        ARG(std::string, "fuel", "droplet fuel"),
        ARG(ablate::eos::EOS, "eos", "chemistrymodel eos"));