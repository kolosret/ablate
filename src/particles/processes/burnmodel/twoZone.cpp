#include "twoZone.hpp"
#include <math.h>

ablate::particles::processes::burnmodel::TwoZone::TwoZone(const std::string fuel, const std::shared_ptr<eos::EOS> &eos){
    fuelstr=fuel;
}

void ablate::particles::processes::burnmodel::TwoZone::ComputeBurnRate(const std::shared_ptr<std::vector<double>> *farfield, PetscReal *burnRate, PetscReal *energySource){

    /** \brief Two Zone burn model
        ref: Vanilla two zone model using transfer numbers
     */

    *burnRate=1e-6;




}

void ablate::particles::processes::burnmodel::TwoZone::solveYFs(double *YFs){



    std::vector<double> reactorT(7);

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
    double MdotF_D_Mdot1;
    double MdotOX_D_Mdot2;
    double Yoxinf;
    double rhoD1;
    double rhoD2;
    double Pi=3.14;
    double Dp;
    double k_l;
    double Apl = Dp*Dp*Pi;
    double Ts;
    double Tp;
    double qlimfac=0.1;
    double Cp1;
    double Cp2;
    double Hc;
    double Tinf;
    double k1;
    double Le2;
    double Le1;



    double LHS = log((MdotF_D_Mdot1-*YFs) / MdotF_D_Mdot1) / log(MdotOX_D_Mdot2 / (MdotOX_D_Mdot2 - Yoxinf));
    double rf_rs = 1 + LHS * rhoD1 / rhoD2;
    double mdot1 = 2 * Pi * Dp * rhoD1 / (1. - 1. /(rf_rs + constSmall)) * log(MdotF_D_Mdot1 / (MdotF_D_Mdot1 - *YFs));
    mdot1 = max(mdot1, 1e-22); // to prevent non-physical solution associated with negative mass flow.

    double Qdot_condl = -k_l * Apl * (Ts - Tp) / (max(0.5*Dp, constSmall));
    double Qdot_limiter = abs(qlimfac * mdot1 * fuelprops.Hvap / abs(Qdot_condl + 1e-22));
    Qdot_limiter = min(Qdot_limiter,1.e+0);
    Qdot_condl = Qdot_condl * Qdot_limiter;

    double QdotF_D_Mdot = Cp1 * Ts - fuelprops.Hvap + Qdot_condl / (mdot1 + 1e-22);		//no radiation
    double QdotO_D_Mdot = QdotF_D_Mdot + MdotF_D_Mdot1 * Hc;
    double Tf = (Cp2 * Tinf - QdotO_D_Mdot) * pow(MdotOX_D_Mdot2 / (MdotOX_D_Mdot2 - Yoxinf), 1 / Le2);
    Tf = (QdotO_D_Mdot + Tf) /Cp2;
    Ts = (QdotF_D_Mdot + (Cp1 * Tf - QdotF_D_Mdot) * exp(-mdot1 * Cp1 * (1. - 1. / rf_rs) / (2. * Pi * Dp * k1))) / Cp1;

    // Recomputing Yfs using Tf and stepping from flame to surface using YFs-T relation.
    //			QdotF_D_Mdot = Cp1 * Ts - fuel.Hvap(Ts) + Qdot_condl / (mdot1 + Const.TINY)+(Qdot_rad)/(mdot1 + Const.TINY); // need to recompute since Ts has been updated.
    QdotF_D_Mdot = Cp1 * Ts - fuelprops.Hvap + Qdot_condl / (mdot1 + 1e-22);		//no radiation
    double YFs_new = MdotF_D_Mdot1*(1.-pow((Cp1 * Ts - QdotF_D_Mdot) / (Cp1 * Tf - QdotF_D_Mdot),Le1));




}


void ablate::particles::processes::burnmodel::TwoZone::vaporpressure(double *Temp, double *Pvap){
     //TODO Change constants for Clausius-Clapeyron

     double Pref=101325;
     double Tref=298;
     double RF=8314;
     double hfg=316000;

     *Temp=1/(1/Tref + RF/hfg*log(Pref/(*Pvap)));
}




#include "registrar.hpp"
REGISTER(ablate::particles::processes::burnmodel::BurnModel, ablate::particles::processes::burnmodel::TwoZone, "Twozone Burnmodel",
        ARG(std::string, "fuel", "droplet fuel"),
        ARG(ablate::eos::EOS, "eos", "chemistrymodel eos"));