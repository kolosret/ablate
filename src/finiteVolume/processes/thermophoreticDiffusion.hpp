#ifndef ABLATELIBRARY_THERMOPHORETICDIFFUSION_HPP
#define ABLATELIBRARY_THERMOPHORETICDIFFUSION_HPP

#include "eos/transport/transportModel.hpp"
#include "finiteVolume/fluxCalculator/fluxCalculator.hpp"
#include "flowProcess.hpp"

namespace ablate::finiteVolume::processes {

/**
 * Thermophoretic diffusion the transport of ndd (ThermoPheretic) and solid carbon (ThermoPheretic).
 */
class ThermophoreticDiffusion : public FlowProcess {
   private:const std::vector<double> CS_Nasa7TLow = {-3.108720720e-01, 4.403536860e-03, 1.903941180e-06, -6.385469660e-09, 2.989642480e-12, -1.086507940e+02, 1.113829530e+00};
   const std::vector<double> CS_Nasa7THigh = {1.455718290e+00, 1.717022160e-03, -6.975627860e-07, 1.352770320e-10, -9.675906520e-15, -6.951388140e+02, -8.525830330e+00};

    /**
     * Store the equation of state to compute pressure
     */
    const std::shared_ptr<eos::transport::TransportModel> transportModel;

    // store the thermodynamicTemperatureFunction to compute viscosity
    eos::ThermodynamicTemperatureFunction viscosityTemperatureFunction;

   public:
    static double CarbonEnthalpy_R_T(const double& Temp) {
        if ((Temp < 200.)) {
            double t200 = CS_Nasa7TLow.at(5) / 200. + CS_Nasa7TLow.at(0) +
                          200. * (CS_Nasa7TLow.at(1) / 2. + 200. * (CS_Nasa7TLow.at(2) / 3. + 200. * (CS_Nasa7TLow.at(3) / 4. + 200. * CS_Nasa7TLow.at(4) / 5.)));
            double t300 = CS_Nasa7TLow.at(5) / 300. + CS_Nasa7TLow.at(0) +
                          300. * (CS_Nasa7TLow.at(1) / 2. + 300. * (CS_Nasa7TLow.at(2) / 3. + 300. * (CS_Nasa7TLow.at(3) / 4. + 300. * CS_Nasa7TLow.at(4) / 5.)));
            return t200 - (t300 - t200) / 100. * (200. - Temp);
        } else if (Temp <= 1000.)
            return CS_Nasa7TLow.at(5) / Temp + CS_Nasa7TLow.at(0) +
                   Temp * (CS_Nasa7TLow.at(1) / 2. + Temp * (CS_Nasa7TLow.at(2) / 3. + Temp * (CS_Nasa7TLow.at(3) / 4. + Temp * CS_Nasa7TLow.at(4) / 5.)));
        else if (Temp <= 5000.)
            return CS_Nasa7THigh.at(5) / Temp + CS_Nasa7THigh.at(0) +
                   Temp * (CS_Nasa7THigh.at(1) / 2. + Temp * (CS_Nasa7THigh.at(2) / 3. + Temp * (CS_Nasa7THigh.at(3) / 4. + Temp * CS_Nasa7THigh.at(4) / 5.)));
        else {
            double t5000 = CS_Nasa7THigh.at(5) / 5000. + CS_Nasa7THigh.at(0) +
                           5000. * (CS_Nasa7THigh.at(1) / 2. + 5000. * (CS_Nasa7THigh.at(2) / 3. + 5000. * (CS_Nasa7THigh.at(3) / 4. + 5000. * CS_Nasa7THigh.at(4) / 5.)));
            double t4900 = CS_Nasa7THigh.at(5) / 4900. + CS_Nasa7THigh.at(0) +
                           4900. * (CS_Nasa7THigh.at(1) / 2. + 4900. * (CS_Nasa7THigh.at(2) / 3. + 4900. * (CS_Nasa7THigh.at(3) / 4. + 4900. * CS_Nasa7THigh.at(4) / 5.)));
            return t5000 + (t5000 - t4900) / 100. * (Temp - 5000.);
        }
    }

    explicit ThermophoreticDiffusion(std::shared_ptr<eos::transport::TransportModel> transportModel);

    /**
     * public function to link this process with the flow
     * @param flow
     */
    void Setup(ablate::finiteVolume::FiniteVolumeSolver& flow) override;

   private:
    /**
     * This computes the energy transfer for species diffusion flux for rhoE
     * f = "euler"
     * u = {"euler", "densityYi"}
     * a = {"T"}
     * ctx = Viscosity Temperature Function
     * @return
     */
    static PetscErrorCode ThermophoreticDiffusionEnergyFlux(PetscInt dim, const PetscFVFaceGeom* fg, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar field[],
                                                            const PetscScalar grad[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar aux[], const PetscScalar gradAux[],
                                                            PetscScalar flux[], void* ctx);

    /**
     * This computes the species transfer for species diffusion flux
     * f = "densityYi" or "progressYi"
     * u = {"euler", "densityYi"} or {"euler", "progressYi"}
     * a = {"T"}
     * ctx = Viscosity Temperature Function
     * @return
     */
    static PetscErrorCode ThermophoreticDiffusionVariableFlux(PetscInt dim, const PetscFVFaceGeom* fg, const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar field[],
                                                              const PetscScalar grad[], const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar aux[], const PetscScalar gradAux[],
                                                              PetscScalar flux[], void* ctx);
};

}  // namespace ablate::finiteVolume::processes
#endif  // ABLATELIBRARY_SPECIESTRANSPORT_HPP
