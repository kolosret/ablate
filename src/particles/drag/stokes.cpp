#include "stokes.hpp"
#include "utilities/mathUtilities.hpp"

void ablate::particles::drag::Stokes::ComputeDragForce(const PetscInt dim, const PetscScalar *partVel, const PetscScalar *flowVel, const PetscScalar muF, const PetscScalar rhoF, const PetscReal partDiam,
                                                          const PetscReal partDens,const PetscReal Re,PetscReal *dragForce) {

    PetscReal relVel[3];
    PetscReal dragForcePrefactor = -0.42 * (PETSC_PI / 8.0) * partDiam * partDiam * rhoF;

    for (int n = 0; n < dim; n++) {
        relVel[n] = partVel[n] - flowVel[n];
    }

    dragForcePrefactor *= ablate::utilities::MathUtilities::MagVector(dim, relVel);


    PetscReal corFactor;
    PetscScalar tauP;

    for (int n = 0; n < dim; n++) {
        dragForce[n] = dragForcePrefactor * relVel[n];
    }
    // Correction factor to account for finite Rep on Stokes drag (see Schiller-Naumann drag closure)
    corFactor = 1.0 + 0.15 * PetscPowReal(PetscSqrtReal(Re), 0.687);
    if (Re < 0.1) {
        corFactor = 1.0;  // returns Stokes drag for low speed particles
    }

    for (int n = 0; n < dim; n++) {
        tauP = partDens * PetscSqr(partDiam) / (18.0 * muF);  // particle relaxation time

        dragForce[n]=corFactor *(flowVel[n] - partVel[n]) / tauP ;
//        + g[n] * (1.0 - rhoF / partDens[n])
    }

};

#include "registrar.hpp"
REGISTER_WITHOUT_ARGUMENTS(ablate::particles::drag::DragModel, ablate::particles::drag::Stokes, "Computes drag according to a high Reynolds number drag model for solid spheres.");