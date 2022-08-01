#include "linear.hpp"

void ablate::particles::drag::Linear::ComputeDragForce(const PetscInt dim, const PetscScalar *partVel, const PetscScalar *flowVel, const PetscScalar muF, const PetscScalar rhoF, const PetscReal partDiam,
                                                       const PetscReal partDens, const PetscReal Re,PetscReal *dragSou) {


//    PetscReal dragForcePrefactor = -3.0 * PETSC_PI * partDiam * muF;

    double Cd = 0.;
    if (Re < 1000) {
        Cd=24.*(1+(pow(Re,2/3)/6))/(Re+0.00001);

    } else{
        Cd=0.424;
    }


    PetscReal dragForce [3];

    for (int n = 0; n < dim; n++) {
        dragForce[n] = 0.5*rhoF*(flowVel[n]-partVel[n])*(flowVel[n]-partVel[n])*Cd*partDiam*partDiam*3.141592/4;
    }

    for (int n = 0; n < dim; n++) {
        dragSou[n]= dragForce[n] /(partDens * PetscPowInt(partDiam,3) *3.1415926/6 + 0.00000001) ;
    }

};

#include "registrar.hpp"
REGISTER_WITHOUT_ARGUMENTS(ablate::particles::drag::DragModel, ablate::particles::drag::Linear, "Computes drag according to Stokes' law.");