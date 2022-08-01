#ifndef ABLATELIBRARY_STOKES_HPP
#define ABLATELIBRARY_STOKES_HPP

#include "dragModel.hpp"

namespace ablate::particles::drag {

class Stokes : public DragModel {
   public:
    void ComputeDragForce(const PetscInt dim, const PetscScalar *partVel, const PetscScalar *flowVel, const PetscScalar muF, const PetscScalar rhoF, const PetscScalar partDiam,const PetscScalar partDens, const PetscReal Re, PetscReal *dragForce)
        override;
};

}  // namespace ablate::particles::drag

#endif  // ABLATELIBRARY_STOKES_HPP
