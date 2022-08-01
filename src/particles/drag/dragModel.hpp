#ifndef ABLATELIBRARY_DRAGMODEL_HPP
#define ABLATELIBRARY_DRAGMODEL_HPP

#include <petsc.h>

namespace ablate::particles::drag {

class DragModel {
   public:
    virtual void ComputeDragForce(const PetscInt dim, const PetscScalar *partVel, const PetscScalar *flowVel, const PetscScalar muF, const PetscScalar rhoF, const PetscScalar partDiam,const PetscScalar partDens,const PetscReal Re,
                                  PetscReal *dragForce) = 0;

    virtual ~DragModel() = default;
};

}  // namespace ablate::particles::drag

#endif  // ABLATELIBRARY_DRAGMODEL_HPP
