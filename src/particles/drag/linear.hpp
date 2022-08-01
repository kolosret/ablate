#ifndef ABLATELIBRARY_LINEAR_HPP
#define ABLATELIBRARY_LINEAR_HPP

#include "dragModel.hpp"

namespace ablate::particles::drag {

class Linear : public DragModel {
   public:
    void ComputeDragForce(const PetscInt dim, const PetscScalar *partVel, const PetscScalar *flowVel, const PetscScalar muF, const PetscScalar rhoF, const PetscScalar partDiam,const PetscScalar partDens, const PetscReal Re, PetscReal *dragSou)
        override;
};

}  // namespace ablate::particles::drag

#endif  // ABLATELIBRARY_LINEAR_HPP
