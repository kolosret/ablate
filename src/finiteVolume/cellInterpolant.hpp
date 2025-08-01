#ifndef ABLATELIBRARY_CELLINTERPOLANT_HPP
#define ABLATELIBRARY_CELLINTERPOLANT_HPP

#include <petsc.h>
#include <vector>
#include "domain/range.hpp"
#include "domain/region.hpp"
#include "domain/subDomain.hpp"
#include "/p/lustre2/kolosret/papi/src/install/include/papi.h"
#include <sys/time.h>
#include <sys/resource.h>

namespace ablate::finiteVolume {
constexpr int numLabel = 7;
class CellInterpolant : private utilities::Loggable<CellInterpolant> {
   public:
    /**
     * Function assumes that the left/right solution and aux variables are discontinuous across the interface
     */
    using DiscontinuousFluxFunction = PetscErrorCode (*)(PetscInt dim, const PetscFVFaceGeom* fg, const PetscInt uOff[], const PetscScalar fieldL[], const PetscScalar fieldR[], const PetscInt aOff[],
                                                         const PetscScalar auxL[], const PetscScalar auxR[], PetscScalar flux[], void* ctx);

    /**
     * Functions that operates on entire cell value.
     */
    using PointFunction = PetscErrorCode (*)(PetscInt dim, PetscReal time, const PetscFVCellGeom* cg, const PetscInt uOff[], const PetscScalar u[], const PetscInt aOff[], const PetscScalar a[],
                                             PetscScalar f[], void* ctx);

    struct DiscontinuousFluxFunctionDescription {
        DiscontinuousFluxFunction function;
        void* context;

        std::vector<PetscInt> updateFields;
        std::vector<PetscInt> inputFields;
        std::vector<PetscInt> auxFields;
    };

    /**
     * struct to describe how to compute RHS finite volume point source terms
     */
    struct PointFunctionDescription {
        PointFunction function;
        void* context;

        std::vector<PetscInt> fields;
        std::vector<PetscInt> inputFields;
        std::vector<PetscInt> auxFields;
    };

    void handle_error(int retval) {
        printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
        exit(1);
    }

    double totalTime = 0.0;
    int callCount = 0;

   private:
    //! use the subDomain to setup the problem
    std::shared_ptr<ablate::domain::SubDomain> subDomain;

    //! store the dmGrad, these are specific to this finite volume solver
    std::vector<DM> gradientCellDms;

    // Maximum value for gradients for the multi-direction flux limiter
    const double maxLimGrad;

    //! Vector to hold all the labels, I dont know the size right now
    std::vector<PetscInt> flowLabelVec;

    /**
     * Function to compute the flux source terms
     */
    void ComputeFluxSourceTerms(DM dm, PetscDS ds, PetscInt totDim, const PetscScalar* xArray, DM dmAux, PetscDS dsAux, PetscInt totDimAux, const PetscScalar* auxArray, DM faceDM,
                                const PetscScalar* faceGeomArray, DM cellDM, const PetscScalar* cellGeomArray, std::vector<DM>& dmGrads, std::vector<const PetscScalar*>& locGradArrays,
                                PetscScalar* locFArray, const std::shared_ptr<domain::Region>& solverRegion, std::vector<CellInterpolant::DiscontinuousFluxFunctionDescription>& rhsFunctions,
                                const ablate::domain::Range& faceRange, const ablate::domain::Range& cellRange);

    /**
     * support call to project to a single face from a side, and also implement face based limiting
     */
    void ProjectToFace(const std::vector<domain::Field>& fields, PetscDS ds, const PetscFVFaceGeom& faceGeom, PetscInt cellId, const PetscFVCellGeom& cellGeom, DM dm, const PetscScalar* xArray,
                       const std::vector<DM>& dmGrads, const std::vector<const PetscScalar*>& gradArrays, PetscScalar* u, PetscScalar* grad, PetscInt neighborCellId,
                       const PetscFVCellGeom& neighborCellGeom, bool projectField = true);

    /**
     * computes the cell gradients
     * @param field
     * @param xLocalVec
     * @param gradLocVec
     * @param dmGrad
     * @param cellGeomVec
     * @param faceGeomVec
     * @param faceRange
     * @param cellRange
     */
    void ComputeFieldGradients(const domain::Field& field, Vec xLocalVec, Vec& gradLocVec, DM& dmGrad, Vec cellGeomVec, Vec faceGeomVec, const ablate::domain::Range& faceRange,
                               const ablate::domain::Range& cellRange);

    /**
     * Helper function to compute the gradient at each cell
     * @param dm
     * @param regionLabel
     * @param regionValue
     * @param fvm
     * @param faceGeometry
     * @param cellGeometry
     * @param dmGrad
     * @return
     */
    static PetscErrorCode ComputeGradientFVM(DM dm, DMLabel regionLabel, PetscInt regionValue, PetscFV fvm, Vec faceGeometry, Vec cellGeometry, DM* dmGrad);

   public:
    /**
     * Create an instance of the cell interpolant for the current solver region
     * @param subDomain
     * @param solverRegion
     * @param faceGeomVec
     * @param cellGeomVec
     */
    CellInterpolant(std::shared_ptr<ablate::domain::SubDomain> subDomain, const std::shared_ptr<domain::Region>& solverRegion, Vec faceGeomVec, Vec cellGeomVec, double maxGradIn,
                    const ablate::domain::Range& faceRange);
    ~CellInterpolant();

    /**
     * Adds in contributions for face based rhs functions
     * @param time
     * @param locXVec
     * @param locFVec
     */
    void ComputeRHS(PetscReal time, Vec locXVec, Vec locAuxVec, Vec locFVec, const std::shared_ptr<domain::Region>& solverRegion,
                    std::vector<CellInterpolant::DiscontinuousFluxFunctionDescription>& rhsFunctions, const ablate::domain::Range& faceRange, const ablate::domain::Range& cellRange, Vec cellGeomVec,
                    Vec faceGeomVec);

    /**
     * Adds in contributions for face based rhs point cell functions
     * @param time
     * @param locXVec
     * @param locFVec
     */
    void ComputeRHS(PetscReal time, Vec locXVec, Vec locAuxVec, Vec locFVec, const std::shared_ptr<domain::Region>& solverRegion, std::vector<CellInterpolant::PointFunctionDescription>& rhsFunctions,
                    const ablate::domain::Range& cellRange, Vec cellGeomVec);
};
}

#endif  // ABLATELIBRARY_CELLINTERPOLANT_HPP