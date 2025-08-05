#include "cellInterpolant.hpp"
#include <petsc/private/dmpleximpl.h>
#include <utility>

ablate::finiteVolume::CellInterpolant::CellInterpolant(std::shared_ptr<ablate::domain::SubDomain> subDomainIn, const std::shared_ptr<domain::Region>& solverRegion, Vec faceGeomVec, Vec cellGeomVec,
                                                       double maxGradIn,const ablate::domain::Range& faceRange)
    : subDomain(std::move(std::move(subDomainIn))), maxLimGrad(maxGradIn),flowLabelVec(numLabel * (faceRange.end - faceRange.start), 0) {
//    StartEvent("FiniteVolumeSolver::CellInterpolant::CompluteRHS::constructor");
    auto getGradientDm = [this, solverRegion, faceGeomVec, cellGeomVec](const domain::Field& fieldInfo, std::vector<DM>& gradDMs) {
        auto petscField = subDomain->GetPetscFieldObject(fieldInfo);
        auto petscFieldFV = (PetscFV)petscField;

        PetscBool computeGradients;
        PetscFVGetComputeGradients(petscFieldFV, &computeGradients) >> utilities::PetscUtilities::checkError;

        if (computeGradients) {
            DM dmGradInt;

            DMLabel regionLabel = nullptr;
            PetscInt regionValue = PETSC_DECIDE;
            domain::Region::GetLabel(solverRegion, subDomain->GetDM(), regionLabel, regionValue);

            ComputeGradientFVM(subDomain->GetFieldDM(fieldInfo), regionLabel, regionValue, petscFieldFV, faceGeomVec, cellGeomVec, &dmGradInt) >> utilities::PetscUtilities::checkError;
            gradDMs.push_back(dmGradInt);
        } else {
            gradDMs.push_back(nullptr);
        }
    };


//    faceCellsVec.resize(2 * (faceRange.end - faceRange.start), 0);


    DMLabel regionLabel = nullptr;
    PetscInt regionValue = PETSC_DECIDE;
    domain::Region::GetLabel(solverRegion, subDomain->GetDM(), regionLabel, regionValue);

    auto dm = subDomain->GetDM();

    DMLabel ghostLabel;
    DMGetLabel(dm, "ghost", &ghostLabel) >> utilities::PetscUtilities::checkError;

    DM faceDM, cellDM;
    VecGetDM(faceGeomVec, &faceDM) >> utilities::PetscUtilities::checkError;
    VecGetDM(cellGeomVec, &cellDM) >> utilities::PetscUtilities::checkError;

    const PetscScalar* cellGeomArray = nullptr;
    const PetscScalar* faceGeomArray = nullptr;
    VecGetArrayRead(cellGeomVec, &cellGeomArray) >> utilities::PetscUtilities::checkError;
    VecGetArrayRead(faceGeomVec, &faceGeomArray) >> utilities::PetscUtilities::checkError;



    for (PetscInt f = faceRange.start; f < faceRange.end; ++f) {
        const PetscInt face = faceRange.points ? faceRange.points[f] : f;

        //! The indicies are in order per face:
        // ghost, nsupp, nchild, leftFlowLabelValue, rightFlowLabelValue


        // make sure that this is a valid face
        PetscInt ghost, nsupp, nchild;
        DMLabelGetValue(ghostLabel, face, &ghost) >> utilities::PetscUtilities::checkError;
        DMPlexGetSupportSize(dm, face, &nsupp) >> utilities::PetscUtilities::checkError;
        DMPlexGetTreeChildren(dm, face, &nchild, nullptr) >> utilities::PetscUtilities::checkError;

        flowLabelVec[numLabel*f]=ghost;
        flowLabelVec[numLabel*f+1]=nsupp;
        flowLabelVec[numLabel*f+2]=nchild;


        // Get the face geometry
        const PetscInt* faceCells;
        PetscFVFaceGeom* fg;
        PetscFVCellGeom *cgL, *cgR;
        DMPlexPointLocalRead(faceDM, face, faceGeomArray, &fg) >> utilities::PetscUtilities::checkError;
        DMPlexGetSupport(dm, face, &faceCells) >> utilities::PetscUtilities::checkError;
        DMPlexPointLocalRead(cellDM, faceCells[0], cellGeomArray, &cgL) >> utilities::PetscUtilities::checkError;
        DMPlexPointLocalRead(cellDM, faceCells[1], cellGeomArray, &cgR) >> utilities::PetscUtilities::checkError;

        PetscInt leftFlowLabelValue = regionValue;
        PetscInt rightFlowLabelValue = regionValue;
        //        start = MPI_Wtime();
        if (regionLabel) {
            DMLabelGetValue(regionLabel, faceCells[0], &leftFlowLabelValue);
            DMLabelGetValue(regionLabel, faceCells[1], &rightFlowLabelValue);
        }
        flowLabelVec[numLabel*f+3]=leftFlowLabelValue;
        flowLabelVec[numLabel*f+4]=rightFlowLabelValue;

        //Check if the cells are ghost cells
        if (ghostLabel) {
            DMLabelGetValue(ghostLabel, faceCells[0], &leftFlowLabelValue);
            DMLabelGetValue(ghostLabel, faceCells[1], &rightFlowLabelValue);
        }
        flowLabelVec[numLabel*f+5]=leftFlowLabelValue;
        flowLabelVec[numLabel*f+6]=rightFlowLabelValue;

    }



    // Compute the gradient dm for each field that supports it
    for (const auto& fieldInfo : subDomain->GetFields()) {
        getGradientDm(fieldInfo, gradientCellDms);
    }
//    EndEvent();
}

ablate::finiteVolume::CellInterpolant::~CellInterpolant() {
    for (auto& dm : gradientCellDms) {
        if (dm) {
            DMDestroy(&dm) >> utilities::PetscUtilities::checkError;
        }
    }
}

void ablate::finiteVolume::CellInterpolant::ComputeRHS(PetscReal time, Vec locXVec, Vec locAuxVec, Vec locFVec, const std::shared_ptr<domain::Region>& solverRegion,
                                                       std::vector<CellInterpolant::DiscontinuousFluxFunctionDescription>& rhsFunctions, const ablate::domain::Range& faceRange,
                                                       const ablate::domain::Range& cellRange, Vec cellGeomVec, Vec faceGeomVec) {
    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::getDMsize");
    auto dm = subDomain->GetDM();
    auto dmAux = subDomain->GetAuxDM();


    /* 1: Get sizes from dm and dmAux */
    PetscSection section = nullptr;
    DMGetLocalSection(dm, &section) >> utilities::PetscUtilities::checkError;

    // Get the ds from he subDomain and required info
    auto ds = subDomain->GetDiscreteSystem();
    PetscInt nf, totDim;
    PetscDSGetNumFields(ds, &nf) >> utilities::PetscUtilities::checkError;
    PetscDSGetTotalDimension(ds, &totDim) >> utilities::PetscUtilities::checkError;


    // Check to see if the dm has an auxVec/auxDM associated with it.  If it does, extract it
    PetscDS dsAux = subDomain->GetAuxDiscreteSystem();
    PetscInt naf = 0, totDimAux = 0;
    if (locAuxVec) {
        PetscDSGetTotalDimension(dsAux, &totDimAux) >> utilities::PetscUtilities::checkError;
        PetscDSGetNumFields(dsAux, &naf) >> utilities::PetscUtilities::checkError;
    }
    EndEvent();

    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::getgeoData");
    /* 2: Get geometric data */
    // We can use a single call for the geometry data because it does not depend on the fv object
    const PetscScalar* cellGeomArray = nullptr;
    const PetscScalar* faceGeomArray = nullptr;
    VecGetArrayRead(cellGeomVec, &cellGeomArray) >> utilities::PetscUtilities::checkError;
    VecGetArrayRead(faceGeomVec, &faceGeomArray) >> utilities::PetscUtilities::checkError;
    DM faceDM, cellDM;
    VecGetDM(faceGeomVec, &faceDM) >> utilities::PetscUtilities::checkError;
    VecGetDM(cellGeomVec, &cellDM) >> utilities::PetscUtilities::checkError;

    // Get raw access to the computed values
    const PetscScalar *xArray, *auxArray = nullptr;
    VecGetArrayRead(locXVec, &xArray) >> utilities::PetscUtilities::checkError;
    if (locAuxVec) {
        VecGetArrayRead(locAuxVec, &auxArray) >> utilities::PetscUtilities::checkError;
    }

    // get raw access to the locF
    PetscScalar* locFArray;
    VecGetArray(locFVec, &locFArray) >> utilities::PetscUtilities::checkError;

    // there must be a separate gradient vector/dm for field because they can be different sizes
    std::vector<Vec> locGradVecs(nf, nullptr);
    EndEvent();


    /* Reconstruct and limit cell gradients */
    // for each field compute the gradient in the localGrads vector
    for (const auto& field : subDomain->GetFields()) {
        // This is logged in the function
        ComputeFieldGradients(field, locXVec, locGradVecs[field.subId], gradientCellDms[field.subId], cellGeomVec, faceGeomVec, faceRange, cellRange);
    }

//    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::setlocalvec");
    std::vector<const PetscScalar*> locGradArrays(nf, nullptr);
    for (const auto& field : subDomain->GetFields()) {
        if (locGradVecs[field.subId]) {
            VecGetArrayRead(locGradVecs[field.subId], &locGradArrays[field.subId]) >> utilities::PetscUtilities::checkError;
        }
    }
//    EndEvent();
    if (time <= 1E-6) std::cout << "The mesh has " << faceRange.end << "  faces \n" ;

    //This is measured inside the function
    ComputeFluxSourceTerms(dm,
                           ds,
                           totDim,
                           xArray,
                           dmAux,
                           dsAux,
                           totDimAux,
                           auxArray,
                           faceDM,
                           faceGeomArray,
                           cellDM,
                           cellGeomArray,
                           gradientCellDms,
                           locGradArrays,
                           locFArray,
                           solverRegion,
                           rhsFunctions,
                           faceRange,
                           cellRange);
    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::cleanup");
    // clean up cell grads
    for (const auto& field : subDomain->GetFields()) {
        if (locGradVecs[field.subId]) {
            VecRestoreArrayRead(locGradVecs[field.subId], &locGradArrays[field.subId]) >> utilities::PetscUtilities::checkError;
            DMRestoreLocalVector(gradientCellDms[field.subId], &locGradVecs[field.subId]) >> utilities::PetscUtilities::checkError;
        }
    }

    // cleanup (restore access to locGradVecs, locAuxGradVecs with DMRestoreLocalVector)
    VecRestoreArrayRead(locXVec, &xArray) >> utilities::PetscUtilities::checkError;
    if (locAuxVec) {
        VecRestoreArrayRead(locAuxVec, &auxArray) >> utilities::PetscUtilities::checkError;
    }

    VecRestoreArray(locFVec, &locFArray) >> utilities::PetscUtilities::checkError;
    VecRestoreArrayRead(faceGeomVec, (const PetscScalar**)&faceGeomArray) >> utilities::PetscUtilities::checkError;
    VecRestoreArrayRead(cellGeomVec, (const PetscScalar**)&cellGeomArray) >> utilities::PetscUtilities::checkError;
    EndEvent();
}

void ablate::finiteVolume::CellInterpolant::ComputeRHS(PetscReal time, Vec locXVec, Vec locAuxVec, Vec locFVec, const std::shared_ptr<domain::Region>& solverRegion,
                                                       std::vector<CellInterpolant::PointFunctionDescription>& rhsFunctions, const ablate::domain::Range& cellRange, Vec cellGeomVec) {
    auto dm = subDomain->GetDM();
    auto dmAux = subDomain->GetAuxDM();

    /* 1: Get sizes from dm and dmAux */
    PetscSection section = nullptr;
    DMGetLocalSection(dm, &section) >> utilities::PetscUtilities::checkError;

    // Get the ds from he subDomain and required info
    auto ds = subDomain->GetDiscreteSystem();
    PetscInt nf, totDim;
    PetscDSGetNumFields(ds, &nf) >> utilities::PetscUtilities::checkError;
    PetscDSGetTotalDimension(ds, &totDim) >> utilities::PetscUtilities::checkError;

    // Check to see if the dm has an auxVec/auxDM associated with it.  If it does, extract it
    PetscDS dsAux = subDomain->GetAuxDiscreteSystem();
    PetscInt naf = 0, totDimAux = 0;
    if (locAuxVec) {
        PetscDSGetTotalDimension(dsAux, &totDimAux) >> utilities::PetscUtilities::checkError;
        PetscDSGetNumFields(dsAux, &naf) >> utilities::PetscUtilities::checkError;
    }

    // We can use a single call for the geometry data because it does not depend on the fv object
    const PetscScalar* cellGeomArray = nullptr;
    VecGetArrayRead(cellGeomVec, &cellGeomArray) >> utilities::PetscUtilities::checkError;
    DM cellDM;
    VecGetDM(cellGeomVec, &cellDM) >> utilities::PetscUtilities::checkError;

    // Get raw access to the computed values
    const PetscScalar *xArray, *auxArray = nullptr;
    VecGetArrayRead(locXVec, &xArray) >> utilities::PetscUtilities::checkError;
    if (locAuxVec) {
        VecGetArrayRead(locAuxVec, &auxArray) >> utilities::PetscUtilities::checkError;
    }

    // get raw access to the locF
    PetscScalar* locFArray;
    VecGetArray(locFVec, &locFArray) >> utilities::PetscUtilities::checkError;

    // Compute the source terms from flux across the interface for cell based gradient functions
    // Precompute the offsets to pass into the rhsFluxFunctionDescriptions
    std::vector<std::vector<PetscInt>> fluxComponentSize(rhsFunctions.size());
    std::vector<std::vector<PetscInt>> fluxComponentOffset(rhsFunctions.size());
    std::vector<std::vector<PetscInt>> uOff(rhsFunctions.size());
    std::vector<std::vector<PetscInt>> aOff(rhsFunctions.size());

    // Get the full set of offsets from the ds
    PetscInt* uOffTotal;
    PetscDSGetComponentOffsets(ds, &uOffTotal) >> utilities::PetscUtilities::checkError;

    for (std::size_t fun = 0; fun < rhsFunctions.size(); fun++) {
        for (std::size_t f = 0; f < rhsFunctions[fun].fields.size(); f++) {
            const auto& field = subDomain->GetField(rhsFunctions[fun].fields[f]);

            PetscInt fieldSize, fieldOffset;
            PetscDSGetFieldSize(ds, field.subId, &fieldSize) >> utilities::PetscUtilities::checkError;
            PetscDSGetFieldOffset(ds, field.subId, &fieldOffset) >> utilities::PetscUtilities::checkError;
            fluxComponentSize[fun].push_back(fieldSize);
            fluxComponentOffset[fun].push_back(fieldOffset);
        }

        for (std::size_t f = 0; f < rhsFunctions[fun].inputFields.size(); f++) {
            uOff[fun].push_back(uOffTotal[rhsFunctions[fun].inputFields[f]]);
        }
    }

    if (dsAux) {
        PetscInt* auxOffTotal;
        PetscDSGetComponentOffsets(dsAux, &auxOffTotal) >> utilities::PetscUtilities::checkError;
        for (std::size_t fun = 0; fun < rhsFunctions.size(); fun++) {
            for (std::size_t f = 0; f < rhsFunctions[fun].auxFields.size(); f++) {
                aOff[fun].push_back(auxOffTotal[rhsFunctions[fun].auxFields[f]]);
            }
        }
    }

    // check to see if there is a ghost label
    DMLabel ghostLabel;
    DMGetLabel(dm, "ghost", &ghostLabel) >> utilities::PetscUtilities::checkError;

    PetscInt dim = subDomain->GetDimensions();

    // Size up a scratch variable
    PetscScalar fScratch[totDim];

    // March over each cell
    for (PetscInt c = cellRange.start; c < cellRange.end; ++c) {
        // if there is a cell array, use it, otherwise it is just c
        const PetscInt cell = cellRange.points ? cellRange.points[c] : c;

        // make sure that this is not a ghost cell
        if (ghostLabel) {
            PetscInt ghostVal;

            DMLabelGetValue(ghostLabel, cell, &ghostVal) >> utilities::PetscUtilities::checkError;
            if (ghostVal > 0) continue;
        }

        // extract the point locations for this cell
        const PetscFVCellGeom* cg;
        const PetscScalar* u;
        PetscScalar* rhs;
        DMPlexPointLocalRead(cellDM, cell, cellGeomArray, &cg) >> utilities::PetscUtilities::checkError;
        DMPlexPointLocalRead(dm, cell, xArray, &u) >> utilities::PetscUtilities::checkError;
        DMPlexPointLocalRef(dm, cell, locFArray, &rhs) >> utilities::PetscUtilities::checkError;

        // if there is an aux field, get it
        const PetscScalar* a = nullptr;
        if (auxArray) {
            DMPlexPointLocalRead(dmAux, cell, auxArray, &a) >> utilities::PetscUtilities::checkError;
        }

        // March over each functionDescriptions
        for (std::size_t fun = 0; fun < rhsFunctions.size(); fun++) {
            rhsFunctions[fun].function(dim, time, cg, uOff[fun].data(), u, aOff[fun].data(), a, fScratch, rhsFunctions[fun].context) >> utilities::PetscUtilities::checkError;

            // copy over each result flux field
            PetscInt r = 0;
            for (std::size_t ff = 0; ff < rhsFunctions[fun].fields.size(); ff++) {
                for (PetscInt d = 0; d < fluxComponentSize[fun][ff]; ++d) {
                    rhs[fluxComponentOffset[fun][ff] + d] += fScratch[r++];
                }
            }
        }
    }

    // cleanup (restore access to locGradVecs, locAuxGradVecs with DMRestoreLocalVector)
    VecRestoreArrayRead(locXVec, &xArray) >> utilities::PetscUtilities::checkError;
    if (locAuxVec) {
        VecRestoreArrayRead(locAuxVec, &auxArray) >> utilities::PetscUtilities::checkError;
    }

    VecRestoreArray(locFVec, &locFArray) >> utilities::PetscUtilities::checkError;
    VecRestoreArrayRead(cellGeomVec, &cellGeomArray) >> utilities::PetscUtilities::checkError;
}

/**
 * This is a duplication of PETSC that we don't have access to
 */
static PetscErrorCode DMPlexApplyLimiter_Internal(DM dm, DM dmCell, PetscLimiter lim, PetscInt dim, PetscInt dof, PetscInt cell, PetscInt field, PetscInt face, PetscInt fStart, PetscInt fEnd,
                                                  PetscReal* cellPhi, const PetscScalar* x, const PetscScalar* cellgeom, const PetscFVCellGeom* cg, const PetscScalar* cx, const PetscScalar* cgrad) {
    const PetscInt* children;
    PetscInt numChildren;

    PetscFunctionBegin;
    PetscCall(DMPlexGetTreeChildren(dm, face, &numChildren, &children));
    if (numChildren) {  // if the tree contains children
        PetscInt c;

        for (c = 0; c < numChildren; c++) {
            PetscInt childFace = children[c];

            if (childFace >= fStart && childFace < fEnd) {
                PetscCall(DMPlexApplyLimiter_Internal(dm, dmCell, lim, dim, dof, cell, field, childFace, fStart, fEnd, cellPhi, x, cellgeom, cg, cx, cgrad));
            }
        }
    } else {                     // if the tree doesn't contain children
        PetscScalar* ncx;        // neighbor cell centered values
        PetscFVCellGeom* ncg;    // neighbor cell geometry
        const PetscInt* fcells;  // cells attached to this face
        PetscInt ncell, d;       // neighbor cell and for loop index
        PetscReal v[3];          // centr

        PetscCall(DMPlexGetSupport(dm, face, &fcells));
        ncell = cell == fcells[0] ? fcells[1] : fcells[0];  // figure out which cell is the neighbor and not this cell
        // Read in the neighbor cell information
        if (field >= 0) {
            PetscCall(DMPlexPointLocalFieldRead(dm, ncell, field, x, &ncx));
        } else {
            PetscCall(DMPlexPointLocalRead(dm, ncell, x, &ncx));
        }
        PetscCall(DMPlexPointLocalRead(dmCell, ncell, cellgeom, &ncg));
        // Calculate the distance between the neighbor cell and this cell
        DMPlex_WaxpyD_Internal(dim, -1, cg->centroid, ncg->centroid, v);  // v_i = NeighborCentroid_i - ThisCentroid_i = dx_i
        for (d = 0; d < dof; ++d) {
            /* We use the symmetric slope limited form of Berger, Aftosmis, and Murman 2005 */
            PetscReal denom = DMPlex_DotD_Internal(dim, &cgrad[d * dim], v);    // denominator = \grad u \cdot dx
            PetscReal phi, flim = 0.5 * PetscRealPart(ncx[d] - cx[d]) / denom;  // f = 1/2 * (u_i+1-u_i)/(\Delta u from cell gradient dot with \delta x)
            // What the above means is that if any cell face has no change, but there was ample enough change close to it such that
            //  the cell gradient dot with delta x is not 0, there is no limiting... i.e f = 0 and all limiters = 0
            PetscCall(PetscLimiterLimit(lim, flim, &phi));
            cellPhi[d] = PetscMin(cellPhi[d], phi);
        }
    }
    PetscFunctionReturn(0);
}

void ablate::finiteVolume::CellInterpolant::ComputeFieldGradients(const domain::Field& field, Vec xLocalVec, Vec& gradLocVec, DM& dmGrad, Vec cellGeomVec, Vec faceGeomVec,
                                                                  const ablate::domain::Range& faceRange, const ablate::domain::Range& cellRange) {
    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::ComputeFieldGradients");
    // get the FVM petsc field associated with this field
    auto fvm = (PetscFV)subDomain->GetPetscFieldObject(field); // 1 read
    auto dm = subDomain->GetFieldDM(field); // 1 read

    // Get the dm for this grad field
    // If there is no grad, return
    if (!dmGrad) { // 1 write
        return;
    }

    // Create a gradLocVec
    DMGetLocalVector(dmGrad, &gradLocVec) >> utilities::PetscUtilities::checkError;

    // Get the correct sized vec (gradient for this field)
    Vec gradGlobVec; // 1 write
    DMGetGlobalVector(dmGrad, &gradGlobVec) >> utilities::PetscUtilities::checkError; // N read
    VecZeroEntries(gradGlobVec) >> utilities::PetscUtilities::checkError; // vector length writes

//    PetscInt size;
//    VecGetSize(gradGlobVec, &size)>> utilities::PetscUtilities::checkError;
//    std::cout << "Global size of gradGlobVec: " << size << std::endl;

    // check to see if there is a ghost label
    DMLabel ghostLabel;
    DMGetLabel(dm, "ghost", &ghostLabel) >> utilities::PetscUtilities::checkError; // 1 read, 1 write

    // Get the face geometry
    DM dmFace;
    const PetscScalar* faceGeometryArray;
    VecGetDM(faceGeomVec, &dmFace) >> utilities::PetscUtilities::checkError;
    VecGetArrayRead(faceGeomVec, &faceGeometryArray);

    // extract the local x array
    const PetscScalar* xLocalArray;
    VecGetArrayRead(xLocalVec, &xLocalArray);

    // extract the global grad array
    PetscScalar* gradGlobArray;
    VecGetArray(gradGlobVec, &gradGlobArray); //This returns just the pointer, no copies.

    // Get the dof and dim
    PetscInt dim = subDomain->GetDimensions();
    PetscInt dof = field.numberComponents;

    for (PetscInt f = faceRange.start; f < faceRange.end; ++f) {
        PetscInt face = faceRange.points ? faceRange.points[f] : f;

        // make sure that this is a face we should use
        PetscBool boundary;
        PetscInt ghost = -1;
        if (ghostLabel) {
            DMLabelGetValue(ghostLabel, face, &ghost); //1 read
        }
        DMIsBoundaryPoint(dm, face, &boundary);
        PetscInt numChildren;
        DMPlexGetTreeChildren(dm, face, &numChildren, nullptr);
        if (ghost >= 0 || boundary || numChildren) continue;

        // Do a sanity check on the number of cells connected to this face
        PetscInt numCells;
        DMPlexGetSupportSize(dm, face, &numCells);
        if (numCells != 2) {
            throw std::runtime_error("face " + std::to_string(face) + " has " + std::to_string(numCells) + " support points (cells): expected 2");
        }

        // add in the contributions from this face
        const PetscInt* cells;
        PetscFVFaceGeom* fg;
        PetscScalar* cx[2];
        PetscScalar* cgrad[2];

        DMPlexGetSupport(dm, face, &cells);
        DMPlexPointLocalRead(dmFace, face, faceGeometryArray, &fg);
        for (PetscInt c = 0; c < 2; ++c) {
            DMPlexPointLocalFieldRead(dm, cells[c], field.id, xLocalArray, &cx[c]) >> utilities::PetscUtilities::checkError;
            DMPlexPointGlobalRef(dmGrad, cells[c], gradGlobArray, &cgrad[c]) >> utilities::PetscUtilities::checkError;
        }
        for (PetscInt pd = 0; pd < dof; ++pd) {
            PetscScalar delta = cx[1][pd] - cx[0][pd];

            for (PetscInt d = 0; d < dim; ++d) {
                if (cgrad[0]) cgrad[0][pd * dim + d] += fg->grad[0][d] * delta;
                if (cgrad[1]) cgrad[1][pd * dim + d] -= fg->grad[1][d] * delta;
//                std::cout << "The index is: " << pd * dim + d << std::endl;

            }
        }
    }
    EndEvent();
    // Check for a limiter the limiter
    PetscLimiter lim;
    PetscFVGetLimiter(fvm, &lim) >> utilities::PetscUtilities::checkError;
    if (false) {
        /* Limit interior gradients (using cell-based loop because it generalizes better to vector limiters) */
        // Get the cell geometry
        DM dmCell;
        const PetscScalar* cellGeometryArray;
        VecGetDM(cellGeomVec, &dmCell) >> utilities::PetscUtilities::checkError;
        VecGetArrayRead(cellGeomVec, &cellGeometryArray);

        // create a temp work array
        PetscReal* cellPhi;                                                                     // Actual Limiter
        DMGetWorkArray(dm, dof, MPIU_REAL, &cellPhi) >> utilities::PetscUtilities::checkError;  // Size it up to be dof length of reals

        for (PetscInt c = cellRange.start; c < cellRange.end; ++c) {
            PetscInt cell = cellRange.points ? cellRange.points[c] : c;

            const PetscInt* cellFaces;  // faces belonging to the cell
            PetscScalar* cx;            // cell solution/field values
            PetscFVCellGeom* cg;        // cell geometry
            PetscScalar* cgrad;         // cell gradient
            PetscInt coneSize;          // cell connectivity

            DMPlexGetConeSize(dm, cell, &coneSize) >> utilities::PetscUtilities::checkError;
            DMPlexGetCone(dm, cell, &cellFaces) >> utilities::PetscUtilities::checkError;
            DMPlexPointLocalFieldRead(dm, cell, field.id, xLocalArray, &cx) >> utilities::PetscUtilities::checkError;
            DMPlexPointLocalRead(dmCell, cell, cellGeometryArray, &cg) >> utilities::PetscUtilities::checkError;
            DMPlexPointGlobalRef(dmGrad, cell, gradGlobArray, &cgrad) >> utilities::PetscUtilities::checkError;

            if (!cgrad) {
                /* Unowned overlap cell, we do not compute */
                continue;
            }
            /* Limiter will be minimum value over all neighbors */
            for (PetscInt d = 0; d < dof; ++d) {
                cellPhi[d] = PETSC_MAX_REAL;
            }
            for (PetscInt f = 0; f < coneSize; ++f) {
                DMPlexApplyLimiter_Internal(dm, dmCell, lim, dim, dof, cell, field.id, cellFaces[f], faceRange.start, faceRange.end, cellPhi, xLocalArray, cellGeometryArray, cg, cx, cgrad) >>
                    utilities::PetscUtilities::checkError;
            }

            /* Apply limiter to gradient */
            PetscBool cancel = PETSC_FALSE;
            for (PetscInt pd = 0; pd < dof; ++pd) {
                if (cellPhi[pd] == 0) {
                    for (PetscInt d = 0; d < dim; d++) {
                        // Due to directional limiting being difficult in unstructured grids,
                        // a strict gradient limiter is introduced here to revert back to cell centered
                        // reconstructions if a component limiter is 0 even though there is a strong
                        // component gradient in a direction (Usually happens at strong shocks seen in
                        // high pressured rocket simulations)
                        if (PetscAbsReal(cgrad[pd * dim + d]) > maxLimGrad) cancel = PETSC_TRUE;
                    }
                }
            }

            for (PetscInt pd = 0; pd < dof; ++pd) {
                /* Scalar limiter applied to each component separately */
                for (PetscInt d = 0; d < dim; ++d) {
                    if (cancel)
                        cgrad[pd * dim + d] *= 0;
                    else
                        cgrad[pd * dim + d] *= cellPhi[pd];
                }
            }
        }

        // clean up the limiter work
        DMRestoreWorkArray(dm, dof, MPIU_REAL, &cellPhi) >> utilities::PetscUtilities::checkError;
        VecRestoreArrayRead(cellGeomVec, &cellGeometryArray);
    }

    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::GradientComm");
    // Communicate gradient values
    VecRestoreArray(gradGlobVec, &gradGlobArray) >> utilities::PetscUtilities::checkError;
    DMGlobalToLocalBegin(dmGrad, gradGlobVec, INSERT_VALUES, gradLocVec) >> utilities::PetscUtilities::checkError;
    DMGlobalToLocalEnd(dmGrad, gradGlobVec, INSERT_VALUES, gradLocVec) >> utilities::PetscUtilities::checkError;

    // cleanup
    VecRestoreArrayRead(xLocalVec, &xLocalArray) >> utilities::PetscUtilities::checkError;
    VecRestoreArrayRead(faceGeomVec, &faceGeometryArray) >> utilities::PetscUtilities::checkError;
    DMRestoreGlobalVector(dmGrad, &gradGlobVec) >> utilities::PetscUtilities::checkError;
    EndEvent();
}

void ablate::finiteVolume::CellInterpolant::ComputeFluxSourceTerms(DM dm, PetscDS ds, PetscInt totDim, const PetscScalar* xArray, DM dmAux, PetscDS dsAux, PetscInt totDimAux,
                                                                   const PetscScalar* auxArray, DM faceDM, const PetscScalar* faceGeomArray, DM cellDM, const PetscScalar* cellGeomArray,
                                                                   std::vector<DM>& dmGrads, std::vector<const PetscScalar*>& locGradArrays, PetscScalar* locFArray,
                                                                   const std::shared_ptr<domain::Region>& solverRegion,
                                                                   std::vector<CellInterpolant::DiscontinuousFluxFunctionDescription>& rhsFunctions, const ablate::domain::Range& faceRange,
                                                                   const ablate::domain::Range& cellRange) {
//    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::ComputeSource::Setup");
    PetscInt dim = subDomain->GetDimensions();


    // Size up the work arrays (uL, uR, gradL, gradR, auxL, auxR, gradAuxL, gradAuxR), these are only sized for one face at a time
    PetscScalar* flux;
    DMGetWorkArray(dm, totDim, MPIU_SCALAR, &flux) >> utilities::PetscUtilities::checkError;

    PetscScalar *uL, *uR;
    DMGetWorkArray(dm, totDim, MPIU_SCALAR, &uL) >> utilities::PetscUtilities::checkError;
    DMGetWorkArray(dm, totDim, MPIU_SCALAR, &uR) >> utilities::PetscUtilities::checkError;

    PetscScalar *gradL, *gradR;
    DMGetWorkArray(dm, dim * totDim, MPIU_SCALAR, &gradL) >> utilities::PetscUtilities::checkError;
    DMGetWorkArray(dm, dim * totDim, MPIU_SCALAR, &gradR) >> utilities::PetscUtilities::checkError;

    // size up the aux variables
    PetscScalar *auxL = nullptr, *auxR = nullptr;

    // Precompute the offsets to pass into the rhsFluxFunctionDescriptions
    std::vector<std::vector<PetscInt>> fluxComponentSize(rhsFunctions.size());
    std::vector<std::vector<PetscInt>> fluxId(rhsFunctions.size());
    std::vector<std::vector<PetscInt>> uOff(rhsFunctions.size());
    std::vector<std::vector<PetscInt>> aOff(rhsFunctions.size());

    // Get the full set of offsets from the ds
    PetscInt* uOffTotal;
    PetscDSGetComponentOffsets(ds, &uOffTotal) >> utilities::PetscUtilities::checkError;

    for (std::size_t fun = 0; fun < rhsFunctions.size(); fun++) {
        for (std::size_t f = 0; f < rhsFunctions[fun].updateFields.size(); f++) {
            const auto& field = subDomain->GetField(rhsFunctions[fun].updateFields[f]);
            fluxComponentSize[fun].push_back(field.numberComponents);
            fluxId[fun].push_back(field.id);
        }
        for (std::size_t f = 0; f < rhsFunctions[fun].inputFields.size(); f++) {
            uOff[fun].push_back(uOffTotal[rhsFunctions[fun].inputFields[f]]);
        }
    }

    if (dsAux) {
        PetscInt* auxOffTotal;
        PetscDSGetComponentOffsets(dsAux, &auxOffTotal) >> utilities::PetscUtilities::checkError;
        for (std::size_t fun = 0; fun < rhsFunctions.size(); fun++) {
            for (std::size_t f = 0; f < rhsFunctions[fun].auxFields.size(); f++) {
                aOff[fun].push_back(auxOffTotal[rhsFunctions[fun].auxFields[f]]);
            }
        }
    }
    // check for ghost cells
    DMLabel ghostLabel;
    DMGetLabel(dm, "ghost", &ghostLabel) >> utilities::PetscUtilities::checkError;

    // get the label for this region
    DMLabel regionLabel = nullptr;
    PetscInt regionValue = 0;
    domain::Region::GetLabel(solverRegion, subDomain->GetDM(), regionLabel, regionValue);

    if (regionValue!=1){
        throw std::runtime_error("The regionValue in cellinterpollant is  " + std::to_string(regionValue) +", the main is assumed to have a regionValue of 1");

    }

//        double totalTime1=0;
//        int callCount1=0;
    //    double totalTime2=0;
    //    int callCount2=0;
    //    double totalTime3=0;
    //    int callCount3=0;
//        int output_it=1000;
//
//        double start = MPI_Wtime();
    //    start = MPI_Wtime();


    // March over each face in this region
//    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::ComputeSource::SourceOut");
    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::ComputeSource::Project");
    for (PetscInt f = faceRange.start; f < faceRange.end; ++f) {
        const PetscInt face = faceRange.points ? faceRange.points[f] : f;

        // make sure that this is a valid facel
        PetscInt ghost, nsupp, nchild;
        DMLabelGetValue(ghostLabel, face, &ghost) >> utilities::PetscUtilities::checkError;
        DMPlexGetSupportSize(dm, face, &nsupp) >> utilities::PetscUtilities::checkError;
        DMPlexGetTreeChildren(dm, face, &nchild, nullptr) >> utilities::PetscUtilities::checkError;

        if (ghost >= 0 || nsupp > 2 || nchild > 0) continue;
        // Get the face geometry
        const PetscInt* faceCells;
        PetscFVFaceGeom* fg;
        PetscFVCellGeom *cgL, *cgR;
        DMPlexPointLocalRead(faceDM, face, faceGeomArray, &fg) >> utilities::PetscUtilities::checkError;
        DMPlexGetSupport(dm, face, &faceCells) >> utilities::PetscUtilities::checkError;
        DMPlexPointLocalRead(cellDM, faceCells[0], cellGeomArray, &cgL) >> utilities::PetscUtilities::checkError;
        DMPlexPointLocalRead(cellDM, faceCells[1], cellGeomArray, &cgR) >> utilities::PetscUtilities::checkError;
//
//        PetscInt leftFlowLabelValue = regionValue;
//        PetscInt rightFlowLabelValue = regionValue;
////        start = MPI_Wtime();
//        if (regionLabel) {
//            DMLabelGetValue(regionLabel, faceCells[0], &leftFlowLabelValue);
//            DMLabelGetValue(regionLabel, faceCells[1], &rightFlowLabelValue);
//        }


//                totalTime1 += MPI_Wtime() - start;
//                ++callCount1;
//                if (callCount1==output_it) {
//                    PetscPrintf(PETSC_COMM_WORLD, "Section1 total %d time: %f s\n", callCount1, totalTime1); }

        // compute the left/right face values
        ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[0], *cgL, dm, xArray, dmGrads, locGradArrays, uL, gradL, 1);
        ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[1], *cgR, dm, xArray, dmGrads, locGradArrays, uR, gradR, 1);

//        ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[0], *cgL, dm, xArray, dmGrads, locGradArrays, uL, gradL,faceCells[1], *cgR, 1);
//        ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[1], *cgR, dm, xArray, dmGrads, locGradArrays, uR, gradR, faceCells[0], *cgL, 1);



    }
    EndEvent();




    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::ComputeSource::SourceTotal");

    bool oldcode=false;
    if (oldcode) {
        for (PetscInt f = faceRange.start; f < faceRange.end; ++f) {
            const PetscInt face = faceRange.points ? faceRange.points[f] : f;

            // make sure that this is a valid facel
            PetscInt ghost, nsupp, nchild;
            DMLabelGetValue(ghostLabel, face, &ghost) >> utilities::PetscUtilities::checkError;
            DMPlexGetSupportSize(dm, face, &nsupp) >> utilities::PetscUtilities::checkError;
            DMPlexGetTreeChildren(dm, face, &nchild, nullptr) >> utilities::PetscUtilities::checkError;

            if (ghost >= 0 || nsupp > 2 || nchild > 0) continue;
            // Get the face geometry
            const PetscInt* faceCells;
            PetscFVFaceGeom* fg;
            PetscFVCellGeom *cgL, *cgR;
            DMPlexPointLocalRead(faceDM, face, faceGeomArray, &fg) >> utilities::PetscUtilities::checkError;
            DMPlexGetSupport(dm, face, &faceCells) >> utilities::PetscUtilities::checkError;

            DMPlexPointLocalRead(cellDM, faceCells[0], cellGeomArray, &cgL) >> utilities::PetscUtilities::checkError;
            DMPlexPointLocalRead(cellDM, faceCells[1], cellGeomArray, &cgR) >> utilities::PetscUtilities::checkError;

            PetscInt leftFlowLabelValue = regionValue;
            PetscInt rightFlowLabelValue = regionValue;

            if (regionLabel) {
                DMLabelGetValue(regionLabel, faceCells[0], &leftFlowLabelValue);
                DMLabelGetValue(regionLabel, faceCells[1], &rightFlowLabelValue);
            }

            ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[0], *cgL, dm, xArray, dmGrads, locGradArrays, uL, gradL, leftFlowLabelValue == regionValue);
            ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[1], *cgR, dm, xArray, dmGrads, locGradArrays, uR, gradR, rightFlowLabelValue == regionValue);

            ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[0], *cgL, dm, xArray, dmGrads, locGradArrays, uL, gradL,faceCells[1], *cgR, leftFlowLabelValue == regionValue);
            ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[1], *cgR, dm, xArray, dmGrads, locGradArrays, uR, gradR, faceCells[0], *cgL, rightFlowLabelValue == regionValue);

            // determine the left/right cells
            if (auxArray) {
                // Get the field values at this cell
                DMPlexPointLocalRead(dmAux, faceCells[0], auxArray, &auxL) >> utilities::PetscUtilities::checkError;
                DMPlexPointLocalRead(dmAux, faceCells[1], auxArray, &auxR) >> utilities::PetscUtilities::checkError;
            }

            for (std::size_t fun = 0; fun < rhsFunctions.size(); fun++) {
                PetscInt fluxOffset = 0;  // Flux offset for the function ( Currently calculated by just adding the number of components of the previous fields)
                PetscArrayzero(flux, totDim) >> utilities::PetscUtilities::checkError;
                const auto& rhsFluxFunctionDescription = rhsFunctions[fun];
                rhsFluxFunctionDescription.function(dim, fg, uOff[fun].data(), uL, uR, aOff[fun].data(), auxL, auxR, flux, rhsFluxFunctionDescription.context) >> utilities::PetscUtilities::checkError;
                // add the fluxes back to the cell
                for (std::size_t updateFieldIdx = 0; updateFieldIdx < rhsFunctions[fun].updateFields.size(); updateFieldIdx++) {
                    PetscInt cellLabelValue = regionValue;
                    PetscScalar *fL = nullptr, *fR = nullptr;
                    DMLabelGetValue(ghostLabel, faceCells[0], &ghost) >> utilities::PetscUtilities::checkError;
                    if (regionLabel) {
                        DMLabelGetValue(regionLabel, faceCells[0], &cellLabelValue) >> utilities::PetscUtilities::checkError;
                    }
                    if (ghost <= 0 && regionValue == cellLabelValue) {
                        DMPlexPointLocalFieldRef(dm, faceCells[0], fluxId[fun][updateFieldIdx], locFArray, &fL) >> utilities::PetscUtilities::checkError;
                    }

                    cellLabelValue = regionValue;
                    DMLabelGetValue(ghostLabel, faceCells[1], &ghost) >> utilities::PetscUtilities::checkError;
                    if (regionLabel) {
                        DMLabelGetValue(regionLabel, faceCells[1], &cellLabelValue) >> utilities::PetscUtilities::checkError;
                    }
                    if (ghost <= 0 && regionValue == cellLabelValue) {
                        DMPlexPointLocalFieldRef(dm, faceCells[1], fluxId[fun][updateFieldIdx], locFArray, &fR) >> utilities::PetscUtilities::checkError;
                    }

                    for (PetscInt d = 0; d < (fluxComponentSize[fun][updateFieldIdx]); ++d) {
                        if (fL) fL[d] -= flux[fluxOffset + d] / cgL->volume;
                        if (fR) fR[d] += flux[fluxOffset + d] / cgR->volume;
                    }
                    fluxOffset += fluxComponentSize[fun][updateFieldIdx];
                }
                //            EndEvent();
            }
            //        ++callCount3;
            //        if (callCount3==output_it) {
            //            PetscPrintf(PETSC_COMM_WORLD, "Section3 total %d time: %f s\n", callCount3, totalTime3); }

            //        EndEvent();
        }
    } else {
        for (PetscInt f = faceRange.start; f < faceRange.end; ++f) {
            const PetscInt face = faceRange.points ? faceRange.points[f] : f;

            // make sure that this is a valid facel
//            PetscInt ghost2, nsupp2, nchild2;
//            DMLabelGetValue(ghostLabel, face, &ghost2) >> utilities::PetscUtilities::checkError;
//            DMPlexGetSupportSize(dm, face, &nsupp2) >> utilities::PetscUtilities::checkError;
//            DMPlexGetTreeChildren(dm, face, &nchild2, nullptr) >> utilities::PetscUtilities::checkError;
//
            PetscInt ghost = flowLabelVec[numLabel*f];
            PetscInt nsupp = flowLabelVec[numLabel*f+1];
            PetscInt nchild = flowLabelVec[numLabel*f+2];

//            if (ghost2!=ghost || nsupp2!=nsupp || nchild2!=nchild){
//                throw std::runtime_error("either ghost or nsupp or child in not matching for face:" + std::to_string(face) );
//            }

            if (ghost >= 0 || nsupp > 2 || nchild > 0) continue;
            // Get the face geometry
            const PetscInt* faceCells;
            PetscFVFaceGeom* fg;
            PetscFVCellGeom *cgL, *cgR;
            DMPlexPointLocalRead(faceDM, face, faceGeomArray, &fg) >> utilities::PetscUtilities::checkError;
            DMPlexGetSupport(dm, face, &faceCells) >> utilities::PetscUtilities::checkError;

            DMPlexPointLocalRead(cellDM, faceCells[0], cellGeomArray, &cgL) >> utilities::PetscUtilities::checkError;
            DMPlexPointLocalRead(cellDM, faceCells[1], cellGeomArray, &cgR) >> utilities::PetscUtilities::checkError;

            PetscInt leftFlowLabelValue = flowLabelVec[numLabel*f+3];
            PetscInt rightFlowLabelValue = flowLabelVec[numLabel*f+4];

//            PetscInt leftFlowLabelValue2 = regionValue;
//            PetscInt rightFlowLabelValue2 = regionValue;
//            if (regionLabel) {
//                DMLabelGetValue(regionLabel, faceCells[0], &leftFlowLabelValue2);
//                DMLabelGetValue(regionLabel, faceCells[1], &rightFlowLabelValue2);
//            }
//            if (leftFlowLabelValue!=leftFlowLabelValue2 || rightFlowLabelValue2!=rightFlowLabelValue ){
//                throw std::runtime_error("one of the face labels dont match matching for face:" + std::to_string(face) );
//            }

//            ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[0], *cgL, dm, xArray, dmGrads, locGradArrays, uL, gradL, leftFlowLabelValue == regionValue);
//            ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[1], *cgR, dm, xArray, dmGrads, locGradArrays, uR, gradR, rightFlowLabelValue == regionValue);

            ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[0], *cgL, dm, xArray, dmGrads, locGradArrays, uL, gradL,faceCells[1], *cgR, leftFlowLabelValue == regionValue);
            ProjectToFace(subDomain->GetFields(), ds, *fg, faceCells[1], *cgR, dm, xArray, dmGrads, locGradArrays, uR, gradR, faceCells[0], *cgL, rightFlowLabelValue == regionValue);


            // determine the left/right cells
            if (auxArray) {
                // Get the field values at this cell
                DMPlexPointLocalRead(dmAux, faceCells[0], auxArray, &auxL) >> utilities::PetscUtilities::checkError;
                DMPlexPointLocalRead(dmAux, faceCells[1], auxArray, &auxR) >> utilities::PetscUtilities::checkError;
            }

            for (std::size_t fun = 0; fun < rhsFunctions.size(); fun++) {
                PetscInt fluxOffset = 0;  // Flux offset for the function ( Currently calculated by just adding the number of components of the previous fields)
                PetscArrayzero(flux, totDim) >> utilities::PetscUtilities::checkError;
                const auto& rhsFluxFunctionDescription = rhsFunctions[fun];
                rhsFluxFunctionDescription.function(dim, fg, uOff[fun].data(), uL, uR, aOff[fun].data(), auxL, auxR, flux, rhsFluxFunctionDescription.context) >> utilities::PetscUtilities::checkError;
                // add the fluxes back to the cell
                for (std::size_t updateFieldIdx = 0; updateFieldIdx < rhsFunctions[fun].updateFields.size(); updateFieldIdx++) {
                    PetscInt cellLabelValue = regionValue;
//
                    PetscScalar *fL = nullptr, *fR = nullptr;

//                    PetscInt cellLabelValue2 = regionValue;
//                    DMLabelGetValue(ghostLabel, faceCells[0], &ghost2) >> utilities::PetscUtilities::checkError;
//                    if (regionLabel) {
//                        DMLabelGetValue(regionLabel, faceCells[0], &cellLabelValue2) >> utilities::PetscUtilities::checkError;
//                    }

                    ghost = flowLabelVec[numLabel*f+5];
                    cellLabelValue = flowLabelVec[numLabel*f+3];

//                    if (ghost!=ghost2 || cellLabelValue2!=cellLabelValue ){
//                        throw std::runtime_error("either ghost  or cellLabelValue2 (left) is not good for source terms, for face:" + std::to_string(face) );
//                    }

                    if (ghost <= 0 && regionValue == cellLabelValue) {
                        DMPlexPointLocalFieldRef(dm, faceCells[0], fluxId[fun][updateFieldIdx], locFArray, &fL) >> utilities::PetscUtilities::checkError;
                    }

                    cellLabelValue = regionValue;
//
//                    DMLabelGetValue(ghostLabel, faceCells[1], &ghost2) >> utilities::PetscUtilities::checkError;
//                    if (regionLabel) {
//                        DMLabelGetValue(regionLabel, faceCells[1], &cellLabelValue2) >> utilities::PetscUtilities::checkError;
//                    }


                    ghost = flowLabelVec[numLabel*f+6];
                    cellLabelValue = flowLabelVec[numLabel*f+4];

//                    if (ghost!=ghost2 || cellLabelValue2!=cellLabelValue ){
//                        throw std::runtime_error("either ghost  or cellLabelValue2 (right) is not good for source terms, for face:" + std::to_string(face) );
//                    }

                    if (ghost <= 0 && regionValue == cellLabelValue) {
                        DMPlexPointLocalFieldRef(dm, faceCells[1], fluxId[fun][updateFieldIdx], locFArray, &fR) >> utilities::PetscUtilities::checkError;
                    }

                    for (PetscInt d = 0; d < (fluxComponentSize[fun][updateFieldIdx]); ++d) {
                        if (fL) fL[d] -= flux[fluxOffset + d] / cgL->volume;
                        if (fR) fR[d] += flux[fluxOffset + d] / cgR->volume;
                    }
                    fluxOffset += fluxComponentSize[fun][updateFieldIdx];
                }
                //            EndEvent();
            }
            //        ++callCount3;
            //        if (callCount3==output_it) {
            //            PetscPrintf(PETSC_COMM_WORLD, "Section3 total %d time: %f s\n", callCount3, totalTime3); }

            //        EndEvent();
        }
    }
    EndEvent();

    // cleanup
//    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::ComputeSource::Cleanup");
    DMRestoreWorkArray(dm, totDim, MPIU_SCALAR, &flux) >> utilities::PetscUtilities::checkError;
    DMRestoreWorkArray(dm, totDim, MPIU_SCALAR, &uL) >> utilities::PetscUtilities::checkError;
    DMRestoreWorkArray(dm, totDim, MPIU_SCALAR, &uR) >> utilities::PetscUtilities::checkError;
    DMRestoreWorkArray(dm, dim * totDim, MPIU_SCALAR, &gradL) >> utilities::PetscUtilities::checkError;
    DMRestoreWorkArray(dm, dim * totDim, MPIU_SCALAR, &gradR) >> utilities::PetscUtilities::checkError;
//    EndEvent();

}

static PetscErrorCode BuildGradientReconstruction_Internal(DM dm, DMLabel regionLabel, PetscInt regionValue, PetscFV fvm, DM dmFace, PetscScalar* fgeom, DM dmCell, PetscScalar* cgeom) {
    DMLabel ghostLabel;
    PetscScalar *dx, *grad, **gref;
    PetscInt dim, cStart, cEnd, c, cEndInterior, maxNumFaces;

    PetscFunctionBegin;
    PetscCall(DMGetDimension(dm, &dim));
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    PetscCall(DMPlexGetCellTypeStratum(dm, DM_POLYTOPE_FV_GHOST, &cEndInterior, nullptr));
    cEndInterior = cEndInterior < 0 ? cEnd : cEndInterior;
    PetscCall(DMPlexGetMaxSizes(dm, &maxNumFaces, nullptr));
    PetscCall(PetscFVLeastSquaresSetMaxFaces(fvm, maxNumFaces));
    PetscCall(DMGetLabel(dm, "ghost", &ghostLabel));
    PetscCall(PetscMalloc3(maxNumFaces * dim, &dx, maxNumFaces * dim, &grad, maxNumFaces, &gref));
    for (c = cStart; c < cEndInterior; c++) {
        const PetscInt* faces;
        PetscInt numFaces, usedFaces, f, d;
        PetscFVCellGeom* cg;
        PetscBool boundary;
        PetscInt ghost;
        PetscInt labelValue;

        // do not attempt to compute a gradient reconstruction stencil in a ghost cell.  It will never be used
        PetscCall(DMLabelGetValue(ghostLabel, c, &ghost));
        if (ghost >= 0) continue;

        if (regionLabel) {
            PetscCall(DMLabelGetValue(regionLabel, c, &labelValue));
            if (labelValue != regionValue) continue;
        }

        PetscCall(DMPlexPointLocalRead(dmCell, c, cgeom, &cg));
        PetscCall(DMPlexGetConeSize(dm, c, &numFaces));
        PetscCall(DMPlexGetCone(dm, c, &faces));
        PetscCheck(!(numFaces < dim), PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Cell %" PetscInt_FMT " has only %" PetscInt_FMT " faces, not enough for gradient reconstruction", c, numFaces);
        for (f = 0, usedFaces = 0; f < numFaces; ++f) {
            PetscFVCellGeom* cg1;
            PetscFVFaceGeom* fg;
            const PetscInt* fcells;
            PetscInt ncell, side;

            if (regionLabel) {
                PetscCall(DMLabelGetValue(regionLabel, faces[f], &labelValue));
                if (labelValue != regionValue) continue;
            }

            PetscCall(DMLabelGetValue(ghostLabel, faces[f], &ghost));
            PetscCall(DMIsBoundaryPoint(dm, faces[f], &boundary));
            if ((ghost >= 0) || boundary) continue;
            PetscCall(DMPlexGetSupport(dm, faces[f], &fcells));
            side = (c != fcells[0]); /* c is on left=0 or right=1 of face */
            ncell = fcells[!side];   /* the neighbor */
            PetscCall(DMPlexPointLocalRef(dmFace, faces[f], fgeom, &fg));
            PetscCall(DMPlexPointLocalRead(dmCell, ncell, cgeom, &cg1));
            for (d = 0; d < dim; ++d) dx[usedFaces * dim + d] = cg1->centroid[d] - cg->centroid[d];
            gref[usedFaces++] = fg->grad[side]; /* Gradient reconstruction term will go here */
        }
        PetscCheck(usedFaces, PETSC_COMM_SELF, PETSC_ERR_USER, "Mesh contains isolated cell (no neighbors). Is it intentional?");
        PetscCall(PetscFVComputeGradient(fvm, usedFaces, dx, grad));
        for (f = 0, usedFaces = 0; f < numFaces; ++f) {
            if (regionLabel) {
                PetscCall(DMLabelGetValue(regionLabel, faces[f], &labelValue));
                if (labelValue != regionValue) continue;
            }
            PetscCall(DMLabelGetValue(ghostLabel, faces[f], &ghost));
            PetscCall(DMIsBoundaryPoint(dm, faces[f], &boundary));
            if ((ghost >= 0) || boundary) continue;
            for (d = 0; d < dim; ++d) gref[usedFaces][d] = grad[usedFaces * dim + d];
            ++usedFaces;
        }
    }
    // Free the memory allocated earlier with PetscMalloc3
    PetscCall(PetscFree3(dx, grad, gref));
    PetscFunctionReturn(0);
}

static PetscErrorCode BuildGradientReconstruction_Internal_Tree(DM dm, DMLabel regionLabel, PetscInt regionValue, PetscFV fvm, DM dmFace, PetscScalar* fgeom, DM dmCell, PetscScalar* cgeom) {
    DMLabel ghostLabel;
    PetscScalar *dx, *grad, **gref;
    PetscInt dim, cStart, cEnd, c, cEndInterior, fStart, fEnd, f, nStart, nEnd, maxNumFaces = 0;
    PetscSection neighSec;
    PetscInt(*neighbors)[2];
    PetscInt* counter;

    PetscFunctionBegin;
    PetscCall(DMGetDimension(dm, &dim));
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    PetscCall(DMPlexGetCellTypeStratum(dm, DM_POLYTOPE_FV_GHOST, &cEndInterior, nullptr));
    if (cEndInterior < 0) cEndInterior = cEnd;
    PetscCall(PetscSectionCreate(PetscObjectComm((PetscObject)dm), &neighSec));
    PetscCall(PetscSectionSetChart(neighSec, cStart, cEndInterior));
    PetscCall(DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd));
    PetscCall(DMGetLabel(dm, "ghost", &ghostLabel));
    for (f = fStart; f < fEnd; f++) {
        const PetscInt* fcells;
        PetscBool boundary;
        PetscInt ghost = -1;
        PetscInt numChildren, numCells, labelValue;

        if (ghostLabel) PetscCall(DMLabelGetValue(ghostLabel, f, &ghost));
        PetscCall(DMIsBoundaryPoint(dm, f, &boundary));
        PetscCall(DMPlexGetTreeChildren(dm, f, &numChildren, nullptr));
        if ((ghost >= 0) || boundary || numChildren) continue;

        if (regionLabel) {
            PetscCall(DMLabelGetValue(regionLabel, f, &labelValue));
            if (labelValue != regionValue) continue;
        }

        PetscCall(DMPlexGetSupportSize(dm, f, &numCells));
        if (numCells == 2) {
            PetscCall(DMPlexGetSupport(dm, f, &fcells));
            for (c = 0; c < 2; c++) {
                PetscInt cell = fcells[c];

                if (cell >= cStart && cell < cEndInterior) {
                    PetscCall(PetscSectionAddDof(neighSec, cell, 1));
                }
            }
        }
    }
    PetscCall(PetscSectionSetUp(neighSec));
    PetscCall(PetscSectionGetMaxDof(neighSec, &maxNumFaces));
    PetscCall(PetscFVLeastSquaresSetMaxFaces(fvm, maxNumFaces));
    nStart = 0;
    PetscCall(PetscSectionGetStorageSize(neighSec, &nEnd));
    PetscCall(PetscMalloc1((nEnd - nStart), &neighbors));
    PetscCall(PetscCalloc1((cEndInterior - cStart), &counter));
    for (f = fStart; f < fEnd; f++) {
        const PetscInt* fcells;
        PetscBool boundary;
        PetscInt ghost = -1;
        PetscInt numChildren, numCells, labelValue;

        if (ghostLabel) PetscCall(DMLabelGetValue(ghostLabel, f, &ghost));
        PetscCall(DMIsBoundaryPoint(dm, f, &boundary));
        PetscCall(DMPlexGetTreeChildren(dm, f, &numChildren, nullptr));
        if ((ghost >= 0) || boundary || numChildren) continue;

        if (regionLabel) {
            PetscCall(DMLabelGetValue(regionLabel, f, &labelValue));
            if (labelValue != regionValue) continue;
        }

        PetscCall(DMPlexGetSupportSize(dm, f, &numCells));
        if (numCells == 2) {
            PetscCall(DMPlexGetSupport(dm, f, &fcells));
            for (c = 0; c < 2; c++) {
                PetscInt cell = fcells[c], off;

                if (regionLabel) {
                    PetscCall(DMLabelGetValue(regionLabel, c, &labelValue));
                    if (labelValue != regionValue) continue;
                }

                if (cell >= cStart && cell < cEndInterior) {
                    PetscCall(PetscSectionGetOffset(neighSec, cell, &off));
                    off += counter[cell - cStart]++;
                    neighbors[off][0] = f;
                    neighbors[off][1] = fcells[1 - c];
                }
            }
        }
    }
    PetscCall(PetscFree(counter));
    PetscCall(PetscMalloc3(maxNumFaces * dim, &dx, maxNumFaces * dim, &grad, maxNumFaces, &gref));
    for (c = cStart; c < cEndInterior; c++) {
        PetscInt numFaces, d, off, labelValue, ghost = -1;
        PetscFVCellGeom* cg;

        PetscCall(DMPlexPointLocalRead(dmCell, c, cgeom, &cg));
        PetscCall(PetscSectionGetDof(neighSec, c, &numFaces));
        PetscCall(PetscSectionGetOffset(neighSec, c, &off));

        if (regionLabel) {
            PetscCall(DMLabelGetValue(regionLabel, c, &labelValue));
            if (labelValue != regionValue) continue;
        }

        // do not attempt to compute a gradient reconstruction stencil in a ghost cell.  It will never be used
        if (ghostLabel) PetscCall(DMLabelGetValue(ghostLabel, c, &ghost));
        if (ghost >= 0) continue;

        PetscCheck(!(numFaces < dim), PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Cell %" PetscInt_FMT " has only %" PetscInt_FMT " faces, not enough for gradient reconstruction", c, numFaces);
        for (f = 0; f < numFaces; ++f) {
            PetscFVCellGeom* cg1;
            PetscFVFaceGeom* fg;
            const PetscInt* fcells;
            PetscInt ncell, side, nface;

            if (regionLabel) {
                PetscCall(DMLabelGetValue(regionLabel, f, &labelValue));
                if (labelValue != regionValue) continue;
            }

            nface = neighbors[off + f][0];
            ncell = neighbors[off + f][1];
            PetscCall(DMPlexGetSupport(dm, nface, &fcells));
            side = (c != fcells[0]);
            PetscCall(DMPlexPointLocalRef(dmFace, nface, fgeom, &fg));
            PetscCall(DMPlexPointLocalRead(dmCell, ncell, cgeom, &cg1));
            for (d = 0; d < dim; ++d) dx[f * dim + d] = cg1->centroid[d] - cg->centroid[d];
            gref[f] = fg->grad[side]; /* Gradient reconstruction term will go here */
        }
        PetscCall(PetscFVComputeGradient(fvm, numFaces, dx, grad));
        for (f = 0; f < numFaces; ++f) {
            for (d = 0; d < dim; ++d) gref[f][d] = grad[f * dim + d];
        }
    }
    PetscCall(PetscFree3(dx, grad, gref));
    PetscCall(PetscSectionDestroy(&neighSec));
    PetscCall(PetscFree(neighbors));
    PetscFunctionReturn(0);
}

PetscErrorCode ablate::finiteVolume::CellInterpolant::ComputeGradientFVM(DM dm, DMLabel regionLabel, PetscInt regionValue, PetscFV fvm, Vec faceGeometry, Vec cellGeometry, DM* dmGrad) {
    DM dmFace, dmCell;
    PetscScalar *fgeom, *cgeom;
    PetscSection sectionGrad, parentSection;
    PetscInt dim, pdim, cStart, cEnd, cEndInterior, c;

    PetscFunctionBegin;
    PetscCall(DMGetDimension(dm, &dim));
    PetscCall(PetscFVGetNumComponents(fvm, &pdim));
    PetscCall(DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd));
    PetscCall(DMPlexGetCellTypeStratum(dm, DM_POLYTOPE_FV_GHOST, &cEndInterior, nullptr));
    /* Construct the interpolant corresponding to each face from the least-square solution over the cell neighborhood */
    PetscCall(VecGetDM(faceGeometry, &dmFace));
    PetscCall(VecGetDM(cellGeometry, &dmCell));
    PetscCall(VecGetArray(faceGeometry, &fgeom));
    PetscCall(VecGetArray(cellGeometry, &cgeom));
    PetscCall(DMPlexGetTree(dm, &parentSection, nullptr, nullptr, nullptr, nullptr));
    if (!parentSection) {
        PetscCall(BuildGradientReconstruction_Internal(dm, regionLabel, regionValue, fvm, dmFace, fgeom, dmCell, cgeom));
    } else {
        PetscCall(BuildGradientReconstruction_Internal_Tree(dm, regionLabel, regionValue, fvm, dmFace, fgeom, dmCell, cgeom));
    }
    PetscCall(VecRestoreArray(faceGeometry, &fgeom));
    PetscCall(VecRestoreArray(cellGeometry, &cgeom));
    /* Create storage for gradients */
    PetscCall(DMClone(dm, dmGrad));
    PetscCall(PetscSectionCreate(PetscObjectComm((PetscObject)dm), &sectionGrad));
    PetscCall(PetscSectionSetChart(sectionGrad, cStart, cEnd));
    for (c = cStart; c < cEnd; ++c) PetscCall(PetscSectionSetDof(sectionGrad, c, pdim * dim));
    PetscCall(PetscSectionSetUp(sectionGrad));
    PetscCall(DMSetLocalSection(*dmGrad, sectionGrad));
    PetscCall(PetscSectionDestroy(&sectionGrad));
    PetscFunctionReturn(0);
}


void ablate::finiteVolume::CellInterpolant::ProjectToFace(const std::vector<domain::Field>& fields, PetscDS ds, const PetscFVFaceGeom& faceGeom, PetscInt cellId,  const PetscFVCellGeom& cellGeom,
                                                          DM dm, const PetscScalar* xArray, const std::vector<DM>& dmGrads, const std::vector<const PetscScalar*>& gradArrays, PetscScalar* u,
                                                          PetscScalar* grad, bool projectField) {
//    StartEvent("FiniteVolumeSolver::CellInterpolant::ComputeRHS::ProjectToFace");

        //Timing
//        double start = MPI_Wtime();


            //Papi low level
//        int EventSet = PAPI_NULL;
//        long long values[1];
//        PAPI_create_eventset(&EventSet);
//        PAPI_add_event(EventSet, PAPI_DP_OPS);
//        PAPI_start(EventSet);


        //PAPI high level
//        int retval;
//        retval = PAPI_hl_region_begin("project");

    const auto dim = subDomain->GetDimensions();
    // [R: 1] — Read subDomain

    // Keep track of derivative offset
    PetscInt* offsets;
    PetscInt* dirOffsets;
    PetscDSGetComponentOffsets(ds, &offsets) >> utilities::PetscUtilities::checkError;
    PetscDSGetComponentDerivativeOffsets(ds, &dirOffsets) >> utilities::PetscUtilities::checkError;


    // March over each field
    for (const auto& field : fields) {
        // [R: 1], [W: 1] per loop variable
        PetscReal dx[3];
        PetscScalar* xCell;
        PetscScalar* gradCell;

        // Get the field values at this cell
        DMPlexPointLocalFieldRead(dm, cellId, field.subId, xArray, &xCell) >> utilities::PetscUtilities::checkError;

        // If we need to project the field
        if (projectField && dmGrads[field.subId]) {
            DMPlexPointLocalRead(dmGrads[field.subId], cellId, gradArrays[field.subId], &gradCell) >> utilities::PetscUtilities::checkError;
            DMPlex_WaxpyD_Internal(dim, -1, cellGeom.centroid, faceGeom.centroid, dx);

            // Project the cell centered value onto the face
            for (PetscInt c = 0; c < field.numberComponents; ++c) {
                u[offsets[field.subId] + c] = xCell[c] + DMPlex_DotD_Internal(dim, &gradCell[c * dim], dx);

                // copy the gradient into the grad vector
                for (PetscInt d = 0; d < dim; d++) {
                    // [R: 1], [W: 1] per loop variable
                    grad[dirOffsets[field.subId] + c * dim + d] = gradCell[c * dim + d];
                    // [R: 1], [W: 1] per dimension
                }
            }

        } else if (dmGrads[field.subId]) {
            // Project the cell centered value onto the face
            DMPlexPointLocalRead(dmGrads[field.subId], cellId, gradArrays[field.subId], &gradCell) >> utilities::PetscUtilities::checkError;
            // Project the cell centered value onto the face
            for (PetscInt c = 0; c < field.numberComponents; ++c) {
                u[offsets[field.subId] + c] = xCell[c];

                // copy the gradient into the grad vector
                for (PetscInt d = 0; d < dim; d++) {
                    grad[dirOffsets[field.subId] + c * dim + d] = gradCell[c * dim + d];
                }
            }

        } else {
            // Just copy the cell centered value on to the face
            for (PetscInt c = 0; c < field.numberComponents; ++c) {
                u[offsets[field.subId] + c] = xCell[c];

                // fill the grad with NAN to prevent use
                for (PetscInt d = 0; d < dim; d++) {
                    grad[dirOffsets[field.subId] + c * dim + d] = NAN;
                }
            }
        }
    }

    //PAPI high level
//        retval = PAPI_hl_region_end("project");
//        if ( retval != PAPI_OK ){
//            handle_error(retval);
//        }

        //Timing function
//            totalTime += MPI_Wtime() - start;
//            ++callCount;
//            PetscPrintf(PETSC_COMM_WORLD, "StaticFunction called %d times, total time: %f s\n", callCount, totalTime);


        //    PAPI low level
        //        PAPI_stop(EventSet, values);
        //        printf("FLOPs counted: %lld\n", values[0]);
        //    if (PAPI_query_event(PAPI_DP_OPS) != PAPI_OK) {
        //        fprintf(stderr, "PAPI_FP_OPS not supported on this architecture\n");
        //        exit(1);
        //    }
        //    // Clean up
        //    PAPI_cleanup_eventset(EventSet);
        //    PAPI_destroy_eventset(&EventSet);

//    EndEvent();


}


void ablate::finiteVolume::CellInterpolant::ProjectToFace(const std::vector<domain::Field>& fields, PetscDS ds, const PetscFVFaceGeom& faceGeom, PetscInt cellId, const PetscFVCellGeom& cellGeom,
                                                          DM dm, const PetscScalar* xArray, const std::vector<DM>& dmGrads, const std::vector<const PetscScalar*>& gradArrays, PetscScalar* u,
                                                          PetscScalar* grad, PetscInt neighborCellId, const PetscFVCellGeom& neighborCellGeom, bool projectField) {

    const auto dim = subDomain->GetDimensions();

    // Keep track of derivative offset
    PetscInt* offsets;
    PetscInt* dirOffsets;
    PetscDSGetComponentOffsets(ds, &offsets) >> utilities::PetscUtilities::checkError;
    PetscDSGetComponentDerivativeOffsets(ds, &dirOffsets) >> utilities::PetscUtilities::checkError;

    // March over each field
    for (const auto& field : fields) {
        PetscReal dx[3];
        PetscScalar* xCell;
        PetscScalar* gradCell;

        // Get the field values at this cell
        DMPlexPointLocalFieldRead(dm, cellId, field.subId, xArray, &xCell) >> utilities::PetscUtilities::checkError;

        // If we need to project the field
        if (projectField && dmGrads[field.subId]) {
            DMPlexPointLocalRead(dmGrads[field.subId], cellId, gradArrays[field.subId], &gradCell) >> utilities::PetscUtilities::checkError;
            DMPlex_WaxpyD_Internal(dim, -1, cellGeom.centroid, faceGeom.centroid, dx);


            // Apply limiter to the gradient
            auto fvm = (PetscFV)subDomain->GetPetscFieldObject(field);
            PetscLimiter lim;
            PetscFVGetLimiter(fvm, &lim) >> utilities::PetscUtilities::checkError;

            PetscReal limitedGrad[dim * field.numberComponents];
            if (lim && neighborCellId >= 0 ) {
                // Get neighbor cell values
                PetscScalar* neighborXCell;
                DMPlexPointLocalFieldRead(dm, neighborCellId, field.subId, xArray, &neighborXCell) >> utilities::PetscUtilities::checkError;

                // Calculate distance vector between cell centers
                PetscReal cellDistance[3];
                // v_i = NeighborCentroid_i - ThisCentroid_i = dx_i
                DMPlex_WaxpyD_Internal(dim, -1, cellGeom.centroid, neighborCellGeom.centroid, cellDistance);

                // Apply limiting per component based on this specific face
                for (PetscInt c = 0; c < field.numberComponents; ++c) {
                    PetscReal phi = 1.0;  // Default to no limiting

                    // Calculate the denom(gradient dot distance)
                    PetscReal denom = 0.0;
                    for (PetscInt d = 0; d < dim; ++d) {
                        denom += gradCell[c * dim + d] * cellDistance[d];
                    }

                    if (PetscAbsReal(denom) > PETSC_SMALL) {
                        // Use symmetric slope limiter form (Berger, Aftosmis, and Murman 2005)
                        PetscReal flim = 0.5 * PetscRealPart(neighborXCell[c] - xCell[c]) / denom;
                        PetscLimiterLimit(lim, flim, &phi) >> utilities::PetscUtilities::checkError;
                    }

//                    // Apply additional maxGradient limiting if needed
//                    PetscBool cancelGrad = PETSC_FALSE;
//                    if (phi == 0.0) {
//                        for (PetscInt d = 0; d < dim; d++) {
//                            if (PetscAbsReal(gradCell[c * dim + d]) > maxLimGrad) {
//                                cancelGrad = PETSC_TRUE;
//                                break;
//                            }
//                        }
//                    }
//                    // Apply limiting to gradient components
//                    for (PetscInt d = 0; d < dim; ++d) {
//                        if (cancelGrad) {
//                            limitedGrad[c * dim + d] = 0.0;
//                        } else {
//                            limitedGrad[c * dim + d] = phi * gradCell[c * dim + d];
//                        }
//                    }
//                }

                    // Apply limiting to gradient components
                    for (PetscInt d = 0; d < dim; ++d) {
                        limitedGrad[c * dim + d] = phi * gradCell[c * dim + d];
                    }
                }

                for (PetscInt c = 0; c < field.numberComponents; ++c) {
                    PetscReal projection = 0.0;
                    for (PetscInt d = 0; d < dim; ++d) {
                        projection += limitedGrad[c * dim + d] * dx[d];
                    }
                    u[offsets[field.subId] + c] = xCell[c] + projection;

                    // Copy the limited gradient into the grad vector
                    for (PetscInt d = 0; d < dim; d++) {
                        grad[dirOffsets[field.subId] + c * dim + d] = limitedGrad[c * dim + d];
                    }
                }
            } else {
                // No limiting - use original gradient directly
                for (PetscInt c = 0; c < field.numberComponents; ++c) {
                    u[offsets[field.subId] + c] = xCell[c] + DMPlex_DotD_Internal(dim, &gradCell[c * dim], dx);

                    // Copy the gradient into the grad vector
                    for (PetscInt d = 0; d < dim; d++) {
                        grad[dirOffsets[field.subId] + c * dim + d] = gradCell[c * dim + d];
                    }
                }
            }


        } else if (dmGrads[field.subId]) {
            // Copy cell centered value and gradient without projection
            DMPlexPointLocalRead(dmGrads[field.subId], cellId, gradArrays[field.subId], &gradCell) >> utilities::PetscUtilities::checkError;
            for (PetscInt c = 0; c < field.numberComponents; ++c) {
                u[offsets[field.subId] + c] = xCell[c];
                for (PetscInt d = 0; d < dim; d++) {
                    grad[dirOffsets[field.subId] + c * dim + d] = gradCell[c * dim + d];
                }
            }

        } else {
            // Just copy the cell centered value on to the face
            for (PetscInt c = 0; c < field.numberComponents; ++c) {
                u[offsets[field.subId] + c] = xCell[c];

                // fill the grad with NAN to prevent use
                for (PetscInt d = 0; d < dim; d++) {
                    grad[dirOffsets[field.subId] + c * dim + d] = NAN;
                }
            }
        }
    }
}