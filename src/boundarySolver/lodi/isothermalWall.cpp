#include "isothermalWall.hpp"
#include "finiteVolume/compressibleFlowFields.hpp"
#include "utilities/mathUtilities.hpp"

using fp = ablate::finiteVolume::CompressibleFlowFields;

ablate::boundarySolver::lodi::IsothermalWall::IsothermalWall(std::shared_ptr<eos::EOS> eos, std::shared_ptr<finiteVolume::processes::PressureGradientScaling> pressureGradientScaling)
    : LODIBoundary(std::move(eos), std::move(pressureGradientScaling)) {}

void ablate::boundarySolver::lodi::IsothermalWall::Setup(ablate::boundarySolver::BoundarySolver &bSolver) {
    ablate::boundarySolver::lodi::LODIBoundary::Setup(bSolver);
    bSolver.RegisterFunction(IsothermalWallFunction, this, fieldNames, fieldNames, {});

    bSolver.RegisterPreRHSFunction(CorrectBoundaryEnergy, this);

    if (nSpecEqs) {
        bSolver.RegisterFunction(
            MirrorSpecies, this, {finiteVolume::CompressibleFlowFields::EULER_FIELD, finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD}, {finiteVolume::CompressibleFlowFields::YI_FIELD});
    }
}

PetscErrorCode ablate::boundarySolver::lodi::IsothermalWall::IsothermalWallFunction(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg,
                                                                                    const PetscFVCellGeom *boundaryCell, const PetscInt uOff[], const PetscScalar *boundaryValues,
                                                                                    const PetscScalar *stencilValues[], const PetscInt aOff[], const PetscScalar *auxValues,
                                                                                    const PetscScalar *stencilAuxValues[], PetscInt stencilSize, const PetscInt stencil[],
                                                                                    const PetscScalar stencilWeights[], const PetscInt sOff[], PetscScalar source[], void *ctx) {
    PetscFunctionBeginUser;
    auto isothermalWall = (IsothermalWall *)ctx;

    // Compute the transformation matrix
    PetscReal transformationMatrix[3][3];
    utilities::MathUtilities::ComputeTransformationMatrix(dim, fg->normal, transformationMatrix);

    // Compute the pressure/values on the boundary
    PetscReal boundaryDensity;
    PetscReal boundaryTemperature;
    PetscReal boundaryVel[3];
    PetscReal boundaryNormalVelocity = 0.0;
    PetscReal boundarySpeedOfSound;
    PetscReal boundaryPressure;

//    boundaryTemperature=300;

    // Get the velocity and pressure on the surface
    {
        boundaryDensity = boundaryValues[uOff[isothermalWall->eulerId] + finiteVolume::CompressibleFlowFields::RHO];
        for (PetscInt d = 0; d < dim; d++) {
            boundaryVel[d] = boundaryValues[uOff[isothermalWall->eulerId] + finiteVolume::CompressibleFlowFields::RHOU + d] / boundaryDensity;
            boundaryNormalVelocity += boundaryVel[d] * fg->normal[d];
            boundaryNormalVelocity += boundaryVel[d] * fg->normal[d];
        }
        PetscCall(isothermalWall->computeTemperature.function(boundaryValues, &boundaryTemperature, isothermalWall->computeTemperature.context.get()));
//        boundaryTemperature=300;
        PetscCall(isothermalWall->computeSpeedOfSound.function(boundaryValues, boundaryTemperature, &boundarySpeedOfSound, isothermalWall->computeSpeedOfSound.context.get()));
        PetscCall(isothermalWall->computePressureFromTemperature.function(boundaryValues, boundaryTemperature, &boundaryPressure, isothermalWall->computePressureFromTemperature.context.get()));
    }

    // Map the boundary velocity into the normal coord system
    PetscReal boundaryVelNormCord[3];
    utilities::MathUtilities::Multiply(dim, transformationMatrix, boundaryVel, boundaryVelNormCord);

    // Compute each stencil point
    std::vector<PetscReal> stencilDensity(stencilSize);
    std::vector<std::vector<PetscReal>> stencilVel(stencilSize, std::vector<PetscReal>(dim));
    std::vector<PetscReal> stencilNormalVelocity(stencilSize);
    std::vector<PetscReal> stencilPressure(stencilSize);

    for (PetscInt s = 0; s < stencilSize; s++) {
        stencilDensity[s] = stencilValues[s][uOff[isothermalWall->eulerId] + finiteVolume::CompressibleFlowFields::RHO];
        for (PetscInt d = 0; d < dim; d++) {
            stencilVel[s][d] = stencilValues[s][uOff[isothermalWall->eulerId] + finiteVolume::CompressibleFlowFields::RHOU + d] / stencilDensity[s];
            stencilNormalVelocity[s] += stencilVel[s][d] * fg->normal[d];
        }
        PetscCall(isothermalWall->computePressure.function(stencilValues[s], &stencilPressure[s], isothermalWall->computePressure.context.get()));
    }

    // Interpolate the normal velocity gradient to the surface
    PetscScalar dVeldNorm;
    BoundarySolver::ComputeGradientAlongNormal(dim, fg, boundaryNormalVelocity, stencilSize, &stencilNormalVelocity[0], stencilWeights, dVeldNorm);
    PetscScalar dPdNorm;
    BoundarySolver::ComputeGradientAlongNormal(dim, fg, boundaryPressure, stencilSize, &stencilPressure[0], stencilWeights, dPdNorm);

    PetscReal boundaryCp, boundaryCv;
    isothermalWall->computeSpecificHeatConstantPressure.function(boundaryValues, boundaryTemperature, &boundaryCp, isothermalWall->computeSpecificHeatConstantPressure.context.get());
    isothermalWall->computeSpecificHeatConstantVolume.function(boundaryValues, boundaryTemperature, &boundaryCv, isothermalWall->computeSpecificHeatConstantVolume.context.get());

    // Compute the enthalpy
    PetscReal boundarySensibleEnthalpy;
    isothermalWall->computeSensibleEnthalpyFunction.function(boundaryValues, boundaryTemperature, &boundarySensibleEnthalpy, isothermalWall->computeSensibleEnthalpyFunction.context.get());

    // get_vel_and_c_prims(PGS, velwall[0], C, Cp, Cv, velnprm, Cprm);
    PetscReal velNormPrim, speedOfSoundPrim;
    isothermalWall->GetVelAndCPrims(boundaryNormalVelocity, boundarySpeedOfSound, boundaryCp, boundaryCv, velNormPrim, speedOfSoundPrim);

    // get_eigenvalues
    std::vector<PetscReal> lambda(isothermalWall->nEqs);
    isothermalWall->GetEigenValues(boundaryNormalVelocity, boundarySpeedOfSound, velNormPrim, speedOfSoundPrim, &lambda[0]);

    // Compute alpha2
    PetscReal alpha2 = 1.0;
    if (isothermalWall->pressureGradientScaling) {
        alpha2 = PetscSqr(isothermalWall->pressureGradientScaling->GetAlpha());
    }

    // Get scriptL
    std::vector<PetscReal> scriptL(isothermalWall->nEqs);
    scriptL[1 + dim] = lambda[1 + dim] * (dPdNorm - boundaryDensity * dVeldNorm * alpha2 * (velNormPrim - boundaryNormalVelocity - speedOfSoundPrim));  // Outgoing
    scriptL[1 + dim] = 0; //L3
    // acoustic
    // wave
    scriptL[0] = scriptL[1 + dim];  // Incoming acoustic wave
    // sL[1][n1][n0] = 0.; // Entropy wave
    scriptL[1] = 0.5e+0 * (boundaryCp / boundaryCv - 1.e+0) * (scriptL[1 + dim] + scriptL[0]) -
                 (boundaryCp / boundaryCv + 1.e+0) * (scriptL[0] - scriptL[1 + dim]) * (velNormPrim - boundaryNormalVelocity) / speedOfSoundPrim;  // Entropy wave
    for (int d = 1; d < dim; d++) {
        scriptL[1 + d] = 0.e+0;  // Tangential velocities
    }
    // Species
    for (int ns = 0; ns < isothermalWall->nSpecEqs; ns++) {
        scriptL[2 + dim + ns] = 0.e+0;
    }
    // Extra variables
    for (int ne = 0; ne < isothermalWall->nEvEqs; ne++) {
        scriptL[2 + dim + isothermalWall->nSpecEqs + ne] = 0.e+0;
    }

    // Get the pointers to the ev fields

    // Directly compute the source terms, note that this may be problem in the future with multiple source terms on the same boundary cell
    isothermalWall->GetmdFdn(sOff,
                             boundaryVelNormCord,
                             boundaryDensity,
                             boundaryTemperature,
                             boundaryCp,
                             boundaryCv,
                             boundarySpeedOfSound,
                             boundarySensibleEnthalpy,
                             velNormPrim,
                             speedOfSoundPrim,
                             boundaryValues,
                             uOff,
                             scriptL.data(),
                             transformationMatrix,
                             source);
//    for (PetscInt d = 0; d < (2+dim+isothermalWall->nSpecEqs); d++) {
//        source[d]=0;
//    }

    PetscFunctionReturn(0);
}
PetscErrorCode ablate::boundarySolver::lodi::IsothermalWall::CorrectBoundaryEnergy(
    ablate::boundarySolver::BoundarySolver& solver,
    TS ts, PetscReal time, bool initialStage, Vec locX, void* ctx) {

    PetscFunctionBeginUser;

    auto isothermalWall = reinterpret_cast<IsothermalWall*>(ctx);
    auto& subDomain = solver.GetSubDomain();

    // Get field information
    const auto& eulerField = subDomain.GetField(finiteVolume::CompressibleFlowFields::EULER_FIELD);
    const auto& densityYiField = subDomain.GetField(finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD);

    // Get the DM and dimension
    DM dm = subDomain.GetDM();
    PetscInt dim = subDomain.GetDimensions();

    // Get array access to the local vector
    PetscScalar* locXArray;
    PetscCall(VecGetArray(locX, &locXArray));

    // March over each boundary cell in the solver region
    ablate::domain::Range cellRange;
    solver.GetCellRange(cellRange);

    for (PetscInt c = cellRange.start; c < cellRange.end; ++c) {
        PetscInt boundaryCell = cellRange.points ? cellRange.points[c] : c;

        // Get pointer to the boundary cell data
        PetscScalar* cellData;
        PetscCall(DMPlexPointLocalRef(dm, boundaryCell, locXArray, &cellData));

        if (cellData) {
            // Extract conserved variables
            PetscReal rho_old = cellData[eulerField.offset + finiteVolume::CompressibleFlowFields::RHO];
//            PetscReal rhoE = cellData[eulerField.offset + finiteVolume::CompressibleFlowFields::RHOE];

            if (rho_old<0){
                std::cout << "The density is negative for cell: " << boundaryCell  << " \n";
            }

            PetscReal rho = PetscMax(0.005, rho_old);
            rho = PetscMin(100, rho);
            // TODO something smarter maybe? ...

            //Reset the density
            cellData[eulerField.offset + finiteVolume::CompressibleFlowFields::RHO] = rho;
//            cellData[eulerField.offset + fp::RHO] = rho;
            //Reset the momentum
            PetscReal mom;
            for (PetscInt d = 0; d < dim; d++) {
                mom=cellData[eulerField.offset + finiteVolume::CompressibleFlowFields::RHOU+d];
                cellData[eulerField.offset + finiteVolume::CompressibleFlowFields::RHOU+d] = mom / rho_old * rho;
            }

            // Reset the Species
            PetscReal yi;
            for (PetscInt sp = 0; sp < isothermalWall->nSpecEqs; sp++) {
                yi=cellData[densityYiField.offset+sp];
                cellData[densityYiField.offset+sp] = yi / rho_old * rho;

            }


            // Calculate kinetic energy
            PetscReal KE = 0.0;
            for (PetscInt d = 0; d < dim; ++d) {
                PetscReal momentum_d = cellData[eulerField.offset + finiteVolume::CompressibleFlowFields::RHOU + d];
                KE += momentum_d * momentum_d;
            }
            KE = 0.5 * KE / rho;

            // Sensible energy
//            PetscReal e_current = rhoE / rho_old - KE;


            PetscReal e_tmp;
            PetscCall(isothermalWall->computeInternalEnergyFromTemperature.function(cellData + eulerField.offset,
                isothermalWall->wallTemperature, &e_tmp, isothermalWall->computeInternalEnergyFromTemperature.context.get()));

//            PetscReal e_corrected = PetscMax(e_tmp, e_current);
//            PetscReal e_corrected_total = e_corrected+ KE;

            e_tmp = e_tmp + KE;

            cellData[eulerField.offset + finiteVolume::CompressibleFlowFields::RHOE] = rho * (e_tmp);





        }
    }

    // Don't forget to restore the range
    solver.RestoreRange(cellRange);

    PetscCall(VecRestoreArray(locX, &locXArray));

    PetscFunctionReturn(0);
}


PetscErrorCode ablate::boundarySolver::lodi::IsothermalWall::MirrorSpecies(PetscInt dim, const ablate::boundarySolver::BoundarySolver::BoundaryFVFaceGeom *fg, const PetscFVCellGeom *boundaryCell,
                                                                           const PetscInt *uOff, PetscScalar *boundaryValues, const PetscScalar *stencilValues, const PetscInt *aOff,
                                                                           PetscScalar *auxValues, const PetscScalar *stencilAuxValues, void *ctx) {
    PetscFunctionBeginUser;
    auto isothermalWall = (IsothermalWall *)ctx;
    const PetscInt EULER_FIELD = 0;
    const PetscInt DENSITY_YI = 1;
    const PetscInt YI = 0;

    PetscScalar boundaryDensity = boundaryValues[uOff[EULER_FIELD] + RHO];
    PetscScalar stencilDensity = stencilValues[uOff[EULER_FIELD] + RHO];
    for (PetscInt sp = 0; sp < isothermalWall->nSpecEqs; sp++) {
        PetscScalar yi = stencilValues[uOff[DENSITY_YI] + sp] / stencilDensity;

        boundaryValues[uOff[DENSITY_YI] + sp] = yi * boundaryDensity;
        auxValues[aOff[YI] + sp] = yi;
    }
//    boundaryValues[uOff[EULER_FIELD] + RHO] = stencilValues[uOff[EULER_FIELD] + RHO];
//    boundaryValues[uOff[EULER_FIELD] + RHO+1] = stencilValues[uOff[EULER_FIELD] + RHO+1];
//    boundaryValues[uOff[EULER_FIELD] + RHO+2] = stencilValues[uOff[EULER_FIELD] + RHO+2];
//    boundaryValues[uOff[EULER_FIELD] + RHO+3] = stencilValues[uOff[EULER_FIELD] + RHO+3];
//    boundaryValues[uOff[EULER_FIELD] + RHO+4] = stencilValues[uOff[EULER_FIELD] + RHO+4];

    PetscFunctionReturn(0);
}
#include "registrar.hpp"
REGISTER(ablate::boundarySolver::BoundaryProcess, ablate::boundarySolver::lodi::IsothermalWall, "Enforces a isothermal wall with fixed velocity/temperature",
         ARG(ablate::eos::EOS, "eos", "The EOS describing the flow field at the wall"),
         OPT(ablate::finiteVolume::processes::PressureGradientScaling, "pgs", "Pressure gradient scaling is used to scale the acoustic propagation speed and increase time step for low speed flows"));
