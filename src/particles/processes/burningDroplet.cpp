#include "burningDroplet.hpp"
#include "particles/particleSolver.hpp"

#include <utility>

ablate::particles::processes::BurningDroplet::BurningDroplet(std::string coupledFieldName, std::shared_ptr<mathFunctions::MathFunction> sourceFunction,const std::shared_ptr<parameters::Parameters>& parameters)
    : coupledFieldName(std::move(coupledFieldName)), sourceFunction(std::move(sourceFunction)) {}

void ablate::particles::processes::BurningDroplet::ComputeEulerianSource(PetscReal startTime, PetscReal endTime, ablate::particles::accessors::SwarmAccessor& swarmAccessorPreStep,
                                                                                  ablate::particles::accessors::SwarmAccessor& swarmAccessorPostStep,
                                                                                  ablate::particles::accessors::EulerianSourceAccessor& eulerianSourceAccessor) {



    auto diam = swarmAccessorPreStep.GetData("ParticleDiameter");

    // Get sizes from the accessors
    const auto np = swarmAccessorPreStep.GetNumberParticles();
    auto coordinates = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleCoordinates];
    auto source = eulerianSourceAccessor[coupledFieldName];
    auto dt = endTime - startTime;

    // Get the function as a petsc function
    auto function = sourceFunction->GetPetscFunction();
    auto context = sourceFunction->GetContext();

    // Store the result in a scratch variable
    PetscReal sourceValues[source.numberComponents];

    // Compute the function
    for (PetscInt p = 0; p < np; ++p) {
        function(coordinates.numberComponents, startTime, coordinates[p], source.numberComponents, sourceValues, context) >> utilities::PetscUtilities::checkError;

        // add to the source accessor and multiply by dt
        for (PetscInt c = 0; c < source.numberComponents; ++c) {
            source(p, c) += sourceValues[c] * dt;
        }
    }
}




void ablate::particles::processes::BurningDroplet::densityYiSource(PetscReal startTime, PetscReal endTime, ablate::particles::accessors::SwarmAccessor& swarmAccessorPreStep,
                                                                   ablate::particles::accessors::SwarmAccessor& swarmAccessorPostStep,
                                                                   ablate::particles::accessors::EulerianSourceAccessor& eulerianSourceAccessor) {

    //TODO add species source

}

void ablate::particles::processes::BurningDroplet::stepaux(PetscReal startTime, PetscReal endTime,  accessors::SwarmAccessor& swarmAccessorPreStep,
                                                           accessors::SwarmAccessor& swarmAccessorPostStep,
                                                           ablate::particles::accessors::EulerianAccessor& eulerianAccessor,
                                                           accessors::EulerianSourceAccessor& eulerianSourceAccessor){


    // Get sizes from the accessors
    const auto dim = eulerianAccessor.GetDimensions();
    const auto np = swarmAccessorPreStep.GetNumberParticles();
    auto coordinates = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleCoordinates];
    auto source = eulerianSourceAccessor[coupledFieldName];
    auto dt = endTime - startTime;


    // Particle

    auto partDiam = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleDiameter];
    auto partVel = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleVelocity];
    auto partDens = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleDensity];
    auto partTemp = swarmAccessorPreStep[ablate::particles::ParticleSolver::ParticleTemperature];

    // Fluid



    //TODO calculate these from the eos for the drag model
    PetscReal muF = 3.178e-5;
    PetscReal rhoF = 1.293;

    PetscReal dragF;

    // Calculate drag and step particle
    dragprocess->ComputeDragForce(partVel.dataSize,partVel[0],partVel[0],muF,rhoF,*partDiam[0],&dragF);


    for (PetscInt p = 0; p < np; ++p) {
        PetscReal rep = 0.0;
        PetscReal corFactor;
        PetscScalar tauP;

        for (PetscInt n = 0; n < dim; n++) {
            rep += rhoF * PetscSqr(fluidVel(p, n) - partVel(p, n)) * partDiam(p) / muF;
        }
        // Correction factor to account for finite Rep on Stokes drag (see Schiller-Naumann drag closure)
        corFactor = 1.0 + 0.15 * PetscPowReal(PetscSqrtReal(rep), 0.687);
        if (rep < 0.1) {
            corFactor = 1.0;  // returns Stokes drag for low speed particles
        }
        // Note: this function assumed that the solution vector order is correct
        tauP = partDens(p) * PetscSqr(partDiam(p)) / (18.0 * muF);  // particle relaxation time
        for (PetscInt n = 0; n < dim; n++) {
            coordinateRhs(p, n) += partVel(p, n);
            velocityRhs(p, n) += corFactor * (fluidVel(p, n) - partVel(p, n)) / tauP + gravityField[n] * (1.0 - rhoF / partDens(p));
        }
    }




    for (PetscInt p = 0; p < np; ++p) {
        burnprocess->ComputeBurnRate(*farfield,  *burnRate,  *energySource)


    }


    // Get the function as a petsc function
    auto function = sourceFunction->GetPetscFunction();
    auto context = sourceFunction->GetContext();

    // Store the result in a scratch variable
    PetscReal sourceValues[source.numberComponents];

    // Compute the function
    for (PetscInt p = 0; p < np; ++p) {
        function(coordinates.numberComponents, startTime, coordinates[p], source.numberComponents, sourceValues, context) >> utilities::PetscUtilities::checkError;

    // add to the source accessor and multiply by dt
    for (PetscInt c = 0; c < source.numberComponents; ++c) {
        source(p, c) += sourceValues[c] * dt;
    }
}





};










#include "registrar.hpp"
REGISTER(ablate::particles::processes::Process, ablate::particles::processes::BurningDroplet, "adds an arbitrary source function for each particle to the Eulerian field",
         ARG(std::string, "coupledField", "the name of the Eulerian coupled field"),
         ARG(ablate::mathFunctions::MathFunction, "sourceFunction", "the function to compute the source"),
         ARG(ablate::parameters::Parameters, "parameters", "fluid parameters for the particles (fluidDensity, fluidViscosity, gravityField)"));
