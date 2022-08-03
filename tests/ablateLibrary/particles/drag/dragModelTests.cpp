#include <functional>
#include <memory>
#include "PetscTestFixture.hpp"
#include "gtest/gtest.h"
#include "particles/drag/linear.hpp"
#include "particles/drag/stokes.hpp"

namespace ablateTesting::particles::drag {
struct DragModelTestParameters {
    std::function<std::shared_ptr<ablate::particles::drag::DragModel>()> createDragModel;
    std::vector<PetscReal> partVel;
    std::vector<PetscReal> flowVel;
    PetscReal muF;
    PetscReal rhoF;
    PetscReal partDiam;
    PetscReal partDens;
    PetscReal Rep;
    std::vector<PetscReal> expectedDragForce;
};

class DragModelTestFixture : public testingResources::PetscTestFixture, public ::testing::WithParamInterface<DragModelTestParameters> {};

TEST_P(DragModelTestFixture, ShouldComputeCorrectDragForce) {
    // arrange
    const auto& param = GetParam();
    auto dragModel = std::make_shared<ablate::particles::drag::Linear>();
    std::vector<PetscReal> computedDragForce(param.expectedDragForce.size());

    // act
    dragModel->ComputeDragForce(param.expectedDragForce.size(), param.partVel.data(), param.flowVel.data(), param.muF, param.rhoF, param.partDiam, param.partDens,param.Rep, computedDragForce.data());

    // assert
    for (std::size_t i = 0; i < param.expectedDragForce.size(); i++) {
        ASSERT_DOUBLE_EQ(param.expectedDragForce[i], computedDragForce[i]);
    }
}

INSTANTIATE_TEST_SUITE_P(DragModelTests, DragModelTestFixture,
                         testing::Values((DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0, 0.0},
                                                                   .flowVel = {10.0, 10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 1.0e-3,
                                                                   .partDens = 800.0,
                                                                   .Rep=490.0,
                                                                   .expectedDragForce = {84.6639,84.6639}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0, 0.0, 0.0},
                                                                   .flowVel = {10.0, 10.0, 10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 1.0e-3,
                                                                   .partDens = 800.0,
                                                                   .Rep=735.0,
                                                                   .expectedDragForce = {84.6639,84.6639,84.6639}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0},
                                                                   .flowVel = {10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 1e-3,
                                                                   .partDens = 800.0,
                                                                   .Rep=245.0,
                                                                   .expectedDragForce = {84.6639}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {10.0},
                                                                   .flowVel = {0.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 1e-3,
                                                                   .partDens = 800.0,
                                                                   .Rep=245.0,
                                                                   .expectedDragForce = {-84.6639}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {10.0},
                                                                   .flowVel = {10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 1e-3,
                                                                   .partDens = 800.0,
                                                                   .Rep=0,
                                                                   .expectedDragForce = {0}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0},
                                                                   .flowVel = {-10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 1e-3,
                                                                   .partDens = 800.0,
                                                                   .Rep=245.0,
                                                                   .expectedDragForce = {-84.6639}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0},
                                                                   .flowVel = {0.000000001},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 1e-3,
                                                                   .partDens = 800.0,
                                                                   .Rep=2.4500e-08,
                                                                   .expectedDragForce = {2.7495e-12}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0},
                                                                   .flowVel = {10.0},
                                                                   .muF = 0.001,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 1e-3,
                                                                   .partDens = 800.0,
                                                                   .Rep=12.2500,
                                                                   .expectedDragForce = {424.2757}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0},
                                                                   .flowVel = {10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 0.01,
                                                                   .partDens = 800.0,
                                                                   .Rep=2450.0,
                                                                   .expectedDragForce = {4.8694}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0},
                                                                   .flowVel = {10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 0.00001,
                                                                   .partDens = 800.0,
                                                                   .Rep=2.450,
                                                                   .expectedDragForce = {1.4658e+05}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0},
                                                                   .flowVel = {10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 0.01,
                                                                   .partDens = 800.0,
                                                                   .Rep=2450.0,
                                                                   .expectedDragForce = {4.8694}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0},
                                                                   .flowVel = {10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 1e-3,
                                                                   .partDens = 500.0,
                                                                   .Rep=245.0,
                                                                   .expectedDragForce = {135.4623}},
                                         (DragModelTestParameters){.createDragModel = []() { return std::make_shared<ablate::particles::drag::Stokes>(); },
                                                                   .partVel = {0.0},
                                                                   .flowVel = {10.0},
                                                                   .muF = 5e-5,
                                                                   .rhoF = 1.225,
                                                                   .partDiam = 0.01,
                                                                   .partDens = 800.0,
                                                                   .Rep=245.0,
                                                                   .expectedDragForce = {8.4664}}));


}  // namespace ablateTesting::particles::drag
