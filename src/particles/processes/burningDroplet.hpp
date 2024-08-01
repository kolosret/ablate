#ifndef ABLATELIBRARY_BURNINGDROPLET_HPP
#define ABLATELIBRARY_BURNINGDROPLET_HPP

#include <string>
#include "coupledProcess.hpp"
#include "particles/processes/burnmodel/burnModel.hpp"
#include "particles/processes/drag/dragModel.hpp"
#include "petsc.h"


namespace ablate::particles::processes {

class BurningDroplet : public CoupledProcess {
   private:
    //! store the coupled field name
    const std::string coupledFieldName;

    std::vector<double> eulerlocal;
    std::vector<double> densityYilocal;


    //! the function for the field, it should be the same size as the field.  This function is integrated in time
    const std::shared_ptr<mathFunctions::MathFunction> sourceFunction;

    // this class has a combustion model
    std::shared_ptr<particles::processes::burnmodel::BurnModel> burnprocess;

    // this class has to include a drag model
    std::shared_ptr<particles::processes::drag::DragModel> dragprocess;


    void stepaux(PetscReal startTime, PetscReal endTime,  accessors::SwarmAccessor& swarmAccessorPreStep,
                 accessors::SwarmAccessor& swarmAccessorPostStep,
                 ablate::particles::accessors::EulerianAccessor& eulerianAccessor,
                 accessors::EulerianSourceAccessor& eulerianSourceAccessor);

   public:
    /**
     * Adds an arbitrary source function for each particle to the eulerian field
     * @param coupledFieldName the name of the eulerian coupled field
     * @param sourceFunction the function to compute the source
     */
    BurningDroplet(std::string coupledFieldName, std::shared_ptr<mathFunctions::MathFunction> sourceFunction, const std::shared_ptr<parameters::Parameters>& parameters);

    /**
     * There is no RHS function for the ArbitraryEulerianSource
     */
    void ComputeRHS(PetscReal time, accessors::SwarmAccessor& swarmAccessor, accessors::RhsAccessor& rhsAccessor, accessors::EulerianAccessor& eulerianAccessor) override {}

    /**
     * Add the arbitrary source to the eulerianSourceAccessor
     * @param startTime
     * @param endTime
     * @param swarmAccessorPreStep
     * @param swarmAccessorPostStep
     * @param eulerianSourceAccessor
     */
    void ComputeEulerianSource(PetscReal startTime, PetscReal endTime, accessors::SwarmAccessor& swarmAccessorPreStep, accessors::SwarmAccessor& swarmAccessorPostStep,
                               accessors::EulerianSourceAccessor& eulerianSourceAccessor) override;

    void densityYiSource(PetscReal startTime, PetscReal endTime, ablate::particles::accessors::SwarmAccessor& swarmAccessorPreStep,
                                                                       ablate::particles::accessors::SwarmAccessor& swarmAccessorPostStep,
                                                                       ablate::particles::accessors::EulerianSourceAccessor& eulerianSourceAccessor);


};

}  // namespace ablate::particles::processes
#endif
