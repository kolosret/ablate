#include "sourceCalculatorZeroRK.hpp"
#include <math.h>
#include <algorithm>
#include "eos/zerork.hpp"
#include "finiteVolume/compressibleFlowFields.hpp"
#include "utilities/mpiUtilities.hpp"
#include "utilities/stringUtilities.hpp"

#include <unordered_map>
#include "domain/range.hpp"
#include <petscdmplex.h>
#include <petscfv.h>


void ablate::eos::zerorkeos::SourceCalculator::ChemistryConstraints::Set(const std::shared_ptr<ablate::parameters::Parameters>& options) {
    if (options) {
        verbose = options->Get("verbose", verbose);
        timinglog = options->Get("timingLog", timinglog);
        sparseJacobian = options->Get("sparseJacobian", sparseJacobian);
        relTolerance = options->Get("relTolerance", relTolerance);
        absTolerance = options->Get("absTolerance", absTolerance);
        thresholdTemperature = options->Get("thresholdTemperature", thresholdTemperature);
        stepLimiter = options->Get("steplimiter", stepLimiter);
        loadBalance = options->Get("loadBalance", loadBalance);
        useSEULEX = options->Get("useSEULEX", useSEULEX);
        iterative = options->Get("iterative", iterative);
        gpu = options->Get("gpu", gpu);
        maxiteration = options->Get("max_steps", maxiteration);
        reactorType = options->Get("reactorType", ReactorType::ConstantVolume);
        errorhandle = options->Get("errorhandle", errorhandle);
        dumpreactor = options->Get("dumpreactor", dumpreactor);
        dumpfailed = options->Get("dumpfailed", dumpfailed);
        n_reactors_max = options->Get("n_reactors_max", n_reactors_max);
        n_reactors_min = options->Get("n_reactors_min", n_reactors_min);
        cvode_num_retries = options->Get("cvode_num_retries", cvode_num_retries);
        cvode_retry_absolute_tolerance_adjustment = options->Get("cvode_retry_absolute_tolerance_adjustment", cvode_retry_absolute_tolerance_adjustment);
        cvode_retry_relative_tolerance_adjustment = options->Get("cvode_retry_relative_tolerance_adjustment", cvode_retry_relative_tolerance_adjustment);
    }
}

ablate::eos::zerorkeos::SourceCalculator::SourceCalculator(const std::vector<domain::Field>& fields, const std::shared_ptr<zerorkEOS> eosIn,
                                                           ablate::eos::zerorkeos::SourceCalculator::ChemistryConstraints constraints, const ablate::domain::Range& cellRange)
    : chemistryConstraints(constraints), eos(eosIn), numberSpecies(eosIn->GetSpeciesVariables().size()) {
    // determine the number of required cells
    std::size_t numberCells = cellRange.end - cellRange.start;

    // determine the source vector size
    sourceZeroRKAtI = std::vector<double>(numberCells * (eosIn->mech->getNumSpecies() + 1));

    auto chemEos = std::dynamic_pointer_cast<ablate::eos::ChemistryModel>(eos);

    auto speciesElementInformation = chemEos->GetSpeciesElementalInformation();
    auto elementInformation = chemEos->GetElementInformation();      // element -> atomic mass
    auto speciesMolecularMass = chemEos->GetSpeciesMolecularMass();    // species -> MW
    const auto& species = chemEos->GetSpeciesVariables();

    const std::vector<std::string> trackingElements{"C", "H"};
    std::map<std::string, double> massFractionsFuel;
    std::map<std::string, double> massFractionsOxidizer;
    for (const auto& spName : species) {
        massFractionsFuel[spName] = 0.0;
        massFractionsOxidizer[spName]  = 0.0;
    }

    const std::string fuelName = "MMA";
    const std::string oxidizerName  = "O2";

    if (!massFractionsFuel.count(fuelName) || !massFractionsOxidizer.count(oxidizerName)) {
        throw std::invalid_argument("SourceCalculator mixture fraction: could not find PMMA or O2 species in mechanism.");
    }

    massFractionsFuel[fuelName] = 1.0;
    massFractionsOxidizer[oxidizerName] = 1.0;

    zMixCoefficients.assign(species.size(), 0.0);
    for (std::size_t s = 0; s < species.size(); ++s) {
        const auto& spName = species[s];
        auto& speciesElement = speciesElementInformation[spName];
        for (const auto& elem : trackingElements) {
            zMixCoefficients[s] += speciesElement[elem] * elementInformation[elem] / speciesMolecularMass[spName];
        }
    }


    zMixFuel     = 0.0;
    zMixOxidizer = 0.0;
    for (std::size_t s = 0; s < species.size(); ++s) {
        const auto& spName = species[s];
        zMixFuel     += zMixCoefficients[s] * massFractionsFuel[spName];
        zMixOxidizer += zMixCoefficients[s] * massFractionsOxidizer[spName];
    }



    // Look for the euler field
    auto eulerField = std::find_if(fields.begin(), fields.end(), [](const auto& field) { return field.name == ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD; });
    if (eulerField == fields.end()) {
        throw std::invalid_argument("ablate::eos::zerorkEOS::BatchSource requires the ablate::finiteVolume::CompressibleFlowFields::EULER_FIELD Field");
    }
    eulerId = eulerField->id;

    auto densityYiField = std::find_if(fields.begin(), fields.end(), [](const auto& field) { return field.name == ablate::finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD; });
    if (densityYiField == fields.end()) {
        throw std::invalid_argument("ablate::eos::zerorkEOS::BatchSource requires the ablate::finiteVolume::CompressibleFlowFields::DENSITY_YI_FIELD Field");
    }
    densityYiId = densityYiField->id;

    int zerork_error_state = 0;
    zrm_handle = zerork_reactor_init();
    // load in mechanism for the plugin
    zerork_status_t zerom_status = zerork_reactor_set_mechanism_files(eos->reactionFile.c_str(), eos->thermoFile.c_str(), zrm_handle);
    if (zerom_status != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    zerork_status_t status_cvode = zerork_reactor_set_int_option("integrator", chemistryConstraints.useSEULEX, zrm_handle);
    if (status_cvode != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    // verbose 0 is no output, max level is 4
    zerork_status_t status_verbose = zerork_reactor_set_int_option("verbosity", chemistryConstraints.verbose, zrm_handle);
    if (status_verbose != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    zerork_status_t status_n_reacmax = zerork_reactor_set_int_option("n_reactors_max", chemistryConstraints.n_reactors_max, zrm_handle);
    if (status_n_reacmax != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    zerork_status_t status_n_reacmin = zerork_reactor_set_int_option("n_reactors_min", chemistryConstraints.n_reactors_min, zrm_handle);
    if (status_n_reacmin != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    // set to  0 to turn it off
    zerork_status_t status_loadbalance = zerork_reactor_set_int_option("load_balance", chemistryConstraints.loadBalance, zrm_handle);
    if (status_loadbalance != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    // Set tolerances
    zerork_status_t status_abstol = zerork_reactor_set_double_option("abs_tol", chemistryConstraints.absTolerance, zrm_handle);
    if (status_abstol != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;
    zerork_status_t status_reltol = zerork_reactor_set_double_option("rel_tol", chemistryConstraints.relTolerance, zrm_handle);
    if (status_reltol != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    switch (chemistryConstraints.reactorType) {
        case ReactorType::ConstantVolume: {
            zerork_status_t status_constvolume = zerork_reactor_set_int_option("constant volume", 1, zrm_handle);
            if (status_constvolume != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;
            break;
        }
        case ReactorType::ConstantPressure: {
            zerork_status_t status_constpress = zerork_reactor_set_int_option("constant_volume", 0, zrm_handle);
            if (status_constpress != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;
            break;
        }
    }
    zerork_status_t status_alwaysSolveTemp = zerork_reactor_set_int_option("always_solve_temperature", 1, zrm_handle);
    if (status_alwaysSolveTemp != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    // Kinetic rate limiter
    zerork_status_t status_steplimiter = zerork_reactor_set_double_option("step_limiter", chemistryConstraints.stepLimiter, zrm_handle);
    if (status_steplimiter != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    // Kinetic rate limiter
    zerork_status_t status_gpu = zerork_reactor_set_int_option("gpu", chemistryConstraints.gpu, zrm_handle);
    if (status_gpu != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    if (chemistryConstraints.timinglog) {
        zerork_reactor_set_string_option("reactor_timing_log_filename", "timing.log", zrm_handle);
    }

    // Use sparse matrix math for the jacobian
    if (!chemistryConstraints.sparseJacobian) {
        zerork_status_t status_sparse = zerork_reactor_set_int_option("dense", 1, zrm_handle);
        if (status_sparse != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;
    }
    zerork_status_t status_iterative = zerork_reactor_set_int_option("iterative", chemistryConstraints.iterative, zrm_handle);
    if (status_iterative != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    zerork_status_t status_maxsteps = zerork_reactor_set_int_option("max_steps", chemistryConstraints.maxiteration, zrm_handle);
    if (status_maxsteps != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    zerork_status_t status_numerrors = zerork_reactor_set_int_option("cvode_num_retries", chemistryConstraints.cvode_num_retries, zrm_handle);
    if (status_numerrors != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    zerork_status_t status_abserroradjust = zerork_reactor_set_double_option("cvode_retry_absolute_tolerance_adjustment", chemistryConstraints.cvode_retry_absolute_tolerance_adjustment, zrm_handle);
    if (status_abserroradjust != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    zerork_status_t status_relerroradjust = zerork_reactor_set_double_option("cvode_retry_relative_tolerance_adjustment", chemistryConstraints.cvode_retry_relative_tolerance_adjustment, zrm_handle);
    if (status_relerroradjust != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    if (chemistryConstraints.dumpreactor) {
        zerork_status_t dump = zerork_reactor_set_int_option("dump_reactors", 1, zrm_handle);
        if (dump != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;
    }

    zerork_status_t status_dumpfailedreactor = zerork_reactor_set_int_option("dump_failed_reactors", chemistryConstraints.dumpfailed, zrm_handle);
    if (status_dumpfailedreactor != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    zerork_status_t status_mech = zerork_reactor_load_mechanism(zrm_handle);  // make sure this call is after gpu setup
    if (status_mech != ZERORK_STATUS_SUCCESS) zerork_error_state += 1;

    if (zerork_error_state != 0) {
        throw std::invalid_argument("ablate::eos::zerork couldnt initialize, something is wrong...");
    }



}

void ablate::eos::zerorkeos::SourceCalculator::ComputeSource(const ablate::domain::Range& cellRange, PetscReal time, PetscReal dt, Vec globFlowVec) {
    StartEvent("zerorkEOS::SourceCalculator::ComputeSource");
    // Get the valid cell range over this region
    auto numberCells = cellRange.end - cellRange.start;

    // Get the solution dm
    DM solutionDm;
    VecGetDM(globFlowVec, &solutionDm) >> utilities::PetscUtilities::checkError;

    // get the rank
    PetscMPIInt rank;
    MPI_Comm_rank(PetscObjectComm((PetscObject)solutionDm), &rank) >> utilities::MpiUtilities::checkError;

    // get the flowSolution
    const PetscScalar* flowArray;
    VecGetArrayRead(globFlowVec, &flowArray) >> utilities::PetscUtilities::checkError;

    PetscInt dim;
    DMGetDimension(solutionDm, &dim) >> utilities::PetscUtilities::checkError;

    std::size_t nCells = cellRange.end - cellRange.start;
    reactorzMix.resize(nCells, 0.0);
    reactorzMixGrad.resize(nCells * dim, 0.0);

    // zerork state load up
    int nSpc = eos->mech->getNumSpecies();  // Number of Species
    int nState = nSpc + 1;
    int nReactors = numberCells;

    // Set up reactor initial states
    std::vector<double> reactorT(nReactors);
    std::vector<double> reactorP(nReactors);
    std::vector<double> density2(nReactors);
    std::vector<double> sensibleenergy(nReactors);
//    std::vector<double> velmag2(nReactors);
    std::vector<double> reactorMassFrac(nReactors * nSpc);
    std::vector<double> enthalpyOfFormation(nSpc);
    std::vector<int> reactorEval(nReactors, 0);

    std::vector<int> reactoriDs(nReactors);
    std::vector<int>    reactorLocalIdx(nReactors, -1);

    // Set up the vectors that are actually evaluated with the temperature threshold
    std::vector<double> reactorTEval(nReactors);
    std::vector<double> reactorPEval(nReactors);
    std::vector<double> reactorMassFracEval(nReactors * nSpc);

    // get the current state from petsc
//    std::vector<PetscInt> cellToLocal(nCells, -1);
    int p = 0; //active reactor index
    for (int i = cellRange.start; i < cellRange.end; ++i) {
        const PetscInt cell = cellRange.points ? cellRange.points[i] : i;
        const std::size_t k = i - cellRange.start;
//        cellToLocal[cell]   = k;

        const PetscScalar* eulerField = nullptr;
        DMPlexPointLocalFieldRead(solutionDm, cell, eulerId, flowArray, &eulerField) >> utilities::PetscUtilities::checkError;
        const PetscScalar* flowDensityField = nullptr;
        DMPlexPointLocalFieldRead(solutionDm, cell, densityYiId, flowArray, &flowDensityField) >> utilities::PetscUtilities::checkError;

        // get the current state at I
        auto density = eulerField[ablate::finiteVolume::CompressibleFlowFields::RHO];
        density2[k] = density;

        double yiSum = 0.0;
        for (int s = 0; s < nSpc - 1; s++) {
            reactorMassFrac[k * nSpc + s] = PetscMax(0.0, flowDensityField[s] / density);
            reactorMassFrac[k * nSpc + s] = PetscMin(1.0, reactorMassFrac[k * nSpc + s]);
            yiSum += reactorMassFrac[k * nSpc + s];
        }
        if (yiSum > 1.0) {
            for (PetscInt s = 0; s < nSpc - 1; s++) {
                // Limit the bounds
                reactorMassFrac[k * nSpc + s] /= yiSum;
            }
            reactorMassFrac[k * nSpc + nSpc - 1] = 0.0;
        } else {
            reactorMassFrac[k * nSpc + nSpc - 1] = 1.0 - yiSum;
        }

        // Compute the internal energy from total energy
        PetscReal speedSquare = 0.0;
        for (PetscInt d = 0; d < dim; d++) {
            speedSquare += PetscSqr(eulerField[ablate::finiteVolume::CompressibleFlowFields::RHOU + d] / density);
        }

        // compute the internal energy needed to compute temperature
        sensibleenergy[k] = eulerField[ablate::finiteVolume::CompressibleFlowFields::RHOE] / density - 0.5 * speedSquare;

        double enthalpymix = eos->mech->getMassEnthalpyFromTY(298.15, &reactorMassFrac[k * nSpc]);

        sensibleenergy[k] += enthalpymix;

        reactorT[k] = eos->mech->getTemperatureFromEY(sensibleenergy[k], &reactorMassFrac[k * nSpc], 2000);
        reactorP[k] = eos->mech->getPressureFromTVY(reactorT[k], 1 / density, &reactorMassFrac[k * nSpc]);

        reactorzMix[k]=ComputeMixtureFraction(&reactorMassFrac[k * nSpc]);
        // Set up the vector that is actually being solved
        if (reactorT[k] > chemistryConstraints.thresholdTemperature) {
            reactorEval[k] = 1;
            reactorTEval[p] = reactorT[k];
            reactorPEval[p] = reactorP[k];
            reactoriDs[p]=cell;
            reactorLocalIdx[p]  = static_cast<int>(k);
            for (int s = 0; s < nSpc; s++) {
                reactorMassFracEval[p * nSpc + s] = reactorMassFrac[k * nSpc + s];
            }
//            reactorzMixEval[p]=ComputeMixtureFraction(&reactorMassFracEval[p * nSpc]);
            p += 1;
        }
    }

    ComputeMixtureFractionGradients(solutionDm, cellRange, reactorzMix, reactorzMixGrad, dim);
//    double maxelem = *std::max_element(reactorzMixGrad.begin(),reactorzMixGrad.end());
//    std::cout << "The max gradient is: "<< maxelem <<std::endl;
//    double minelem = *std::min_element(reactorzMixGrad.begin(),reactorzMixGrad.end());
//    std::cout << "The min gradient is: "<< minelem <<std::endl;
    double diff=0;

//    if (!reactorzMixGrad.empty()) {
//        double maxelem = *std::max_element(reactorzMixGrad.begin(),reactorzMixGrad.end());
//        double minelem = *std::min_element(reactorzMixGrad.begin(),reactorzMixGrad.end());
//        std::cout << "The max gradient is: "<< maxelem <<std::endl;
//        std::cout << "The min gradient is: "<< minelem <<std::endl;
//    } else {
//        std::cout << "No gradients on this rank.\n";
//    }


    //Constants for method 4
    const double p2 = 0.01560047;
    const double p1 = 0.30251833;
    const double p0 = 1.19877839;

    //Constants for method 5
    const double pt2 = 0.01363758;
    const double pt1 = 0.20908996;
    const double pt0 = -2.98799396;

    //Constants for method 6 for CPU
//    const double pt2_max = 0.01097835677406711;
//    const double pt1_max = 0.20394039678277337;
//    const double pt0_max = -2.7264374443298514;

    //Constants for method 6 for GPU
    const double pt2_max = 0.010940142202676377;
    const double pt1_max = 0.17848686803192607;
    const double pt0_max = -4.295246112850871;

    char* weight = getenv("SET_WEIGHT");

    std::vector<double> reactorzMixGradMag(p, 0.0);
    std::vector<double> reactorChiEval(p, -1.0);

    int OHind=12; //Eventually this shouldnt be hardcoded
    double OHcutoff=0.001;
    double Xin_lo=-3.;
    double Xin_high=1.;
    double Xout_low=0.0;
    double Xout_high=1.0;

    for (int i = 0; i<p;++i){

//        const PetscInt cell  = reactoriDs[i];        // DMPlex cell point
//        const PetscInt k     = cellToLocal[cell];
        const int      k    = reactorLocalIdx[i];

        if (k < 0 || k >= static_cast<int>(nCells)) {
            throw std::runtime_error("ComputeSource: invalid local cell index for reactor.");
        }
        diff = SutherlandDiff(reactorTEval[i], density2[k]);
        //Calculate the square of the gradient
//        diff= SutherlandDiff(reactorTEval[i],density2[reactoriDs[i]]);
//        for (int s = 0; s < dim; s++) {
//            reactorzMixGradMag[i] += reactorzMixGrad[i * dim + s]*reactorzMixGrad[i * dim + s];
//        }
        for (int s = 0; s < dim; s++) {
            double g = reactorzMixGrad[static_cast<std::size_t>(k) * dim + s];
            reactorzMixGradMag[i] += g * g;
        }

        if (weight && atoi(weight) == 1) {
            reactorChiEval[i] = 2 * diff * reactorzMixGradMag[i];
        }else if (weight && atoi(weight) == 2) {
            double chi = 2 * diff * reactorzMixGradMag[i];
            reactorChiEval[i] = chi*chi;
        }else if (weight && atoi(weight) == 3) {
            reactorChiEval[i] = 2 * diff * reactorzMixGradMag[i];
        }else if (weight && atoi(weight) == 4) {
            double chi =2 * diff * reactorzMixGradMag[i];
            chi=std::max(chi,0.0001);
            //Using the fitting constants
            double u = std::log10(chi);
            double v = p0 + p1 * u + p2 * u * u;
            reactorChiEval[i] = std::pow(10.0, v);
        }else if (weight && atoi(weight) == 5) {
            double chi =2 * diff * reactorzMixGradMag[i];
            chi=std::max(chi,0.0001);
            //Using the fitting constants
            double u = std::log10(chi);
            double v = pt0 + pt1 * u + pt2 * u * u;
            reactorChiEval[i] = std::pow(10.0, v);
        }else if (weight && atoi(weight) == 6) {

            double chi = 2 * diff * reactorzMixGradMag[i];
            //limit chi
            chi=std::max(chi,0.0001);

            double Yi_OH = reactorMassFracEval[i * nSpc + OHind];
            //Using the fitting constants
            double u = std::log10(chi);
            double logchi = u;
                if (u >= Xin_lo && u <= Xin_high && Yi_OH > OHcutoff){
                logchi= Xout_low + (u - Xin_lo)*(Xout_high - Xout_low)/(Xin_high - Xin_lo);
//                std::cout<< "Chi was adjusted from: "<<u<<" to "<<  logchi <<" and YiOH is "<< Yi_OH << std::endl;
            }

            double v = pt0_max + pt1_max * logchi + pt2_max * logchi * logchi;
            reactorChiEval[i] = std::pow(10.0, v);
        }else if (weight && atoi(weight) == 7) {
            double chi =2 * diff * reactorzMixGradMag[i];
            chi=std::max(chi,0.0001);
            //Using the fitting constants
            double u = std::log10(chi);
            double v = pt0 + pt1 * u + pt2 * u * u;
            reactorChiEval[i] = std::pow(10.0, v);
        }else if (weight && atoi(weight) == 8) {

            double chi = 2 * diff * reactorzMixGradMag[i];
            //limit chi
            chi=std::max(chi,0.0001);

            double Yi_OH = reactorMassFracEval[i * nSpc + OHind];
            //Using the fitting constants
            double u = std::log10(chi);
            double logchi = u;
            if (u >= Xin_lo && u <= Xin_high && Yi_OH > OHcutoff){
                logchi= Xout_low + (u - Xin_lo)*(Xout_high - Xout_low)/(Xin_high - Xin_lo);
                //                std::cout<< "Chi was adjusted from: "<<u<<" to "<<  logchi <<" and YiOH is "<< Yi_OH << std::endl;
            }

            double v = pt0_max + pt1_max * logchi + pt2_max * logchi * logchi;
            reactorChiEval[i] = std::pow(10.0, v);
        }else if (weight && atoi(weight) == 9) {
            double chi =2 * diff * reactorzMixGradMag[i];
            chi=std::max(chi,0.0001);
            //Using the fitting constants
            double u = std::log10(chi);
            double v = pt0 + pt1 * u + pt2 * u * u;
            reactorChiEval[i] = std::pow(10.0, v);
        }
        else{// Dont set anything by default
            reactorChiEval[i] = 2 * diff * reactorzMixGradMag[i];
        }

//        std::cout << "The Chi for cell: " << cell <<" are: " << reactorChiEval[i] <<std::endl;

    }

    if (p > 0) {
//        double maxChi = *std::max_element(reactorChiEval.begin(),reactorChiEval.end());
//        double minChi = *std::min_element(reactorChiEval.begin(),reactorChiEval.end());
//        std::cout << "The max chi is: "<< maxChi <<std::endl;
//        std::cout << "The min chi is: "<< minChi <<std::endl;
    } else {
        std::cout << "No active reactors (p == 0) on this rank; Chi not computed.\n";
    }

//    double maxChi = *std::max_element(reactorChiEval.begin(),reactorChiEval.end());
//    std::cout << "The max chi is: "<< maxChi <<std::endl;
//    double minChi = *std::min_element(reactorChiEval.begin(),reactorChiEval.end());
//    std::cout << "The min chi is: "<< minChi <<std::endl;



    int N = reactorChiEval.size();
    std::vector<double> reactorWeights(N);

    // save Yi for source terms
    std::vector<double> ys = reactorMassFracEval;

    if (weight && atoi(weight) == 1) {
        std::cout <<"Setting Chi as the weight"<<std::endl;
        if (N > 0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorChiEval[0], zrm_handle);
        }
    }
    else if(weight && atoi(weight) == 2){
        std::cout <<"Using option 2 for the weights"<<std::endl;
        if (N > 0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorChiEval[0], zrm_handle);
        }
    }else if(weight && atoi(weight) == 3){
        std::cout <<"Using option 3 for the weights"<<std::endl;
        reactorWeights = reactorChiEval;

        for(int i = 0; i < N; ++i){
            reactorWeights[i] = std::max(reactorChiEval[i], 0.0);
        }

        // Build index list
        std::vector<int> idx(N);
        std::iota(idx.begin(), idx.end(), 0);

        // Sort indices by descending Chi
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            return reactorChiEval[a] > reactorChiEval[b];
        });

        // Max Chi value
        double chi_max = reactorChiEval[idx[0]];

        // Assign MAX Chi to top-N hardest (batch size)
        int n = std::min(chemistryConstraints.n_reactors_max, N);
        for(int i = 0; i < n; ++i){
            int id = idx[i];
            reactorWeights[id] = chi_max;
        }
        if (N > 0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorWeights[0], zrm_handle);
        }

    }else if(weight && atoi(weight) == 4){
        std::cout <<"Using option 4 for the weights"<<std::endl;
        if (N > 0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorChiEval[0], zrm_handle);
        }
    }else if(weight && atoi(weight) == 5){
        std::cout <<"Using option 5 for the weights"<<std::endl;
        if (N > 0) {
            double minVal = *std::min_element(reactorChiEval.begin(), reactorChiEval.end());
            minVal = std::max(minVal,0.00001);
            for (auto& v : reactorChiEval) {
                v /= minVal;
            }
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorChiEval[0], zrm_handle);
        }
    }else if(weight && atoi(weight) == 6){
        std::cout <<"Using option 6 for the weights"<<std::endl;
        if (N > 0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorChiEval[0], zrm_handle);
        }
    }else if(weight && atoi(weight) == 7){
        std::cout <<"Using option 6 for the weights"<<std::endl;
        if (N > 0) {
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorChiEval[0], zrm_handle);
        }
    }else if(weight && atoi(weight) == 8){
        std::cout <<"Using option 5 for the weights"<<std::endl;
        if (N > 0) {
            double minVal = *std::min_element(reactorChiEval.begin(), reactorChiEval.end());
            minVal = std::max(minVal,0.00001);
            for (auto& v : reactorChiEval) {
                v /= minVal;
            }
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorChiEval[0], zrm_handle);
        }
    }else if(weight && atoi(weight) == 9){
        std::cout <<"Using option 5 for the weights"<<std::endl;
        if (N > 0) {
            double minVal = *std::min_element(reactorChiEval.begin(), reactorChiEval.end());
            minVal = std::max(minVal,0.00001);
            for (auto& v : reactorChiEval) {
                v /= minVal;
            }
            zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorChiEval[0], zrm_handle);
        }
    }
    else{
        std::cout <<"No weigths are set"<<std::endl;
    }

//    zerork_reactor_set_aux_field_pointer(ZERORK_FIELD_COST, &reactorChiEval[0], zrm_handle);

    // Solve for all the reactors, this calls zerork_cfd_plugin.cpp, zerork_reactor_manager.cpp
    auto nReactorsEval = std::reduce(reactorEval.begin(), reactorEval.end());
    zerork_status_t flag = ZERORK_STATUS_SUCCESS;
    flag = zerork_reactor_set_reactor_ids(&reactoriDs[0],zrm_handle);

    flag = zerork_reactor_solve(1, time, dt, nReactorsEval, &reactorTEval[0], &reactorPEval[0], &reactorMassFracEval[0], zrm_handle);

    if (flag != ZERORK_STATUS_SUCCESS) {
        std::cout << "Integration failed on some of the ranks, even after reducing tolerances in ZeroRK."
                  << "\n";
        if (chemistryConstraints.errorhandle == 1) {
            int ii = 0;
            // For now try to manually decrease the tolerances and recalculate every rank!
            // Zerork already reduced the tolerances
            while (flag != ZERORK_STATUS_SUCCESS) {
                ++ii;
                std::cout << "Manually tightening tolerances further."
                          << "\n";
                zerork_reactor_set_double_option(
                    "abs_tol", chemistryConstraints.absTolerance * pow(chemistryConstraints.cvode_retry_absolute_tolerance_adjustment, ii * chemistryConstraints.cvode_num_retries), zrm_handle);
                flag = zerork_reactor_solve(2, time, dt, nReactorsEval, &reactorTEval[0], &reactorPEval[0], &reactorMassFracEval[0], zrm_handle);
                // Try tightening the tolerances
                if (ii == 2) {
                    std::cout << "At this point the tolerances are probably too tight."
                              << "\n"
                              << "Consider dumping the state, by setting dumpfailed = 1 in the input file and try to understand why is it failing."
                              << "\n";
                    break;
                }
            }
            // Integration error usually only occur for certain specific states, which will hopefully be advected away for the next step...
            // Resetting the tolerances to the original inputs
            zerork_reactor_set_double_option("abs_tol", chemistryConstraints.absTolerance, zrm_handle);
            zerork_reactor_set_double_option("rel_tol", chemistryConstraints.relTolerance, zrm_handle);
        }
        if (chemistryConstraints.errorhandle == 2) {
            try {
                // For errorhandle 2 stop the simualtion
                if (flag != ZERORK_STATUS_SUCCESS) {
                    std::cout << "Warning: Could not integrate chemistry after reducing the tolerances multiple times."
                              << "\n";
                    std::cout << "Option 2 was selected for error handling, the simulations exits now. "
                              << "\n";
                    throw std::runtime_error("ablate::eos::zerorkEOS::Computesource zerork couldn't integrate the simulation.");
                }
            } catch (const runtime_error& e) {
                exit(1);
            }
        }
    }

    // Set all the sources to 0
    sourceZeroRKAtI.assign(nState * nReactors, 0);

    for (int s = 0; s < nSpc - 1; s++) {
        std::vector<double> tempvec(nSpc, 0.);
        tempvec[s] = 1;
        enthalpyOfFormation[s] = eos->mech->getMassEnthalpyFromTY(298.15, &tempvec[0]);
    }

    // Here we should recompute density for constant pressure reactors for tighter coupling
    // however this can get to unexpected non-physical results.
    // Think about using source terms for density for different coupling...
    //    if(chemistryConstraints.reactorType==ReactorType::ConstantPressure){
    //        int q=0;
    //        for (int i=0;i<nReactors;i++){
    //            if (reactorEval[i]==1){
    //                density2[i]=eos->mech->getDensityFromTPY(reactorTEval[q], reactorPEval[q],&reactorMassFracEval[q]);
    //                q+=1;
    //            }
    //        }
    //    }

    int q = 0;
    for (int i = 0; i < nReactors; ++i) {
        if (reactorEval[i] != 1) {
            for (int j = 0; j < nState; ++j) {
                // Set sourceterms to 0 if the reactor was evaluated
                sourceZeroRKAtI[i * nState + j] = 0;
            }
            q += 1;
        } else {
            for (int s = 0; s < nSpc; s++) {
                sourceZeroRKAtI[i * nState] += (ys[(i - q) * nSpc + s] - reactorMassFracEval[(i - q) * nSpc + s]) * enthalpyOfFormation[s];
            }

            for (int s = 0; s < nSpc; ++s) {
                // for constant density problem, d Yi rho/dt = rho * d Yi/dt + Yi*d rho/dt = rho*dYi/dt ~~ rho*(Yi+1 - Y1)/dt
                sourceZeroRKAtI[i * nState + s + 1] = reactorMassFracEval[(i - q) * nSpc + s] - ys[(i - q) * nSpc + s];
            }

            // Now scale everything by density/dt
            for (int j = 0; j < nState; ++j) {
                // for constant density problem, d Yi rho/dt = rho * d Yi/dt + Yi*d rho/dt = rho*dYi/dt ~~ rho*(Yi+1 - Y1)/dt
                sourceZeroRKAtI[i * nState + j] *= density2[i] / dt;
            }
        }
    }
    VecRestoreArrayRead(globFlowVec, &flowArray) >> utilities::PetscUtilities::checkError;
    EndEvent();
}
void ablate::eos::zerorkeos::SourceCalculator::AddSource(const ablate::domain::Range& cellRange, Vec, Vec locFVec) {
    StartEvent("zerorkEOS::SourceCalculator::AddSource");

    // get access to the fArray
    PetscScalar* fArray;
    VecGetArray(locFVec, &fArray) >> utilities::PetscUtilities::checkError;

    // Get the solution dm
    DM dm;
    VecGetDM(locFVec, &dm) >> utilities::PetscUtilities::checkError;

    int nSpc = eos->mech->getNumSpecies();
    int nState = nSpc + 1;
    for (int i = cellRange.start; i < cellRange.end; ++i) {
        const PetscInt cell = cellRange.points ? cellRange.points[i] : i;
        const std::size_t k = i - cellRange.start;

        // Get the current state variables for this cell
        PetscScalar* eulerSource = nullptr;
        DMPlexPointLocalFieldRef(dm, cell, eulerId, fArray, &eulerSource) >> utilities::PetscUtilities::checkError;
        PetscScalar* densityYiSource = nullptr;
        DMPlexPointLocalFieldRef(dm, cell, densityYiId, fArray, &densityYiSource) >> utilities::PetscUtilities::checkError;

        eulerSource[ablate::finiteVolume::CompressibleFlowFields::RHOE] += sourceZeroRKAtI[k * nState];
        for (std::size_t sp = 0; sp < numberSpecies; sp++) {
            densityYiSource[sp] += sourceZeroRKAtI[k * nState + sp + 1];
        }
    }

    // cleanup
    VecRestoreArray(locFVec, &fArray) >> utilities::PetscUtilities::checkError;
    EndEvent();
}

double ablate::eos::zerorkeos::SourceCalculator::ComputeMixtureFraction(const double* yi) const {
    double zMix = 0.0;
    for (std::size_t s = 0; s < zMixCoefficients.size(); ++s) {
        zMix += zMixCoefficients[s] * yi[s];
    }
    if (zMix < 0.0) zMix = 0.0;
    if (zMix > 1.0) zMix = 1.0;
    return (zMix - zMixOxidizer) / (zMixFuel - zMixOxidizer);
}

std::ostream& ablate::eos::zerorkeos::operator<<(std::ostream& os, const ablate::eos::zerorkeos::SourceCalculator::ReactorType& v) {
    switch (v) {
        case ablate::eos::zerorkeos::SourceCalculator::ReactorType::ConstantPressure:
            return os << "ConstantPressure";
        case ablate::eos::zerorkeos::SourceCalculator::ReactorType::ConstantVolume:
            return os << "ConstantVolume";
        default:
            return os;
    }
}

std::istream& ablate::eos::zerorkeos::operator>>(std::istream& is, ablate::eos::zerorkeos::SourceCalculator::ReactorType& v) {
    std::string enumString;
    is >> enumString;

    // make the comparisons easier to converting to lower
    ablate::utilities::StringUtilities::ToLower(enumString);

    if (enumString == "constantvolume") {
        v = ablate::eos::zerorkeos::SourceCalculator::ReactorType::ConstantVolume;
    } else if (enumString == "constantpressure") {
        // default to constant pressure
        v = ablate::eos::zerorkeos::SourceCalculator::ReactorType::ConstantPressure;
    } else {
        throw std::invalid_argument(
            " Unknown reactor type set. \n"
            " Acceptable reactor types: ConstantPressure, ConstantVolume. \n"
            " Default is Contstant volume");
    }
    return is;
}


double ablate::eos::zerorkeos::SourceCalculator::SutherlandDiff(double T,double rho){
    double Diff=0;
    double mu0=1.716e-5;  // Pa·s at T0 (air)
    double T0=273.15;     // K
    double S=110.4;       // Sutherland constant for air [K]
    double Sc=0.7;       // Schmidt number (approx for air)

    double mu = mu0 * std::pow((T / T0), 1.5) * (T0 + S) / (T + S);
    return Diff = mu / (rho * Sc);
}









void ablate::eos::zerorkeos::SourceCalculator::ComputeMixtureFractionGradients(
    DM dm,
    const ablate::domain::Range& cellRange,
    const std::vector<double>& zMixValues,
    std::vector<double>& zMixGrad,
    PetscInt dim) {

    // number of local cells in this range
    const PetscInt nCells = cellRange.end - cellRange.start;
    if ((PetscInt)zMixValues.size() != nCells) {
        throw std::runtime_error("ComputeMixtureFractionGradients: zMixValues size must equal number of cells in cellRange.");
    }

    // get spatial dimension if dim <= 0
    PetscInt dmDim;
    DMGetDimension(dm, &dmDim) >> utilities::PetscUtilities::checkError;
    if (dim <= 0) dim = dmDim;

    // geometry: cell and face
    Vec cellGeomVec = nullptr, faceGeomVec = nullptr;
    DMPlexComputeGeometryFVM(dm, &cellGeomVec, &faceGeomVec) >> utilities::PetscUtilities::checkError;

    DM dmFace = nullptr, dmCell = nullptr;
    const PetscScalar* faceGeomArray = nullptr;
    const PetscScalar* cellGeomArray = nullptr;
    VecGetDM(faceGeomVec, &dmFace) >> utilities::PetscUtilities::checkError;
    VecGetDM(cellGeomVec, &dmCell) >> utilities::PetscUtilities::checkError;
    VecGetArrayRead(faceGeomVec, &faceGeomArray) >> utilities::PetscUtilities::checkError;
    VecGetArrayRead(cellGeomVec, &cellGeomArray) >> utilities::PetscUtilities::checkError;

    // ghost label (for faces)
    DMLabel ghostLabel = nullptr;
    DMGetLabel(dm, "ghost", &ghostLabel) >> utilities::PetscUtilities::checkError;

    // build map from DMPLEX cell point -> local index [0..nCells)
    PetscInt cStart, cEnd;
    DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd) >> utilities::PetscUtilities::checkError;
    std::vector<PetscInt> cellToLocal(cEnd, -1);
    for (PetscInt i = cellRange.start; i < cellRange.end; ++i) {
        const PetscInt cell = cellRange.points ? cellRange.points[i] : i;
        const PetscInt localIdx = i - cellRange.start;
        cellToLocal[cell] = localIdx;
    }

    // zero output gradients
    zMixGrad.assign((std::size_t)nCells * dim, 0.0);

    // loop over cells in the range
    for (PetscInt i = cellRange.start; i < cellRange.end; ++i) {
        const PetscInt cell     = cellRange.points ? cellRange.points[i] : i;
        const PetscInt localIdx = i - cellRange.start;

        // cell geometry
        const PetscFVCellGeom* cg = nullptr;
        DMPlexPointLocalRead(dmCell, cell, cellGeomArray, &cg) >> utilities::PetscUtilities::checkError;
        if (!cg || cg->volume <= 0.0) continue;

        const double zCell = zMixValues[(std::size_t)localIdx];

        // accumulator for ∑ z_f * n_f A_f
        double gradAccum[3] = {0.0, 0.0, 0.0};

        // faces of this cell
        PetscInt        numFaces = 0;
        const PetscInt* faces    = nullptr;
        DMPlexGetConeSize(dm, cell, &numFaces) >> utilities::PetscUtilities::checkError;
        DMPlexGetCone(dm, cell, &faces) >> utilities::PetscUtilities::checkError;

        for (PetscInt fi = 0; fi < numFaces; ++fi) {
            const PetscInt face = faces[fi];

            // skip ghost / boundary / children
            PetscInt  ghost = -1;
            PetscBool boundary;
            PetscInt  numChildren;
            if (ghostLabel) {
                DMLabelGetValue(ghostLabel, face, &ghost) >> utilities::PetscUtilities::checkError;
            }
            DMIsBoundaryPoint(dm, face, &boundary) >> utilities::PetscUtilities::checkError;
            DMPlexGetTreeChildren(dm, face, &numChildren, nullptr) >> utilities::PetscUtilities::checkError;
            if (ghost >= 0 || boundary || numChildren > 0) continue;

            // get neighbors of this face
            const PetscInt* cells = nullptr;
            PetscInt        supportSize = 0;
            DMPlexGetSupportSize(dm, face, &supportSize) >> utilities::PetscUtilities::checkError;
            DMPlexGetSupport(dm, face, &cells) >> utilities::PetscUtilities::checkError;
            if (supportSize == 0) continue; // should not happen

            const PetscInt c0 = cells[0];
            const PetscInt c1 = (supportSize > 1) ? cells[1] : -1;

            // orientation of normal: use sign to orient fg->normal wrt this cell
            double sign = 0.0;
            if (cell == c0) {
                sign = 1.0;
            } else if (cell == c1) {
                sign = -1.0;
            } else {
                continue; // this face is not attached to this cell (paranoia)
            }

            // face geometry
            const PetscFVFaceGeom* fg = nullptr;
            DMPlexPointLocalRead(dmFace, face, faceGeomArray, &fg) >> utilities::PetscUtilities::checkError;
            if (!fg) continue;

            // mixture fraction at face: central average if neighbor exists in our range
            double zFace = zCell;
            if (c1 >= 0) {
                const PetscInt otherCell = (cell == c0) ? c1 : c0;
                const PetscInt otherLocal = (otherCell >= 0 && otherCell < cEnd) ? cellToLocal[otherCell] : -1;
                if (otherLocal >= 0) {
                    const double zOther = zMixValues[(std::size_t)otherLocal];
                    zFace = 0.5 * (zCell + zOther);
                }
            }

            // accumulate: z_f * (outward normal * area)
            for (PetscInt d = 0; d < dim; ++d) {
                gradAccum[d] += zFace * sign * static_cast<double>(fg->normal[d]);
            }
        }

        const double invVol = 1.0 / static_cast<double>(cg->volume);
        for (PetscInt d = 0; d < dim; ++d) {
            zMixGrad[(std::size_t)localIdx * dim + d] = gradAccum[d] * invVol;
        }
    }

    // cleanup geometry vecs
    VecRestoreArrayRead(faceGeomVec, &faceGeomArray) >> utilities::PetscUtilities::checkError;
    VecRestoreArrayRead(cellGeomVec, &cellGeomArray) >> utilities::PetscUtilities::checkError;
    VecDestroy(&faceGeomVec) >> utilities::PetscUtilities::checkError;
    VecDestroy(&cellGeomVec) >> utilities::PetscUtilities::checkError;
}
