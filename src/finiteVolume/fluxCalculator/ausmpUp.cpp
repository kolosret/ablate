#include "ausmpUp.hpp"

double ablate::finiteVolume::fluxCalculator::AusmpUp::totalTime = 0.0;
int ablate::finiteVolume::fluxCalculator::AusmpUp::callCount = 0;

ablate::finiteVolume::fluxCalculator::AusmpUp::AusmpUp(double mInf, std::shared_ptr<ablate::finiteVolume::processes::PressureGradientScaling> pgs) : pgs(pgs), mInf(mInf) {}

//ablate::finiteVolume::fluxCalculator::Direction ablate::finiteVolume::fluxCalculator::AusmpUp::AusmpUpFunction(void* ctx, PetscReal uL, PetscReal aL, PetscReal rhoL, PetscReal pL, PetscReal uR,
//                                                                                                                PetscReal aR, PetscReal rhoR, PetscReal pR, PetscReal* massFlux, PetscReal* p12) {
//
//    // extract pgs/minf if provided
//    auto ausmUp = (ablate::finiteVolume::fluxCalculator::AusmpUp*)ctx;
//    PetscReal pgsAlpha = ausmUp->pgs ? ausmUp->pgs->GetAlpha() : 1.0;
//    PetscReal mInf = ausmUp->mInf;
//
//    // Compute the density at the interface
//    PetscReal rho12 = (0.5) * (rhoL + rhoR);
//
//    // compute the speed of sound at a12
//    PetscReal a12 = 0.5 * (aL + aR) / pgsAlpha;  // Simple average of aL and aR.  This can be replaced with eq. 30;
//
//    // Compute the left and right mach numbers
//    PetscReal mL = uL / a12;
//    PetscReal mR = uR / a12;
//
//    // Compute mBar2 (eq 70)
//    PetscReal mBar2 = (PetscSqr(uL) + PetscSqr(uR)) / (2.0 * a12 * a12);
//
//    // compute mInf2 or set fa to unity
//    PetscReal fa = 1.0;
//    if (mInf > 0) {
//        PetscReal mInf2 = PetscSqr(mInf);
//
//        PetscReal mO2 = PetscMin(1.0, PetscMax(mBar2, mInf2));
//        PetscReal mO = PetscSqrtReal(mO2);
//        fa = mO * (2.0 - mO);
//    }
//    // compute the mach number on the interface
//    PetscReal m12 = M4Plus(mL) + M4Minus(mR) - (Kp / fa) * PetscMax(1.0 - (sigma * mBar2), 0) * (pR - pL) / (rho12 * a12 * a12 * pgsAlpha * pgsAlpha);
//
//    // store the mass flux;
//    Direction direction;
//    if (m12 > 0) {
//        direction = LEFT;
//        *massFlux = a12 * m12 * rhoL;
//    } else {
//        direction = RIGHT;
//        *massFlux = a12 * m12 * rhoR;
//    }
//
//    // Pressure
//    if (p12) {
//        double p5Plus = P5Plus(mL, fa);
//        double p5Minus = P5Minus(mR, fa);
//
//        *p12 = p5Plus * pL + p5Minus * pR - Ku * p5Plus * p5Minus * rho12 * fa * a12 * a12 * pgsAlpha * pgsAlpha * (mR - mL);
//        *p12 /= PetscSqr(pgsAlpha);
//    }
//    return direction;
//}
ablate::finiteVolume::fluxCalculator::Direction ablate::finiteVolume::fluxCalculator::AusmpUp::AusmpUpFunction(void* ctx, PetscReal uL, PetscReal aL, PetscReal rhoL, PetscReal pL, PetscReal uR,

                                                                                                               PetscReal aR, PetscReal rhoR, PetscReal pR, PetscReal* massFlux, PetscReal* p12) {


    //PAPI high level
//        int retval;
//        retval = PAPI_hl_region_begin("ausmup");

//        //PAPI low level
//        int EventSet = PAPI_NULL;
//        long long values[1];
//        // Create event set
//        PAPI_create_eventset(&EventSet);
//        PAPI_add_event(EventSet, PAPI_DP_OPS);
//        PAPI_start(EventSet);


    //Timing
//    double start = MPI_Wtime();


    auto ausmUp = (ablate::finiteVolume::fluxCalculator::AusmpUp*)ctx;  // 1 read (ctx), 1 write (ausmUp pointer assignment)
    PetscReal pgsAlpha = ausmUp->pgs ? ausmUp->pgs->GetAlpha() : 1.0;   // 1 read (ausmUp->pgs), 1 conditional read (GetAlpha), 1 write
    PetscReal mInf = ausmUp->mInf;                                      // 1 read, 1 write


    // Compute the density at the interface
    PetscReal rho12 = 0.5 * (rhoL + rhoR);  // 2 reads (rhoL, rhoR), 1 write

    // compute the speed of sound at a12
    PetscReal a12 = 0.5 * (aL + aR) / pgsAlpha;  // 3 reads (aL, aR, pgsAlpha), 1 write

    // Compute the left and right mach numbers
    PetscReal mL = uL / a12;  // 2 reads (uL, a12), 1 write
    PetscReal mR = uR / a12;  // 2 reads (uR, a12), 1 write

    //R 12,  W 7 so far
    // Compute mBar2 (eq 70)
    PetscReal mBar2 = (PetscSqr(uL) + PetscSqr(uR)) / (2.0 * a12 * a12);  // 3 reads (uL, uR, a12), 1 write

    // compute mInf2 or set fa to unity
    PetscReal fa = 1.0;                    // 1 write
    if (mInf > 0) {                        // 1 read (mInf)
        PetscReal mInf2 = PetscSqr(mInf);  // 1 read, 1 write

        PetscReal mO2 = PetscMin(1.0, PetscMax(mBar2, mInf2));  // 2 reads (mBar2, mInf2), 1 write
        PetscReal mO = PetscSqrtReal(mO2);                      // 1 read, 1 write
        fa = mO * (2.0 - mO);                                   // 1 read, 1 write
    }

    // compute the mach number on the interface
    PetscReal m12 = M4Plus(mL) + M4Minus(mR) - (Kp / fa) * PetscMax(1.0 - (sigma * mBar2), 0) * (pR - pL) / (rho12 * a12 * a12 * pgsAlpha * pgsAlpha);
    // ~12 reads (mL, mR, fa, mBar2, pR, pL, rho12, a12 x2, pgsAlpha x2, Kp, sigma), 1 write

    // store the mass flux;
    // doesnt matter which branch for memory
    Direction direction;               // 1 write
    if (m12 > 0) {                     // 1 read
        direction = LEFT;              // 1 write
        *massFlux = a12 * m12 * rhoL;  // 3 reads, 1 write
    } else {
        direction = RIGHT;             // 1 write
        *massFlux = a12 * m12 * rhoR;  // 3 reads, 1 write
    }

    //32 R + 2*M4_R, 13 W+2*M4_R so far + 1 if loop

    // Pressure
    if (p12) {                             // 1 read
        double p5Plus = P5Plus(mL, fa);    // 2 reads, 1 write
        double p5Minus = P5Minus(mR, fa);  // 2 reads, 1 write

        *p12 = p5Plus * pL + p5Minus * pR - Ku * p5Plus * p5Minus * rho12 * fa * a12 * a12 * pgsAlpha * pgsAlpha * (mR - mL);
        // ~11 reads (p5Plus, p5Minus, pL, pR, Ku, rho12, fa, a12 x2, pgsAlpha x2, mR, mL), 1 write
        *p12 /= PetscSqr(pgsAlpha);  // 2 reads (p12, pgsAlpha), 1 write
    }


//    PAPI low level
//        PAPI_stop(EventSet, values);
//        printf("FLOPs counted: %lld\n", values[0]);


    //PAPI high level
//    retval = PAPI_hl_region_end("ausmup");
//    if ( retval != PAPI_OK ){
//        std::cout << "abc" << std::endl;
//    }



    //Timing function
//    totalTime += MPI_Wtime() - start;
//    ++callCount;
//    if (callCount == 100000) {
//        PetscPrintf(
//            PETSC_COMM_WORLD, "StaticFunction called %d times, total time: %f s\n", ablate::finiteVolume::fluxCalculator::AusmpUp::callCount, ablate::finiteVolume::fluxCalculator::AusmpUp::totalTime);
//    }
    return direction;  // 1 read
}

PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M1Plus(PetscReal m) {
    return 0.5 * (m + PetscAbs(m)); // 2 reads (m), 1 write (return)
}

PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M2Plus(PetscReal m) {
    return 0.25 * PetscSqr(m + 1); // 1 read (m), 1 write (return)
}

PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M1Minus(PetscReal m) {
    return 0.5 * (m - PetscAbs(m)); // 2 reads (m), 1 write (return)
}

PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M2Minus(PetscReal m) {
    return -0.25 * PetscSqr(m - 1); // 1 read (m), 1 write (return)
}

PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M4Plus(PetscReal m) {
    if (PetscAbs(m) >= 1.0) { // 1 read (m)
        return M1Plus(m); // 1 read (m), 1 write (return)
    } else {
        return M2Plus(m) * (1.0 - 16.0 * beta * M2Minus(m));
        // 2 reads (m), 1 read (beta), 1 write (return)
    }
}

PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M4Minus(PetscReal m) {
    if (PetscAbs(m) >= 1.0) { // 1 read (m)
        return M1Minus(m); // 1 read (m), 1 write (return)
    } else {
        return M2Minus(m) * (1.0 + 16.0 * beta * M2Plus(m));
        // 2 reads (m), 1 read (beta), 1 write (return)
    }
}

PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::P5Plus(PetscReal m, double fa) {
    if (PetscAbs(m) >= 1.0) { // 1 read (m)
        return (M1Plus(m) / (m + 1E-30)); // 2 reads (m), 1 write (return)
    } else {
        double alpha = 3.0 / 16.0 * (-4.0 + 5 * fa * fa); // 2 reads (fa), 1 write (alpha)

        return (M2Plus(m) * ((2.0 - m) - 16. * alpha * m * M2Minus(m)));
        // 3 reads (m, alpha, fa), 1 write (return)
    }
}

PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::P5Minus(PetscReal m, double fa) {
    if (PetscAbs(m) >= 1.0) { // 1 read (m)
        return (M1Minus(m) / (m + 1E-30)); // 2 reads (m), 1 write (return)
    } else {
        double alpha = 3.0 / 16.0 * (-4.0 + 5 * fa * fa); // 2 reads (fa), 1 write (alpha)
        return (M2Minus(m) * ((-2.0 - m) + 16. * alpha * m * M2Plus(m)));
        // 3 reads (m, alpha, fa), 1 write (return)
    }
}

//PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M1Plus(PetscReal m) { return 0.5 * (m + PetscAbs(m)); }
//
//PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M2Plus(PetscReal m) { return 0.25 * PetscSqr(m + 1); }
//
//PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M1Minus(PetscReal m) { return 0.5 * (m - PetscAbs(m)); }
//PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M2Minus(PetscReal m) { return -0.25 * PetscSqr(m - 1); }
//
//PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M4Plus(PetscReal m) {
//    if (PetscAbs(m) >= 1.0) {
//        return M1Plus(m);
//    } else {
//        return M2Plus(m) * (1.0 - 16.0 * beta * M2Minus(m));
//    }
//}
//PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::M4Minus(PetscReal m) {
//    if (PetscAbs(m) >= 1.0) {
//        return M1Minus(m);
//    } else {
//        return M2Minus(m) * (1.0 + 16.0 * beta * M2Plus(m));
//    }
//}
//PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::P5Plus(PetscReal m, double fa) {
//    if (PetscAbs(m) >= 1.0) {
//        return (M1Plus(m) / (m + 1E-30));
//    } else {
//        // compute alpha
//        double alpha = 3.0 / 16.0 * (-4.0 + 5 * fa * fa);
//
//        return (M2Plus(m) * ((2.0 - m) - 16. * alpha * m * M2Minus(m)));
//    }
//}
//PetscReal ablate::finiteVolume::fluxCalculator::AusmpUp::P5Minus(PetscReal m, double fa) {
//    if (PetscAbs(m) >= 1.0) {
//        return (M1Minus(m) / (m + 1E-30));
//    } else {
//        double alpha = 3.0 / 16.0 * (-4.0 + 5 * fa * fa);
//        return (M2Minus(m) * ((-2.0 - m) + 16. * alpha * m * M2Plus(m)));
//    }
//}


#include "registrar.hpp"
REGISTER(ablate::finiteVolume::fluxCalculator::FluxCalculator, ablate::finiteVolume::fluxCalculator::AusmpUp, "A sequel to AUSM, Part II: AUSM+-up for all speeds, Meng-Sing Liou, Pages 137-170, 2006",
         OPT(double, "mInf", "the reference mach number"),
         OPT(ablate::finiteVolume::processes::PressureGradientScaling, "pgs", "Pressure gradient scaling is used to scale the acoustic propagation speed and increase time step for low speed flows"));