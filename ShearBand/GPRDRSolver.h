#ifndef __GPRDRSolver_CLASS_HEADER__
#define __GPRDRSolver_CLASS_HEADER__

// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================

#include <string>

#include "exahype/solvers/LimitingADERDGSolver.h"
#include "GPRDRSolver_ADERDG.h"
#include "GPRDRSolver_FV.h"

namespace GPRDR{
  class GPRDRSolver;
}

class GPRDR::GPRDRSolver: public exahype::solvers::LimitingADERDGSolver {  
  public:
    static constexpr int NumberOfVariables         = GPRDR::AbstractGPRDRSolver_ADERDG::NumberOfVariables;
    static constexpr int NumberOfParameters        = GPRDR::AbstractGPRDRSolver_ADERDG::NumberOfParameters;
    static constexpr int Order                     = GPRDR::AbstractGPRDRSolver_ADERDG::Order;
    static constexpr int NumberOfGlobalObservables = GPRDR::AbstractGPRDRSolver_ADERDG::NumberOfGlobalObservables;
    static constexpr int NumberOfDMPObservables    = GPRDR::AbstractGPRDRSolver_ADERDG::NumberOfDMPObservables;
    static constexpr int PatchSize                 = GPRDR::AbstractGPRDRSolver_FV::PatchSize;
    static constexpr int GhostLayerWidth           = GPRDR::AbstractGPRDRSolver_FV::GhostLayerWidth;
      
    // limiter projection matrices
    double dg2fv[(Order+1)*PatchSize];
    double fv2dg[(Order+1)*PatchSize];
    double leg2lob[(Order+1)*(Order+1)];
      
    GPRDRSolver(
        const double maximumMeshSize,
        const int maximumMeshDepth,
        const int haloCells,
        const int haloBufferCells,
        const int limiterBufferCells,
        const int regularisedFineGridLevels,
        const exahype::solvers::Solver::TimeStepping timeStepping,
        const int DMPObservables,
        const double DMPRelaxationParameter,
        const double DMPDifferenceScaling
);
    
    void projectOnFVLimiterSpace(const double* const luh, double* const lim) const override;
    void projectOnDGSpace(const double* const lim, double* const luh) const override;
    bool discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, double* const boundaryMinPerVariables, double* const boundaryMaxPerVariables) override;
    void findCellLocalMinAndMax(const double* const luh, double* const localMinPerVariables, double* const localMaxPerVariable) override;
    void findCellLocalLimiterMinAndMax(const double* const lim, double* const localMinPerObservable, double* const localMaxPerObservable) override;
};

#endif // __GPRDRSolver_CLASS_HEADER__