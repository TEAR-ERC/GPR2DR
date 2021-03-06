// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
// ========================
//   www.exahype.eu
// ========================
#ifndef POSTPROCESSING_ProbeWriter4_CLASS_HEADER_
#define POSTPROCESSING_ProbeWriter4_CLASS_HEADER_

#include "exahype/plotters/Plotter.h"

namespace GPRDR {
  class GPRDRSolver;
  class ProbeWriter4;
}

class GPRDR::ProbeWriter4 : public exahype::plotters::Plotter::UserOnTheFlyPostProcessing {
public:
  ProbeWriter4(GPRDR::GPRDRSolver& solver);
  virtual ~ProbeWriter4();

  void startPlotting(double time) override;
  void finishPlotting() override;
  void mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp) override;
};

#endif /* POSTPROCESSING_ProbeWriter4_CLASS_HEADER_ */