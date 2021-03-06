// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ProbeWriter4.h"
#include "PDE.h"
#include "C2P-GPRDR.h"
#include <algorithm>
#include <stdlib.h>


GPRDR::ProbeWriter4::ProbeWriter4(GPRDR::GPRDRSolver& solver) {
  // @TODO Please insert your code here.
}

GPRDR::ProbeWriter4::~ProbeWriter4() {
}

void GPRDR::ProbeWriter4::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void GPRDR::ProbeWriter4::finishPlotting() {
  // @TODO Please insert your code here.
}

void GPRDR::ProbeWriter4::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  const int writtenUnknowns = 28;
  double V[28];

  pdecons2prim_(V,Q);

  for (int i=0; i<writtenUnknowns; i++){ 
    outputQuantities[i] = std::max(-1e36,std::min(1e36,V[i]));
  }
  double aux[16];
  double x_3[3];
  x_3[2]=0;
  std::copy_n(&x[0],DIMENSIONS,&x_3[0]);

  pdeauxvar_(aux,Q,x_3,&timeStamp);

  for (int i=writtenUnknowns; i<writtenUnknowns+16; i++){
    outputQuantities[i] = std::max(-1e36,std::min(1e36,aux[i-writtenUnknowns]));

}
}
