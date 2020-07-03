// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ProbeWriter1.h"
#include "PDE.h"
#include "C2P-GPRDR.h"
#include <algorithm>

GPRDR::ProbeWriter1::ProbeWriter1(GPRDR::GPRDRSolver& solver) {
  // @TODO Please insert your code here.
}

GPRDR::ProbeWriter1::~ProbeWriter1() {
}

void GPRDR::ProbeWriter1::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void GPRDR::ProbeWriter1::finishPlotting() {
  // @TODO Please insert your code here.
}

void GPRDR::ProbeWriter1::mapQuantities(
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
    //outputQuantities[i] = std::max(-1e36,std::min(1e36,V[i]));
     outputQuantities[i] = V[i];
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
