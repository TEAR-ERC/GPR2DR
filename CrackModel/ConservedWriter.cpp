// This file was generated by the ExaHyPE toolkit.
// It will not be overwritten.
//
//
// ========================
//   www.exahype.eu
// ========================
#include "ConservedWriter.h"
#include "PDE.h"
#include "C2P-GPRDR.h"
#include <algorithm>
#include <stdlib.h> 

GPRDR::ConservedWriter::ConservedWriter(GPRDR::GPRDRSolver& solver) {
  // @TODO Please insert your code here.
}

GPRDR::ConservedWriter::~ConservedWriter() {
}

void GPRDR::ConservedWriter::startPlotting( double time) {
  // @TODO Please insert your code here.
}


void GPRDR::ConservedWriter::finishPlotting() {
  // @TODO Please insert your code here.
}

void GPRDR::ConservedWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* const Q,
    double* const outputQuantities,
    double timeStamp
) {
  const double epsi=1.0e-36;
  double V[28];

  pdecons2prim_(V,Q);

  const int writtenUnknowns = 28;
//  for (int i=0; i<writtenUnknowns; i++){ 
//    if (std::abs(V[i])<epsi){
//      outputQuantities[i]=0.0;
//    }else {
//      outputQuantities[i]=std::max(-1e36,std::min(1e36,V[i]));
//    }
//  }

  double aux[16];
  outputQuantities[0]=V[1]; // u
  outputQuantities[1]=V[24]; // U
  outputQuantities[2]=V[27]; //gradU
  double x_3[3];
  x_3[2]=0;
  std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
  pdeauxvar_(aux,Q,x_3,&timeStamp);

   outputQuantities[3]=V[20]; //xi
   outputQuantities[4]=aux[3]; //sxy
//  for (int i=writtenUnknowns; i<writtenUnknowns+16; i++){ 
//    if(std::abs(aux[i-writtenUnknowns])<epsi){
//       outputQuantities[i]=0.0;
//    }else {
//    outputQuantities[i] = std::max(-1e36,std::min(1e36,aux[i-writtenUnknowns]));
//    }
//}
}
