#include "GPRDRSolver_FV.h"
#include "kernels/finitevolumes/musclhancock/c/musclhancock.h"
#include "GPRDRSolver_FV_Variables.h"
#include "InitialData.h"
#include "PDE.h"
#include "ODE.h"


tarch::logging::Log GPRDR::GPRDRSolver_FV::_log( "GPRDR::GPRDRSolver_FV" );

void GPRDR::GPRDRSolver_FV::init(const std::vector<std::string>& cmdlineargs,const exahype::parser::ParserView& constants) {
}

void GPRDR::GPRDRSolver_FV::adjustSolution(const double* const x,const double t,const double dt, double* const Q) {
   const int nVar = GPRDR::GPRDRSolver_FV::NumberOfVariables;

    if (tarch::la::equals(t,0.0)) {
    int md = exahype::solvers::Solver::getMaximumAdaptiveMeshDepth();
    double cms = exahype::solvers::Solver::getCoarsestMeshSize();
    const int order = GPRDR::GPRDRSolver_FV::PatchSize;

    std::fill_n(Q,nVar,0.0);
    
    //    initialdata_(x, &ti, Qgp,&md,&cms,&order);
    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
    
    initialdata_(x_3, &t, Q);
  }

//    dynamicrupture_(x,&t,Q);//Duo April 10

    for(int i = 0; i< nVar; i++){
      assert(std::isfinite(Q[i]));
    }
}

void GPRDR::GPRDRSolver_FV::eigenvalues(const double* const Q, const int dIndex, double* const lambda) {
  // Tip: You find documentation for this method in header file "GPRDR::GPRDRSolver_FV.h".
  // Tip: See header file "GPRDR::AbstractGPRDRSolver_FV.h" for toolkit generated compile-time 
  //      constants such as PatchSize, NumberOfVariables, and NumberOfParameters.
  
  // @todo Please implement/augment if required
  double nv[3] = {0.};
  nv[dIndex] = 1;
  pdeeigenvalues_(lambda, Q, nv);
}

void GPRDR::GPRDRSolver_FV::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int direction,
    const double* const stateInside,
    double* const stateOutside) {
  const int nVar = GPRDR::GPRDRSolver_FV::NumberOfVariables;
 
	double Qgp[nVar];

	double ti = t + 0.5 * dt;
	// Compute the outer state according to the initial condition
    double x_3[3];
    x_3[2]=0;
    std::copy_n(&x[0],DIMENSIONS,&x_3[0]);
    
    initialdata_(x_3, &t, Qgp);
	// Assign the proper outer state
	for(int m=0; m < nVar; m++) {
        stateOutside[m] = Qgp[m];
	}
	 //std::copy_n(stateInside,nVar,stateOutside);

}

//***********************************************************
//*********************** PDE *******************************
//***********************************************************

//to add new PDEs specify them in the specification file, delete this file and its header and rerun the toolkit


void GPRDR::GPRDRSolver_FV::flux(const double* const Q,double** const F) {
  const int nVar = GPRDR::GPRDRSolver_FV::NumberOfVariables;
  if(DIMENSIONS == 2){
    double F_3[nVar];
    pdeflux_(F[0], F[1],F_3, Q);
  }else{
    pdeflux_(F[0], F[1],F[2], Q);
  }
}




//You can either implement this method or modify fusedSource
void GPRDR::GPRDRSolver_FV::algebraicSource(const tarch::la::Vector<DIMENSIONS, double>& x, double t, const double *const Q, double *S) {
  pdesource_(S, Q);
}

void  GPRDR::GPRDRSolver_FV::nonConservativeProduct(const double* const Q,const double* const gradQ,double* const BgradQ) {
  pdencp_(BgradQ, Q, gradQ);
}


void GPRDR::GPRDRSolver_FV::solutionUpdate(double* luh,const tarch::la::Vector<DIMENSIONS,double>& cellCenter,const tarch::la::Vector<DIMENSIONS,double>& cellSize,const double t, const double dt,double& maxAdmissibleDt) {
  GPRDR::AbstractGPRDRSolver_FV::solutionUpdate(luh, cellCenter, cellSize, t, dt, maxAdmissibleDt);
  constexpr int patchSize          = GPRDR::GPRDRSolver_FV::PatchSize;
  constexpr int ghostLayerWidth    = GPRDR::GPRDRSolver_FV::GhostLayerWidth;
  constexpr int patchBegin         = ghostLayerWidth; // patchBegin cell is inside domain
  constexpr int patchEnd           = patchBegin+patchSize; // patchEnd cell is outside domain
  constexpr int wholePatchSize      = patchSize+2*ghostLayerWidth;
  double x[3];
  double slp_p[3],slp_m[3];
  double slp;
  slp = 0.0;
  x[2] = 0.0;

#if DIMENSIONS==3
  kernels::idx4 idx(wholePatchSize,wholePatchSize,wholePatchSize,NumberOfVariables);
  for (int k = patchBegin; k < patchEnd; k++) {
#else
    kernels::idx3 idx(wholePatchSize,wholePatchSize,NumberOfVariables);
#endif    
    for (int j = patchBegin; j < patchEnd; j++) {
      for (int i = patchBegin; i < patchEnd; i++) {
#if DIMENSIONS==3
	double* luh_cell = luh + idx(k,j,i,0);
#else
	double* luh_cell = luh + idx(j,i,0);
    double* luh_m = luh+idx(patchEnd-j,i,0);

    slp = std::abs( (luh_cell[24]-luh_m[24]) );
//    std::cout<<"slip="<<slp<<"\n";
#endif
	double luh_cell_new[NumberOfVariables];
    //std::cout<<"patchSize"<<patchSize<<"\n";

      x[0] = cellSize[0]/patchSize*(0.5+i-patchBegin)+cellCenter[0]-0.5*cellSize[0];
      x[1] = cellSize[1]/patchSize*(0.5+j-patchBegin)+cellCenter[1]-0.5*cellSize[1];
      dynamicrupture_(x,&t,&luh_cell[0],&slp);   //Duo April 28


	updatesolutionode_(&luh_cell_new[0],luh_cell,&dt);
	std::copy_n(luh_cell_new,NumberOfVariables,luh_cell);

      }
    }
#if DIMENSIONS==3
  }
#endif  
  
}


#include "kernels/finitevolumes/riemannsolvers/c/riemannsolvers.h"

double GPRDR::GPRDRSolver_FV::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, const double* gradQL, const double* gradQR, const double* cellSize, int direction) {


    //return kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, false, GRMHDbSolver_FV>(*static_cast<GRMHDbSolver_FV*>(this), fL,fR,qL,qR,gradQL, gradQR, cellSize, direction);
	constexpr int numberOfVariables = AbstractGPRDRSolver_FV::NumberOfVariables;

	//printf("SONO QUI IN riemannSolver");
	/* HLLEM */
	
	//const int numberOfVariables = GRMHDb::AbstractGRMHDbSolver_FV::NumberOfVariables;
	const int numberOfParameters = GPRDR::AbstractGPRDRSolver_FV::NumberOfParameters;
	const int numberOfData = numberOfVariables + numberOfParameters;
	const int order = 0;  // for finite volume we use one single d.o.f., i.e. the cell average.
	const int basisSize = order + 1;
	// Compute the average variables and parameters from the left and the right
	double QavL[numberOfData] = { 0.0 }; // ~(numberOfVariables+numberOfParameters)
	double QavR[numberOfData] = { 0.0 }; // ~(numberOfVariables+numberOfParameters)

										 // std::cout << "opened ---------------------"<< std::endl;
         
        // printf("\n******* RIEMANN SOLVER FV*****************");

	//kernels::idx2 idx_QLR(basisSize, numberOfData);
	//for (int j = 0; j < basisSize; j++) {
	//	const double weight = kernels::legendre::weights[order][j];
    //
	//	for (int k = 0; k < numberOfData; k++) {
	//		QavL[k] += weight * qL[idx_QLR(j, k)];
	//		QavR[k] += weight * qR[idx_QLR(j, k)];
	//	}
	//}
	
        // printf("\n***DONE*****");

	//	double lambda = 2.0;
	
	//hllemfluxfv_(fL, fR, qL, qR, QavL, QavR, &direction);

	hllemfluxfv_(fL, fR, qL, qR, &direction);
	//double lambda = kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, false, GPRDRSolver_FV>(*static_cast<GPRDRSolver_FV*>(this), fL, fR, qL, qR, gradQL, gradQR, cellSize, direction);
	/* OSHER */
	//double lambda = kernels::finitevolumes::riemannsolvers::c::generalisedOsherSolomon<false, true, false, 3, EulerSolver_FV>(*static_cast<EulerSolver_FV*>(this), fL, fR, qL, qR, direction);
	/* RUSANOV */
	
	// avoid spurious numerical diffusion (ony for Cowling approximation)
	//for (int m = 17; m < numberOfVariables; m++) {
	//	fL[m] = 0.0;
	//	fR[m] = 0.0;
	//}
	double lambda = 6000.0;
	return lambda; 
}

//double GPRDR::GPRDRSolver_FV::riemannSolver(double* fL, double *fR, const double* qL, const double* qR, const double* gradQL, const double* gradQR, const double* cellSize, int direction) {
//  const int numberOfVariables  = GPRDR::AbstractGPRDRSolver_FV::NumberOfVariables;
//  const int numberOfParameters = GPRDR::AbstractGPRDRSolver_FV::NumberOfParameters;
//  const int numberOfData       = numberOfVariables+numberOfParameters;
//  const int order              = 0;
//  const int basisSize          = order+1;
//  // Compute the average variables and parameters from the left and the right
//  double QavL[numberOfData] = {0.0}; // ~(numberOfVariables+numberOfParameters)
//  double QavR[numberOfData] = {0.0}; // ~(numberOfVariables+numberOfParameters)
//  
//  // std::cout << "opened ---------------------"<< std::endl;
//  
//    kernels::idx2 idx_QLR(basisSize, numberOfData);
//    for (int j = 0; j < basisSize; j++) {
//      const double weight = kernels::legendre::weights[order][j];
//
//      for (int k = 0; k < numberOfData; k++) {
//        QavL[k] += weight * qL[idx_QLR(j, k)];
//        QavR[k] += weight * qR[idx_QLR(j, k)];
//      }
//    }
//double lambda = kernels::finitevolumes::riemannsolvers::c::rusanov<true, true, false, GPRDRSolver_FV>(*static_cast<GPRDRSolver_FV*>(this), fL, fR, qL, qR, gradQL, gradQR, cellSize, direction);
//// Call the Fortran routine
//hllemriemannsolver_(&basisSize, &direction, fL,fR,qL, qR,QavL, QavR);
////testriemannsolver_(&basisSize, &normalNonZero, fL,fR,qL, qR,QavL, QavR);
//return lambda;
//}
