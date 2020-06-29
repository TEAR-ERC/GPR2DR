// Tools.h  
// Fortran functions:
extern "C" {
//
void getnumericalsolution_(double* V,double* Q);
void getexactsolution_(double* x,double* timestep,double* V);
void inittecplot_(const int* N_in,const int* M_in, const int* basisSize, const int* Ghostlayers);
//void inittecplot_(const int* N_in,const int* M_in);
//
}/* extern "C" */
  