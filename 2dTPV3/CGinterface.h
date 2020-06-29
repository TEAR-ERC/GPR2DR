#ifndef __CGINTERFACE_GPRDR__
#define __CGINTERFACE_GPRDR__


extern "C" {

// FORTRAN functions called by C
void loadcgfile_(const double* const MyOffset, const double* const MyDomain,const double* const cx,const double* const cy, const int * const isbinary);
//void loadcgfile_(const double* const MyOffset, const double* const MyDomain,const int * length, const char* parsetup,const double* const cx,const double* const cy, const int * const isbinary);
}/* extern "C" */
#endif /* __INITIAL_DATA_ADAPTER_CPP_FORTRAN_MHD__ */