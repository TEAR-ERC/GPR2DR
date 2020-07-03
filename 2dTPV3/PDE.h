#ifndef __EXAHYPE_USER_PDE__
#define __EXAHYPE_USER_PDE__

// Fortran functions:
extern "C" {
void minimumtreedepth_(int* depth);
void hastoadjustsolution_(double* t, bool* refine);
void adjustedsolutionvalues_(const double* const x,const double* w,const double* t,const double* dt,double* Q);
void pdeflux_(double* Fx, double* Fy, double* Fz, const double* const Q);
void pdesource_(double* S, const double* const Q);
void pdencp_(double* BgradQ, const double* const Q, const double* const gradQ);
void pdeeigenvalues_(double* lambda, const double* const Q, double* nv);
void pdevarname_(char* MyNameOUT, int* ind);
void registerinitialdata_(const char* const id_name, int* id_name_len);
//void getnumericalsolution_(double* V,double* Q);
//void getexactsolution_(double* V,double* pos,double* timeStamp);

//void inittecplot_(const int* N_in,const int* M_in);
void pdeauxvar_(double* aux, const double* const Q,double* x, const double* const time);
void pdecritialstress_(double* CS, const double* const Q);
void dynamicrupture_(const double* x, const double* const t, double* Q,const double* slp);


//void hllemfluxfv_(double* FL, double* FR, const double* const  QL, const double* const  QR, const double* const  QavL, const double* const  QavR, const int* normalNonZeroIndex);
void hllemfluxfv_(double* FL, double* FR, const double* const  QL, const double* const  QR, const int* normalNonZeroIndex);
void hllemriemannsolver_(const int* basisSize, const int* normalNonZeroIndex, double* FL, double* FR, const double* const  QL, const double* const  QR, const double* const  QavL, const double* const  QavR);

}/* extern "C" */

#endif /* __EXAHYPE_USER_PDE__ */
