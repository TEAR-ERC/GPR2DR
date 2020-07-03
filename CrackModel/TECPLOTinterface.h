#ifndef __EXAHYPE_USER_TECINT__
#define __EXAHYPE_USER_TECINT__


// Fortran functions:
extern "C" {
	//void inittecplot_(const int* N_in, const int* M_in, const int* basisSize, const int* Ghostlayers);
    void elementcalltecplotplotter_(const double *wh, const double* lx0, const double* ldx, const int* limiter);
	void elementcalltecplotaderdgplotter_(const double *wh, const double* lx0, const double* ldx, const int* limiter);
	void elementcalltecplotfvplotter_(const double *wh, const double* lx0, const double* ldx, const int* limiter);
void finishtecplotplotter_(const int* Myrank);
void initializetecplotplotter_(const double* time);
}/* extern "C" */

#endif /* __EXAHYPE_USER_PDE__ */
