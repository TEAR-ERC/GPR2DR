#ifndef __EXAHYPE_USER_ODE__
#define __EXAHYPE_USER_ODE__

// Fortran functions:
extern "C" {
void updatesolutionode_(double* Qnew, const double* const Q0,const double* const dt);
}/* extern "C" */

#endif /* __EXAHYPE_USER_PDE__ */
