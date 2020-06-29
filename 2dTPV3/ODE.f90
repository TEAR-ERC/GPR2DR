#ifndef HEADER_ODE
#define HEADER_ODE

!#include "MainVariables.f90"
!#include "expintegrator_type.f90"
!#include "expintegrator.f90" 
!#include "expintegrator_ode.f90"


RECURSIVE SUBROUTINE UpdateSolutionODE(Qnew,Qold,loc_dt)
  USE MainVariables, ONLY : nVar, nDim, EQN
#ifdef ODESOLVER
    use expintegrator_type, only : ODESETTINGS, nVarODE                                     , t_ode_parameters
    use expintegrator, only : expintegrator_adaptive
    USE expintegrator_ode, only : export_expintegrator_state,extract_expintegrator_state     , convert_parameters
#endif
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: Qold(nVar), Qnew(nVar),loc_dt
  INTENT(IN)  :: Qold,loc_dt
  INTENT(OUT) :: Qnew
  ! Local varialbes
#ifdef ODESOLVER
    real :: AC(nVarODE),AC0(nVarODE)
    INTEGER :: LOCiter, substep_count
    
    real ::V0(nVar),Q0(nVar)
    type(t_ode_parameters) :: ODE
#endif
	REAL :: x00(3),time
	INTEGER :: iErr,i
	
	x00=0.
	time=0.
	Q0=Qold
	
	if(maxval(abs(Q0))<1.e-3) then
		Qnew=Qold
		return
	end if
    !Q0 = subuh1(:,ii,jj,kk)    
    call extract_expintegrator_state(AC0, Q0)
    
    CALL PDECons2Prim(V0,Q0,x00,time,iErr)
    CALL convert_parameters(V0, ODE, EQN)
	!if(AC0(1)>1.e-3) then
	!	print *, AC0,loc_dt
	!	
	!end if
	!AC0(1)=0.271809817913919
	!
	!AC0(2)=0.999711778985515
	!AC0(3)= -1.424312079886052E-003
	!AC0(4)= -1.403901704045918E-020
	!AC0(5)= 1.336400383137943E-003
	!AC0(6)= 0.999633322663082
	!AC0(7)= 2.077140788241882E-018
	!AC0(8)= 1.133080272349694E-020
	!AC0(9)= -2.076754583610124E-018
	!AC0(10)= 0.999679787557563
	!print *, ' -INPUT- ',Qold
	!LOCiter=1
    CALL expintegrator_adaptive(AC0, EQN, loc_dt, ODESETTINGS,AC, substep_count, LOCiter, ODE)
	!AC=AC0
	!if(AC0(1)>1.e-3) then
	!print *, ' -OUTPUT- ' ,AC,' --------- ',loc_dt,' --------- ', LOCiter, substep_count
		
	!end if
	if(LOCiter<0) then
		print *, ' -OUTPUT- ' ,LOCiter,loc_dt, nVar
		print *, AC0
		print *,'ODESETTINGS%delta_max           '  , ODESETTINGS%delta_max
		print *,'ODESETTINGS%eps_fpi             '  , ODESETTINGS%eps_fpi
		print *,'ODESETTINGS%eps_ccc             '  , ODESETTINGS%eps_ccc
		print *,'ODESETTINGS%increment_ccfl      '  , ODESETTINGS%increment_ccfl
		print *,'ODESETTINGS%decrement_accuracy  '  , ODESETTINGS%decrement_accuracy
		print *,'ODESETTINGS%decrement_positivity'  , ODESETTINGS%decrement_positivity
		print *,'ODESETTINGS%maxiter_accept      '  , ODESETTINGS%maxiter_accept
		print *,'ODESETTINGS%maxiter_reject      '  , ODESETTINGS%maxiter_reject
		print *,'ODE%KK                          '  , ODE%KK
		print *,'ODE%tau0                        '  , ODE%tau0 
		print *,'ODE%Y0                          '  , ODE%Y0
		print *,'ODE%Y1                          '  , ODE%Y1
		print *,'ODE%tau2                        '  , ODE%tau2
		print *,'ODE%EOS                         '  , ODE%EOS
		print *,'ODE%T0                          '  , ODE%T0
		print *,'ODE%lam1                        '  , ODE%lam1
		print *,'ODE%mu1                         '  , ODE%mu1
		print *,'ODE%lam2                        '  , ODE%lam2
		print *,'ODE%mu2                         '  , ODE%mu2
		print *,'ODE%K1                          '  , ODE%K1
		print *,'ODE%K2                          '  , ODE%K2
		
		stop
	end if
	!AC=AC0
    call export_expintegrator_state(Q0, AC,(/ substep_count, LOCiter /))
    !subuh1(:,ii,jj,kk) = Q0
	!stop
	
    Qnew=Q0
    ! Duo dQ24/dt = -Q23
    ! Qnew(24)=Qold(24)- loc_dt*(Qnew(23)+Qold(23))/2

	!Qnew(24)=24
	!Qnew(23)=23
	!do i=22,nVar
	!	Qnew(i)=i*10.0
	!end do
  END SUBROUTINE UpdateSolutionODE

#endif
