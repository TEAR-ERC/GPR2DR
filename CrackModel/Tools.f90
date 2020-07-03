! Tools.f90

RECURSIVE SUBROUTINE InitTECPLOT(N_in,M_in,SubLim_in,Ghostlayers_in)
	USE TECPLOTPLOTTERmod
	implicit none
	INTEGER :: N_in,M_in,SubLim_in,Ghostlayers_in
	CALL SetMainParameters(N_in,M_in,SubLim_in,Ghostlayers_in)
END SUBROUTINE InitTECPLOT


!RECURSIVE SUBROUTINE InitTECPLOT(N_in,SubLim_in,Ghostlayers_in)
!	USE, INTRINSIC :: ISO_C_BINDING
!	USE TECPLOTPLOTTERmod !, only : SetMainParameters
!	implicit none
!	INTEGER :: N_in,SubLim_in,Ghostlayers_in
!	CALL SetMainParameters(N_in,SubLim_in,Ghostlayers_in)
!END SUBROUTINE InitTECPLOT

    

RECURSIVE SUBROUTINE getNumericalSolution(V,Q) 
  USE MainVariables, ONLY: nVar  
  IMPLICIT NONE     
  REAL				:: V(nVar), Q(nVar)
  INTEGER 			:: iErr
  !
  CALL PDECons2Prim(V,Q,iErr)
  !
END SUBROUTINE getNumericalSolution

RECURSIVE SUBROUTINE getExactSolution(x,timeStamp,V) 
  USE MainVariables, ONLY: nAux,nVar,nDim
  IMPLICIT NONE     
  REAL				:: V(nVar), Q(nVar), x(nDim), timeStamp
  INTEGER 			:: iErr
  !
  call InitialData(x, timeStamp, Q)
  CALL PDECons2Prim(V,Q,iErr)
  ! 
END SUBROUTINE getExactSolution


