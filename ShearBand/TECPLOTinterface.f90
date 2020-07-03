RECURSIVE SUBROUTINE ElementCallTECPLOTPLOTTER(wh_in,lx0,ldx,limiter)
	USE TECPLOTPLOTTERmod
	USE MainVariables, only : nVar, nDim
	implicit none
	REAL, INTENT(IN) :: wh_in(nVar,nDOFm),lx0(nDim),ldx(nDim)
	real :: wh(nVar,nDOFm)
	real :: lx0_3(3),ldx_3(3)
	integer :: limiter
	wh=wh_in
	lx0_3=0.
	ldx_3=0.
	lx0_3(1:nDim)=lx0
	ldx_3(1:nDim)=ldx
	CALL ElementTECPLOTPLOTTER(wh,lx0_3,ldx_3,limiter)
END SUBROUTINE ElementCallTECPLOTPLOTTER

RECURSIVE SUBROUTINE ElementCallTECPLOTFVPLOTTER(wh,lx0,ldx,limiter)
	USE TECPLOTPLOTTERmod
	USE MainVariables, only : nVar, nDim
	implicit none
	REAL, INTENT(IN) :: wh(nVar,(nSubLim+2*nSubLim_GL)**nDIM),lx0(nDim),ldx(nDim)
	real :: lx0_3(3),ldx_3(3)
	integer :: limiter,i,j,k,Stencil
	lx0_3(1:nDim)=lx0
	ldx_3(1:nDim)=ldx
	CALL ElementTECPLOTPLOTTER_FV(wh,lx0_3,ldx_3,limiter)
END SUBROUTINE ElementCallTECPLOTFVPLOTTER


RECURSIVE SUBROUTINE InitializeTECPLOTPLOTTER(time)
	USE TECPLOTPLOTTERmod
	implicit none
	real :: time
	CALL InitTECPLOTPLOTTER(time)
END SUBROUTINE InitializeTECPLOTPLOTTER

RECURSIVE SUBROUTINE FinishTECPLOTPLOTTER(Myrank)
	USE TECPLOTPLOTTERmod
	implicit none
	integer :: Myrank
	CALL FinalizeTECPLOTPLOTTER(Myrank)
END SUBROUTINE FinishTECPLOTPLOTTER


