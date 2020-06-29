#ifndef HEADER_MAINVARIABLES
#define HEADER_MAINVARIABLES

! DIM Parameters
  
  MODULE MainVariables 
    IMPLICIT NONE  
    PUBLIC  

    ! This typesDef.f90 is a relict still from the old Fortran interface in Exahype.
    ! However, the following two variables are needed in the Fortran code. They
    ! should provided by the glue code generator in an appropriate way.

    ! If you modify the SRHD.exahype, please do the following mapping by hand:
    !
    ! solver ADER-DG SRHDSolver / unknowns   ->  goes to ->  nVar
    ! computational-domain / dimension       ->  goes to ->  nDim
    !

    ! Here, we obtain DIMENSIONS from the C++ preprocessor
#if defined(Dim3)
    INTEGER, PARAMETER             	:: nDim = 3                   ! The number of space dimensions
#elif defined(Dim2)
    INTEGER, PARAMETER             	:: nDim = 2                   ! The number of space dimensions
#endif
	INTEGER, PARAMETER             	:: nAux = 16
    INTEGER, PARAMETER             	:: nVar = 27                           ! The number of variables of the PDE system 
    INTEGER, PARAMETER 				:: nLin = 7
    !CHARACTER(LEN=20), PARAMETER	:: ICType='NLOPRUPTURE'
	CHARACTER(LEN=20)				:: ICType
	LOGICAL							:: USECG
    TYPE tEquation 
        ! COnsistency variables
        REAL    :: PI, LAMBDA, MU,C0,cv2,gamma2,k1,k2, phi, Cf, Cr, theta , g, ch
        REAL    :: Su, Pr, Sc, kappa, R, eta, taul, taus, cleaning_a, cleaning_eps
        REAL    :: gamma1, pi1, pi2, sigma
        INTEGER :: SubType                          ! Equation subtype
        ! Specific for rupture processes
        REAL :: EPSILON1
!	REAL :: CoulombFLpar(10)    		! add new variables for Coulomb friction
!	CoulombFLpar=0.
!	CoulombFLpar(1)=0.6
!	CoulombFLpar(2)=0.2
!	CoulombFLpar(3)=0.5
        REAL    :: lam1,lam2, lambda1,mu1,mu2
        REAL(8), DIMENSION(9) :: F0
        REAL(8), DIMENSION(3) :: F0J
        REAL(8) :: tau0, KK, Y0, Y1, aexp, xieps, gamma, p0, cv, rho0, taumin,tau2
        REAL(8) :: alpha1, beta1, alpha2, beta2, CL1, CS1, CL2, CS2
        REAL(8) :: Yeq_A, Yeq_B, Yeq_C, Yeq_s0, Yeq_smoothing
        INTEGER :: Yeq_mode, EOS, JACOBIAN_MODE
        REAL    :: T0
        REAL(8) :: mu1mu2ratio,lam1lam2ratio
        REAL    :: gvec(3),Y0lim
        INTEGER :: nMATs
        CHARACTER(LEN=50), ALLOCATABLE :: MATERIALS(:)
        LOGICAL :: FLAG
    END TYPE tEquation
	
	TYPE(tEquation) :: EQN
	
  ! 3-point Gaussian quadrature 
  REAL, PARAMETER     :: sGP3(3) = (/ 0.5-sqrt(15.)/10.,0.5,0.5+sqrt(15.)/10. /) 
  REAL, PARAMETER     :: wGP3(3) = (/ 5./18., 8./18., 5./18. /) 
  INTEGER, PARAMETER  :: nGP3 = 3 

  END MODULE MainVariables  
#endif
