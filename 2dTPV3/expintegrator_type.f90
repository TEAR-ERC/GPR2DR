#ifndef HEADER_EXPINTEGRATOR_TYPE
#define HEADER_EXPINTEGRATOR_TYPE
#if defined(ODESOLVER) || defined(EQNTYPEC99) || defined(EQNTYPED99)
    module expintegrator_type
#if defined(EQNTYPEB99)
        INTEGER, PARAMETER :: nvarode = 10
        INTEGER, PARAMETER :: nccc = 14
        INTEGER, PARAMETER :: ncc  = 1
#endif 

#if defined(EQNTYPEC99) || defined(EQNTYPED99)
        INTEGER, PARAMETER :: nvarode = 10
        INTEGER, PARAMETER :: nccc = 14
        INTEGER, PARAMETER :: ncc  = 1
#endif 

#ifdef EQNTYPE7
        !INTEGER, PARAMETER :: nvarode = ...
        !INTEGER, PARAMETER :: nccc = ...
        !INTEGER, PARAMETER :: ncc  =.. 
#endif

    TYPE T_odeexp_settings
        REAL(8) :: delta_max, eps_fpi, eps_ccc(nccc), increment_ccfl, decrement_accuracy, decrement_positivity
        INTEGER :: maxiter_accept, maxiter_reject
    END TYPE T_odeexp_settings
    TYPE(T_odeexp_settings) :: ODESETTINGS    
    
    TYPE T_ode_parameters 
        REAL(8), DIMENSION(9) :: F0
        REAL(8), DIMENSION(3) :: F0J
        REAL(8) :: tau0, KK, Y0, Y1, aexp, xieps, rho0, taumin
        REAL(8) :: alpha1, beta1, alpha2, beta2, lam1, mu1, lam2, mu2, K1, K2
        REAL(8) :: Yeq_A, Yeq_B, Yeq_C, Yeq_s0, Yeq_smoothing
        INTEGER :: Yeq_mode, JACOBIAN_MODE
        REAL(8) :: ALPHA
		REAL(8) :: mu1mu2ratio,lam1lam2ratio
        !
        REAL(8) :: gamma,p0,cv,tau2,S
        integer :: EOS
        real(8) :: T0
    END TYPE T_ode_parameters
        
    end module expintegrator_type
#endif    
#endif
