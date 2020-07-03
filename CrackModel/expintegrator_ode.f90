#ifndef HEADER_EXPINTEGRATOR_ODE
#define HEADER_EXPINTEGRATOR_ODE
#if defined(ODESOLVER) || defined(EQNTYPEC99) || defined(EQNTYPEB99) || defined(EQNTYPED99)
!#include "MainVariables.f90"
!#include "expintegrator_type.f90"
!#include "expintegrator_linalg.f90"
!#include "expintegrator_linalg_fast33.f90"

MODULE expintegrator_ode
    use MainVariables, only: tEquation
    use expintegrator_type, only : T_ode_parameters,nvarode,ncc,nccc
    USE expintegrator_linalg, only: compute_matrix_exponential, solve, determinant_lu, build_jacobian_matrix_mask, &
        build_reduced_jacobian_matrix, compute_first_order_linear_ode_solution, select_solution_from_mask, halven_rank42_rows, &
        flatten_rows, build_by_rows, cross_product, diadic_product, trace
    USE expintegrator_linalg_fast33, only: fast33_compute_matrix_exponential, fast33_build_jacobian_matrix_mask, &
        fast33_build_reduced_jacobian_matrix, fast33_compute_first_order_linear_ode_solution, fast33_select_solution_from_mask
    USE, INTRINSIC :: ieee_arithmetic

    IMPLICIT NONE



    

    PRIVATE
    !
    !
    PUBLIC :: convert_parameters
    PUBLIC :: system_linearization
    PUBLIC :: linearized_system_solution
    PUBLIC :: timescale
    PUBLIC :: test_positivity
    PUBLIC :: compute_source_jacobian_approximate
    PUBLIC :: compute_source_jacobian_A
    PUBLIC :: compute_source_jacobian_A_approximate
    PUBLIC :: source
    !PUBLIC :: compute_l_m_dedc
    !PUBLIC :: compute_total_stress
    PUBLIC :: extract_expintegrator_state
    PUBLIC :: export_expintegrator_state
    !PUBLIC :: compute_hydrodynamic_pressure
CONTAINS

    PURE SUBROUTINE extract_expintegrator_state(Vode, Vpde)
        ! This subroutine extracts the variables on which the exponential integrator should act
        REAL(8), DIMENSION(:),  INTENT(OUT)   :: Vode
        REAL(8), DIMENSION(:),  INTENT(IN)    :: Vpde
        
#if defined(EQNTYPEB99)
        Vode(1)=max(0.0,min(1.0,Vpde(18))) 
        Vode(2:10)=Vpde(5:13)
#endif

#if defined(EQNTYPEC99) || defined(EQNTYPED99)
        Vode(1)=max(0.0,min(1.0,Vpde(21))) 
        Vode(2:10)=Vpde(9:17)
        !Vode(11:13)=Vpde(6:8)
#endif

#ifdef EQNTYPE7
        !Vode(1)=...
#endif
    END SUBROUTINE extract_expintegrator_state

    PURE SUBROUTINE export_expintegrator_state(Vpde, Vode, Param)
        ! This subroutine exports the variables on which the exponential integrator acted
        REAL(8), DIMENSION(:),  INTENT(OUT) :: Vpde
        REAL(8), DIMENSION(:),  INTENT(IN)  :: Vode
        INTEGER, DIMENSION(:),  INTENT(IN)  :: Param
        
#if defined(EQNTYPEB99)
        Vpde(18)    = Vode(1)
        Vpde(5:13)  = Vode(2:10)
#endif

#if defined(EQNTYPEC99)
        Vpde(21)    = Vode(1)
        Vpde(9:17)  = Vode(2:10)
        !Vpde(6:8)   = Vode(11:13)
        Vpde(36)    = Param(2)
#endif
#if defined(EQNTYPED99)
        Vpde(21)    = Vode(1)
        Vpde(9:17)  = Vode(2:10)
        !Vpde(6:8)   = Vode(11:13)
        !Vpde(24)    = real(Param(2),8)
#endif
#ifdef EQNTYPE7

#endif
    END SUBROUTINE export_expintegrator_state

    PURE SUBROUTINE convert_parameters(Q0, ODE, EQN)
#ifdef EQNTYPED99
    USE GPRmaterials, only : AssignMaterialPropertiesMix
#endif
        ! An interface between the T_ode_parameters and T_equation_parameters data structures.
        REAL(8), DIMENSION(:), INTENT(IN)  :: Q0
        TYPE(tEquation), INTENT(IN)              :: EQN
        TYPE(T_ode_parameters),      INTENT(OUT) :: ODE
        REAL(8), DIMENSION(14)                  :: RUPTURE_PARAMS
#if defined(EQNTYPEB99)
        ODE%rho0   = EQN%rho0
        ODE%tau0   = EQN%tau0    
        ODE%KK     = EQN%KK      
        ODE%Y0     = EQN%Y0      
        ODE%Y1     = EQN%Y1      
        ODE%aexp   = EQN%aexp    
        ODE%alpha1 = EQN%alpha1  
        ODE%beta1  = EQN%beta1   
        ODE%alpha2 = EQN%alpha2  
        ODE%beta2  = EQN%beta2   
        ODE%xieps  = EQN%xieps
        ODE%taumin = EQN%taumin
        ODE%F0     = EQN%F0
        ODE%mu1    = EQN%rho0*EQN%CS1**2
        ODE%lam1   = EQN%rho0*(EQN%CL1**2 - 2.0d0*EQN%CS1**2)
        ODE%mu1    = EQN%rho0*EQN%CS1**2
        !ODE%lam2   = EQN%rho0*(EQN%CL2**2 - 2.0d0*EQN%CS2**2)
        !ODE%mu2    = EQN%rho0* EQN%CS2**2
		ODE%mu2    =EQN%mu2  
		ODE%lam2   =EQN%lam2
        ODE%K1     = ODE%lam1 + 2.0d0*ODE%mu1/3.0d0
        ODE%K2     = ODE%lam2 + 2.0d0*ODE%mu2/3.0d0
        ODE%Yeq_mode = EQN%Yeq_mode
        ODE%Yeq_A    = EQN%Yeq_A
        ODE%Yeq_B    = EQN%Yeq_B
        ODE%Yeq_C    = EQN%Yeq_C
        ODE%Yeq_s0   = EQN%Yeq_s0
        ODE%jacobian_mode = EQN%jacobian_mode
        ODE%Yeq_smoothing =EQN%Yeq_smoothing
#endif  
#ifdef EQNTYPEC99
        ODE%F0J     = EQN%F0J
        ODE%mu1    = Q0(20)
        ODE%lam1   = Q0(19)

		ODE%mu2    = Q0(20)*Q0(22)                      ! Warning, this was /
		ODE%lam2   = Q0(19)+Q0(23)*(ODE%mu1-ODE%mu2) 
        
        ODE%K1     = ODE%lam1 + 2.0d0*ODE%mu1/3.0d0
        ODE%K2     = ODE%lam2 + 2.0d0*ODE%mu2/3.0d0

        ODE%S      = Q0(5)
        ODE%rho0   = Q0(1)
        ODE%tau2   = EQN%tau2 
        ODE%gamma  = EQN%gamma     
        ODE%cv     = EQN%cv      
        ODE%p0     = EQN%p0 
        ODE%EOS    = EQN%EOS      
        ODE%T0     = EQN%T0
        
        
        ODE%KK     = Q0(24)      
        ODE%Y0     = Q0(25) 
        ODE%Y1     = Q0(26) 
        ODE%aexp   = Q0(27) 
        ODE%alpha1 = Q0(28) 
        ODE%beta1  = Q0(29) 
        ODE%alpha2 = Q0(30) 
        ODE%beta2  = Q0(31) 
        
        ODE%Yeq_mode = EQN%Yeq_mode
        ODE%Yeq_A  = Q0(32)
        ODE%Yeq_B  = Q0(33)
        ODE%Yeq_C  = Q0(34)
        ODE%Yeq_s0 = Q0(35)
        
        ODE%tau0   = EQN%tau0
        ODE%jacobian_mode = EQN%jacobian_mode
        ODE%Yeq_smoothing = ODE%Y0*1.e-3
        
        ODE%taumin = EQN%taumin
        ODE%xieps  = EQN%xieps
#endif
#ifdef EQNTYPED99
        call AssignMaterialPropertiesMix(RUPTURE_PARAMS,EQN%MATERIALS,(/Q0(22), 1.-Q0(22)/),EQN%nMATs,.TRUE.)
        ODE%F0J    = 0.! EQN%F0J
		ODE%F0    = 0.!EQN%F0J
        ODE%mu1    = Q0(20)
        ODE%lam1   = Q0(19)

        ODE%lam1lam2ratio=RUPTURE_PARAMS(2)
        ODE%mu1mu2ratio=RUPTURE_PARAMS(1)
        
		ODE%mu2    = Q0(20)*RUPTURE_PARAMS(1)                     ! Warning, this was /
		ODE%lam2   = Q0(19)+RUPTURE_PARAMS(2)*(ODE%mu1-ODE%mu2) 
        
        ODE%K1     = ODE%lam1 + 2.0d0*ODE%mu1/3.0d0
        ODE%K2     = ODE%lam2 + 2.0d0*ODE%mu2/3.0d0

        ODE%S      = Q0(5)
        ODE%rho0   = Q0(1)
        ODE%tau2   = EQN%tau2 
        ODE%gamma  = EQN%gamma     
        ODE%cv     = EQN%cv      
        ODE%p0     = EQN%p0 
        ODE%EOS    = EQN%EOS      
        ODE%T0     = EQN%T0
        
        ODE%alpha     = Q0(18)
        
        
        ODE%KK     = RUPTURE_PARAMS(3)*abs(Q0(18))**2    
        ODE%Y0     = RUPTURE_PARAMS(4)
        ODE%Y1     = RUPTURE_PARAMS(5)
        ODE%aexp   = RUPTURE_PARAMS(6)
        ODE%alpha1 = RUPTURE_PARAMS(7)
        ODE%beta1  = RUPTURE_PARAMS(8)
        ODE%alpha2 = RUPTURE_PARAMS(9)
        ODE%beta2  = RUPTURE_PARAMS(10)*min(1.0,max(0.0,1.0-Q0(23))) ! Duo April 10
        
        ODE%Yeq_mode = EQN%Yeq_mode
        ODE%Yeq_A  =  RUPTURE_PARAMS(11)
        ODE%Yeq_B  =  RUPTURE_PARAMS(12)
        ODE%Yeq_C  =  RUPTURE_PARAMS(13)
        ODE%Yeq_s0 =  RUPTURE_PARAMS(14)
        
        ODE%tau0   = EQN%tau0!/(abs(Q0(18))+1.e-20)
        ODE%jacobian_mode = EQN%jacobian_mode
        ODE%Yeq_smoothing = ODE%Y0*1.e-3
        
        ODE%taumin = EQN%taumin
        ODE%xieps  = EQN%xieps
#endif
#ifdef EQNTYPE7

#endif
        continue
    END SUBROUTINE convert_parameters

    PURE SUBROUTINE system_linearization(Ccc, Cc, V, ODE)
        ! This subroutine computes all or some of the linearization coefficients: 
        ! with numerical differentiation of the jacobian matrix, one could store 
        ! all the coefficients of the jacobian matrix in Ccc, but currently we 
        ! prefer to just compute some indicator variables for adaptive timestepping 
        ! and call this subroutine when evaluating the system source. 
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V
        REAL(8), DIMENSION(nccc),    INTENT(OUT) :: Ccc
        REAL(8), DIMENSION(ncc),     INTENT(OUT) :: Cc
        TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
#if defined(EQNTYPEB99) || defined(EQNTYPEC99) || defined(EQNTYPED99)
        CALL system_linearization_gprdr(Ccc, Cc, V, ODE)
#endif
#ifdef EQNTYPE7

#endif
    END SUBROUTINE system_linearization

    PURE FUNCTION source(t, Q, ODE)
        ! Source term of the full system
        REAL(8),                     INTENT(IN) :: t
        REAL(8), DIMENSION(nvarode), INTENT(IN) :: Q
        TYPE(T_ode_parameters),      INTENT(IN) :: ODE
        REAL(8), DIMENSION(nvarode)             :: source
        REAL(8), DIMENSION(nccc)                :: Ccc
        REAL(8), DIMENSION(ncc)                 :: Cc
        ! ode-dependant variables
        REAL(8), DIMENSION(3,3)                 :: A
        REAL(8)                                 :: alpha1, alpha0

        CALL system_linearization(Ccc, Cc, Q, ODE)
#if defined(EQNTYPEB99) || defined(EQNTYPEC99)  || defined(EQNTYPED99)
        A = build_by_rows(Q(2:10))
        CALL evaluate_source_A(source(2:10), A, Ccc(13)**2, ODE)
        alpha1 = sign(1.0d0, Ccc(11))*Ccc(11)**2
        alpha0 = sign(1.0d0, Ccc(12))*Ccc(12)**2*alpha1
        source(1) = alpha1*Q(1) + alpha0
#endif
#ifdef EQNTYPEC99
        !source(11:13)=Ccc(15)*Q(11:13)+Ccc(16:18)
        !source(11:13)=0.
#endif

#ifdef EQNTYPE7

#endif
    END FUNCTION source

    PURE SUBROUTINE test_positivity(V, ODE, posstatus)
        ! This subroutine detects positivity violations in the state vector V.
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V
        TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
        INTEGER,                     INTENT(OUT) :: posstatus
        REAL(8), PARAMETER                       :: eps = 0.0d0
        
#if defined(EQNTYPEB99) || defined(EQNTYPEC99) || defined(EQNTYPED99)
        IF (V(1) .ge. -eps .and. V(1) .le. 1.0d0 + eps) THEN
            posstatus = 1
        ELSE
            posstatus = -1
        END IF
#endif
#ifdef EQNTYPE7

#endif
    END SUBROUTINE test_positivity

    PURE FUNCTION timescale(V, ODE)
        ! This subroutine should be a smart estimate for the size of the first timestep.
        ! For the moment I am using a not-so-smart estimate.
        REAL(8), DIMENSION(nvarode), INTENT(IN) :: V
        TYPE(T_ode_parameters),      INTENT(IN) :: ODE
        REAL(8), DIMENSION(nccc)                :: Ccc
        REAL(8), DIMENSION(ncc)                 :: Cc
        REAL(8)                                 :: timescale
        CALL system_linearization(Ccc, Cc, V, ODE)
        timescale = 1.0d-9
    END FUNCTION timescale

    ! -------------------------------------------------------------------------------------------- !

    PURE SUBROUTINE compute_source_jacobian_approximate(Jm, Q, t, ODE)
        ! Approximate computation of the jacobian of source for the fully coupled system
        REAL(8), DIMENSION(nvarode,nvarode), INTENT(OUT) :: Jm
        REAL(8),                             INTENT(IN)  :: t
        REAL(8), DIMENSION(nvarode),         INTENT(IN)  :: Q
        TYPE(T_ode_parameters),              INTENT(IN)  :: ODE
        REAL(8), DIMENSION(nvarode)                      :: Sl, Sr, Ql, Qr
        REAL(8)                                          :: delta
        INTEGER                                          :: k
        REAL(8), PARAMETER                               :: eps = 0.0d0
        REAL(8), PARAMETER                               :: delta0 = 1.0d-6
        INTEGER, PARAMETER                               :: finite_difference_mode = 2
        REAL(8), DIMENSION(nvarode,2)                    :: bounds
        IF (finite_difference_mode .eq. 0) THEN
            Ql = Q
            Qr = Q
            ! first differentiate wrt A
            delta = delta0
            DO k = 1, nvarode
                ! delta = delta0
                delta = delta0*max(1.0d0, abs(Q(k)))
                ! delta = delta0*max(1.0d-6, abs(Q(k)))
                Ql(k) = Q(k) - delta
                Qr(k) = Q(k) + delta
                Sl = source(t, Ql, ODE)
                Sr = source(t, Qr, ODE)
                Jm(:,k) = (Sr - Sl)/(Qr(k) - Ql(k) + 1.0d-100)
                Ql(k) = Q(k)
                Qr(k) = Q(k)
            END DO          
        ELSE IF (finite_difference_mode .eq. 1) THEN
            ! not tested!
            bounds(:,1) = -1.0d-300
            bounds(:,2) = 1.0d+300
            Ql = min(bounds(:,2), max(bounds(:,1), Q))
            Qr = Ql
            ! first differentiate wrt A
            delta = delta0
            DO k = 1, nvarode
                ! delta = delta0
                delta = delta0*max(1.0d0, abs(Q(k)))
                ! delta = delta0*max(1.0d-6, abs(Q(k)))
                Ql(k) = max(Q(k) - delta, bounds(k,1))
                Qr(k) = min(Q(k) + delta, bounds(k,2))
                Sl = source(t, Ql, ODE)
                Sr = source(t, Qr, ODE)
                Jm(:,k) = (Sr - Sl)/(Qr(k) - Ql(k) + 1.0d-100)
                Ql(k) = Q(k)
                Qr(k) = Q(k)
            END DO   
        ELSE IF (finite_difference_mode .eq. 2) THEN
            CALL compute_source_jacobian_approximate_gprdr(Jm, Q, t, ODE)
        END IF
    END SUBROUTINE compute_source_jacobian_approximate

    PURE SUBROUTINE linearized_system_solution(V, t, V0, Ccc, Cc, ODE)
        ! This subroutine computes the exact solution of the linearized ODE
        REAL(8), DIMENSION(nvarode), INTENT(OUT) :: V
        REAL(8),                     INTENT(IN)  :: t
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V0
        REAL(8), DIMENSION(nccc),    INTENT(IN)  :: Ccc
        REAL(8), DIMENSION(ncc),     INTENT(IN)  :: Cc
        TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
        IF (ODE%jacobian_mode .eq. 0) THEN
            ! compute a numerical approximation of the jacobian of the full system
            CALL linearized_system_solution_auto(V, t, V0, Ccc, Cc, ODE)
        ELSE IF (ODE%jacobian_mode .eq. 19) THEN
            ! Compute the jacobian for the strain relaxation part of the system analitically
            ! while decoupling the damage parameter evolution equation. The picard iterations
            ! will provide the coupling for the two subsystems
            ! CALL linearized_system_solution_gprdr_split_1_333(V, t, V0, Ccc, Cc, ODE)
            CALL linearized_system_solution_gprdr_split_1_9(V, t, V0, Ccc, Cc, ODE)
        ELSE IF (ODE%jacobian_mode .eq. 300) THEN
            ! Switch between block diagonal jacobian and better approximations (fully coupled)
            CALL linearized_system_solution_gprdr_split_adaptive_jacobian(V, t, V0, Ccc, Cc, ODE, .true.)
        ELSE IF (ODE%jacobian_mode .eq. 319) THEN
            ! Switch between block diagonal jacobian and better approximations (1-9 split)
            CALL linearized_system_solution_gprdr_split_adaptive_jacobian(V, t, V0, Ccc, Cc, ODE, .false.)
        ELSE
            ! compute a numerical approximation of the jacobian of the full system
            CALL linearized_system_solution_auto(V, t, V0, Ccc, Cc, ODE)
        END IF
    END SUBROUTINE linearized_system_solution

    ! -------------------------------------------------------------------------------------------- !

    PURE SUBROUTINE linearized_system_solution_auto(V, t, V0, Ccc, Cc, ODE)
        ! This subroutine computes the exact solution of the linearized ODE, when
        ! the system jacobian is approximated by finite differences.
        REAL(8), DIMENSION(nvarode), INTENT(OUT) :: V
        REAL(8),                     INTENT(IN)  :: t
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V0
        REAL(8), DIMENSION(nccc),    INTENT(IN)  :: Ccc
        REAL(8), DIMENSION(ncc),     INTENT(IN)  :: Cc
        TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
        REAL(8), DIMENSION(nvarode)              :: Vm, S0, Ve, S0_reduced, V0_reduced
        INTEGER, DIMENSION(nvarode)              :: mask
        REAL(8), DIMENSION(nvarode,nvarode)      :: Jm, Jm_reduced
        INTEGER                                  :: nv2
        REAL(8), PARAMETER                       :: eps = 1.0d-14 ! let us try 1.0d-14! maybe to be set to 1.0d-12, not more!
        ! here we extract the state vector for linearization
        Vm = Ccc(1:nvarode)
        ! note: in the following line we use t instead of tm (not important, the system is autonomous)
        CALL compute_source_jacobian_approximate(Jm, Vm, t, ODE)
        S0 = source(t, Vm, ODE) ! we use t instead of tm (not important, the system is autonomous)
        S0 = S0 - matmul(Jm, Vm)
        CALL build_jacobian_matrix_mask(mask, nv2, Jm, eps)
        CALL build_reduced_jacobian_matrix(Jm_reduced, S0_reduced, V0_reduced, Jm, S0, V0, mask)
        CALL compute_first_order_linear_ode_solution(Ve, Jm_reduced, V0_reduced, S0_reduced, nv2, t)
        CALL select_solution_from_mask(V, V0, S0, Ve, mask, t)
    END SUBROUTINE linearized_system_solution_auto

    ! -------------------------------------------------------------------------------------------- !











    ! -------------------------------------------------------------------------------------------- !
    ! ---- ODE-specific helper functions --------------------------------------------------------- !
    ! -------------------------------------------------------------------------------------------- !

    PURE SUBROUTINE linearized_system_solution_gprdr_split_adaptive_jacobian(V, t, V0, Ccc, Cc, ODE, fully_coupled_jacobian)
        ! adaptively chooses between a linearized_system_solution_gprdr_split_1_9_block_diagonal
        ! when AT A \simeq A AT and linearized_system_solution_gprdr_split_1_9 (or fully coupled) otherwise
        
        REAL(8), DIMENSION(nvarode), INTENT(OUT) :: V
        REAL(8),                     INTENT(IN)  :: t
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V0
        REAL(8), DIMENSION(nccc),    INTENT(IN)  :: Ccc
        REAL(8), DIMENSION(ncc),     INTENT(IN)  :: Cc
        LOGICAL,                     INTENT(IN)  :: fully_coupled_jacobian
        TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
        REAL(8), PARAMETER                       :: max_rel_off_diagonal = 1.0d-3 ! 1.0d-6
        REAL(8)                                  :: max_A_off_diagonal
        REAL(8)                                  :: Aeps(3,3)                   ! 2D test
        max_A_off_diagonal = maxval(abs([V0(3), V0(4), V0(5), V0(7), V0(8), V0(9)]))
        max_A_off_diagonal = MaxAeps(V0)
        IF (max_A_off_diagonal .lt. 1.e-5) THEN
        !IF (max_A_off_diagonal .lt. max_rel_off_diagonal*(abs(V0(2)) + abs(V0(6)) + abs(V0(10)))) THEN
            ! here I use a very safe criterion for allowing the simplification to be considered acceptable (maxval)
            CALL linearized_system_solution_gprdr_split_1_9_block_diagonal(V, t, V0, Ccc, Cc, ODE)
        ELSE
            IF (fully_coupled_jacobian) THEN
                CALL linearized_system_solution_auto(V, t, V0, Ccc, Cc, ODE) 
            ELSE
                CALL linearized_system_solution_gprdr_split_1_9(V, t, V0, Ccc, Cc, ODE)
            END IF
        END IF        
    END SUBROUTINE linearized_system_solution_gprdr_split_adaptive_jacobian

    PURE FUNCTION MaxAeps(V)
        IMPLICIT NONE
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V
        REAL(8), DIMENSION(3,3)                  :: Aeps
        REAL(8)                                  :: MaxAeps
        REAL(8) :: t1,t2,t3,t4,t5,t12,t21,t25,t26,t27,t35,t40,T13
        
        t1 = V(3) ** 2
        t2 = V(4) ** 2
        t3 = V(5) ** 2
        t4 = V(8) ** 2
        t5 = t1 + t2 - t3 - t4
        t13 = -V(2) * V(3) + V(2) * V(5) + V(3) * V(6) + V(4) * V(7) - V(5) * V(6) - V(8) * V(9)
        t21 = -V(2) * V(4) + V(2) * V(8) + V(3) * V(9) + V(4) * V(10) - V(5) * V(7) - V(8) * V(10)
        t25 = V(7) ** 2
        t26 = V(9) ** 2
        t27 = -t1 + t3 + t25 - t26
        t35 = -V(3) * V(4) + V(5) * V(8) - V(6) * V(7) + V(6) * V(9) + V(7) * V(10) - V(9) * V(10)
        t40 = -t2 - t25 + t4 + t26
        Aeps(1,1) = V(3) * t13 + V(4) * t21 + V(2) * t5
        Aeps(1,2) = V(2) * t13 + V(3) * t27 + V(4) * t35
        Aeps(1,3) = V(2) * t21 + V(3) * t35 + V(4) * t40
        Aeps(2,1) = V(6) * t13 + V(7) * t21 + V(5) * t5
        Aeps(2,2) = V(5) * t13 + V(6) * t27 + V(7) * t35
        Aeps(2,3) = V(5) * t21 + V(6) * t35 + V(7) * t40
        Aeps(3,1) = V(9) * t13 + V(10) * t21 + V(8) * t5
        Aeps(3,2) = V(8) * t13 + V(9) * t27 + V(10) * t35
        Aeps(3,3) = V(8) * t21 + V(9) * t35 + V(10) * t40
        MaxAeps=MAXVAL(Aeps)
    END FUNCTION MaxAeps
    
    PURE SUBROUTINE linearized_system_solution_gprdr_split_1_9(V, t, V0, Ccc, Cc, ODE)
        ! This subroutine computes the exact solution of the linearized ODE, when
        ! we consider separate subsystems to be coupled through picard iterations
        ! and adaptive timestepping.
        REAL(8), DIMENSION(nvarode), INTENT(OUT) :: V
        REAL(8),                     INTENT(IN)  :: t
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V0
        REAL(8), DIMENSION(nccc),    INTENT(IN)  :: Ccc
        REAL(8), DIMENSION(ncc),     INTENT(IN)  :: Cc
        TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
        REAL(8), DIMENSION(3,3)                  :: Am
        REAL(8), DIMENSION(9,9)                  :: Jm
        REAL(8), DIMENSION(9)                    :: S0, Ve
        INTEGER, DIMENSION(9)                    :: mask
        REAL(8), PARAMETER                       :: eps = 1.0d-8 ! here eps must be chosen slightly larger
        REAL(8), PARAMETER                       :: eps_c2 = 1.0d-12
        REAL(8), DIMENSION(9)                    :: S0_reduced, V0_reduced
        REAL(8), DIMENSION(9,9)                  :: Jm_reduced
        INTEGER                                  :: nv2
        REAL(8)                                  :: alpha1, alphar, itaum
        ! Compute the jacobian for the strain relaxation part of the system analitically
        ! while decoupling the damage parameter evolution equation. The picard iterations
        ! will provide the coupling for the two subsystems
        Am = build_by_rows(Ccc(2:10))
        alpha1 = sign(1.0d0, Ccc(11))*Ccc(11)**2
        alphar = sign(1.0d0, Ccc(12))*Ccc(12)**2
        itaum  = Ccc(13)**2
        CALL compute_source_jacobian_A(Jm, Am, itaum)
        ! CALL compute_source_jacobian_A_approximate(Jm, Am, itaum, ODE)
        CALL evaluate_source_A(S0, Am, itaum, ODE)
        S0 = S0 - matmul(Jm, flatten_rows(Am))
        ! build a mask for empty rows of the jacobian matrix
        CALL build_jacobian_matrix_mask(mask, nv2, Jm, eps)
        CALL build_reduced_jacobian_matrix(Jm_reduced, S0_reduced, V0_reduced, Jm, S0, V0(2:10), mask)
        CALL compute_first_order_linear_ode_solution(Ve, Jm_reduced, V0_reduced, S0_reduced, nv2, t)
        CALL select_solution_from_mask(V(2:10), V0(2:10), S0, Ve, mask, t)
        ! compute the solution to the damage parameter equation
        IF (V0(1) .le. 1.0d0 - eps_c2) THEN 
            V(1) = (V0(1) + alphar)*exp(alpha1*t) - alphar
        ELSE
            V(1) = 1.0d0
        END IF
#ifdef EQNTYPEC99
        ! Add the source terms for the thermal impulse
        !V(11:13) = V0(11:13)! No thermal contribution in the ODESolver ! (V0(11:13) + Ccc(16:18))*exp(Ccc(15)*t) - Ccc(16:18)
#endif
    END SUBROUTINE linearized_system_solution_gprdr_split_1_9

    PURE SUBROUTINE linearized_system_solution_gprdr_split_1_9_block_diagonal(V, t, V0, Ccc, Cc, ODE)
        ! This subroutine computes the exact solution of the linearized ODE, when
        ! we consider separate subsystems to be coupled through picard iterations
        ! and adaptive timestepping.
        REAL(8), DIMENSION(nvarode), INTENT(OUT) :: V
        REAL(8),                     INTENT(IN)  :: t
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V0
        REAL(8), DIMENSION(nccc),    INTENT(IN)  :: Ccc
        REAL(8), DIMENSION(ncc),     INTENT(IN)  :: Cc
        TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
        REAL(8), DIMENSION(3,3)                  :: Am
        REAL(8), DIMENSION(3,3,3)                :: Jm
        REAL(8), DIMENSION(9)                    :: S0
        REAL(8), DIMENSION(9)                    :: Amp
        REAL(8), DIMENSION(9)                    :: V0p
        INTEGER, DIMENSION(3)                    :: mask
        REAL(8), PARAMETER                       :: eps = 1.0d-8 ! here eps must be chosen slightly larger
        REAL(8), PARAMETER                       :: eps_c2 = 1.0d-12
        REAL(8), DIMENSION(3)                    :: S0_reduced, V0_reduced, Ve
        REAL(8), DIMENSION(3,3)                  :: Jm_reduced
        INTEGER                                  :: nv2, i, j
        REAL(8)                                  :: alpha1, alphar, itaum
        LOGICAL, PARAMETER                       :: use_faster_3by3_exponential = .true.
        ! Compute the jacobian for the strain relaxation part of the system analitically
        ! while decoupling the damage parameter evolution equation. The picard iterations
        ! will provide the coupling for the two subsystems
        Amp = Ccc(2:10)
        V0p = V0(2:10)
        Am = build_by_rows(Ccc(2:10))
        alpha1 = sign(1.0d0, Ccc(11))*Ccc(11)**2
        alphar = sign(1.0d0, Ccc(12))*Ccc(12)**2
        itaum  = Ccc(13)**2
        CALL compute_source_jacobian_A_block_diagonal(Jm, Am, itaum)
        CALL evaluate_source_A(S0, Am, itaum, ODE)
        CALL pack_vector_diagfirst(S0)
        CALL pack_vector_diagfirst(Amp)
        CALL pack_vector_diagfirst(V0p)
        ! solve the three simplified_subsystem
        IF (use_faster_3by3_exponential) THEN
            DO j = 0, 2       
                i = j*3 + 1     
                S0(i:i+2) = S0(i:i+2) - matmul(Jm(:,:,j+1), Amp(i:i+2))
                CALL fast33_build_jacobian_matrix_mask(mask, nv2, Jm(:,:,j+1), eps)
                CALL fast33_build_reduced_jacobian_matrix(Jm_reduced, S0_reduced, V0_reduced, Jm(:,:,j+1), S0(i:i+2), V0p(i:i+2), mask)
                CALL fast33_compute_first_order_linear_ode_solution(Ve, Jm_reduced, V0_reduced, S0_reduced, nv2, t)
                CALL fast33_select_solution_from_mask(V(i+1:i+3), V0p(i:i+2), S0(i:i+2), Ve, mask, t)
            END DO
        ELSE
            DO j = 0, 2       
                i = j*3 + 1     
                S0(i:i+2) = S0(i:i+2) - matmul(Jm(:,:,j+1), Amp(i:i+2))
                CALL build_jacobian_matrix_mask(mask, nv2, Jm(:,:,j+1), eps)
                CALL build_reduced_jacobian_matrix(Jm_reduced, S0_reduced, V0_reduced, Jm(:,:,j+1), S0(i:i+2), V0p(i:i+2), mask)
                CALL compute_first_order_linear_ode_solution(Ve, Jm_reduced, V0_reduced, S0_reduced, nv2, t)
                CALL select_solution_from_mask(V(i+1:i+3), V0p(i:i+2), S0(i:i+2), Ve, mask, t)
            END DO
        END IF

        CALL unpack_vector_diagfirst(V(2:10))
        !
        ! compute the solution to the damage parameter equation
        IF (V0(1) .le. 1.0d0 - eps_c2) THEN 
            V(1) = (V0(1) + alphar)*exp(alpha1*t) - alphar
        ELSE
            V(1) = 1.0d0
        END IF
#ifdef EQNTYPEC99
        ! Add the source terms for the thermal impulse
        !V(11:13) =V0(11:13)! No thermal contribuion in the ode solver ! (V0(11:13) + Ccc(16:18))*exp(Ccc(15)*t) - Ccc(16:18)
#endif
    END SUBROUTINE linearized_system_solution_gprdr_split_1_9_block_diagonal

    PURE SUBROUTINE pack_vector_diagfirst(A)
        ! This subroutine sorts the vector of nine components of A in such a way that
        ! the diagonal components are the first three and are followed by the six off-diagonal elements.
        REAL(8), DIMENSION(9), INTENT(INOUT) :: A
        A(2:9) = [A(5), A(9), A(2), A(3), A(4), A(6), A(7), A(8)]
    END SUBROUTINE pack_vector_diagfirst

    PURE SUBROUTINE unpack_vector_diagfirst(A)
        ! This subroutine is the inverse of pack_vector
        REAL(8), DIMENSION(9), INTENT(INOUT) :: A
        A(2:9) = [A(4), A(5), A(6), A(2), A(7), A(8), A(9), A(3)]
    END SUBROUTINE unpack_vector_diagfirst

    PURE SUBROUTINE compute_source_jacobian_approximate_gprdr(Jm, Q, t, ODE)
        ! Approximate computation of the jacobian of source for the fully coupled system
        REAL(8), DIMENSION(nvarode,nvarode), INTENT(OUT) :: Jm
        REAL(8),                             INTENT(IN)  :: t
        REAL(8), DIMENSION(nvarode),         INTENT(IN)  :: Q
        TYPE(T_ode_parameters),              INTENT(IN)  :: ODE
        REAL(8), DIMENSION(nvarode)                      :: Sl, Sr, Ql, Qr
        REAL(8)                                          :: delta
        INTEGER                                          :: k
        REAL(8), PARAMETER                               :: eps = 1.0d-10
        REAL(8), PARAMETER                               :: delta0 = 1.0d-6
        ! REAL(8), PARAMETER                               :: eps_c2 = 1.0d-18
        Ql = Q
        Qr = Q
        ! first differentiate wrt A
        delta = delta0
        DO k = 2, 10
            ! delta = delta0
            delta = delta0*max(1.0d0, abs(Q(k)))
            ! delta = delta0*max(1.0d-6, abs(Q(k)))
            Ql(k) = Q(k) - delta
            Qr(k) = Q(k) + delta
            Sl = source(t, Ql, ODE)
            Sr = source(t, Qr, ODE)
            Jm(:,k) = (Sr - Sl)/(Qr(k) - Ql(k) + 1.0d-100)
            Ql(k) = Q(k)
            Qr(k) = Q(k)
        END DO
        ! then differentiate wrt c
        IF (1.0d0 - Q(1) .gt. eps) THEN
            ! always special care for differentiating wrt c near zero
            delta = max(1.0d-10, min(delta0, 0.5d0*Q(1), 0.5d0*(1.0d0 - Q(1)))) !!!!
            ! differentiate wrt c
            Ql(1) = max(Q(1) - delta, 0.0d0)
            Qr(1) = min(Q(1) + delta, 1.0d0)
            Ql(1) = min(Ql(1), 1.0d0)
            Qr(1) = min(Qr(1), 1.0d0)
            Sl = source(t, Ql, ODE)
            Sr = source(t, Qr, ODE)
            Jm(:,1) = (Sr - Sl)/(Qr(1) - Ql(1) + 1.0d-100)
            Ql(1) = Q(1)
            Qr(1) = Q(1)
            Ql(1) = min(Ql(1), 1.0d0)
            Qr(1) = min(Qr(1), 1.0d0)
        ELSE
            Jm(1,:) = 0.0d0 !!! ok
            Jm(:,1) = 0.0d0
            Jm(1,1) = -1.0d0 
            ! Jm(1,1) = 1.0d0 - Q(1) 
            ! Jm(1,1) = 1.0d0 - Q(1) 
            ! Jm(1,2:) = 0.0d0
        END IF
        ! IF (Q(1) .gt. 1.0d0 - eps_c2) THEN 
        !     Jm(1,:) = 0.0d0
        !     Jm(1,1) = -1.0d0
        ! END IF
    END SUBROUTINE compute_source_jacobian_approximate_gprdr

    PURE SUBROUTINE system_linearization_gprdr(Ccc, Cc, V, ODE)
        USE SpecificVarEqn99, only : compute_l_m_dedc, compute_total_stress, compute_temperature
        ! This subroutine computes all or some of the linearization coefficients: 
        ! with numerical differentiation of the jacobian matrix, one could store 
        ! all the coefficients of the jacobian matrix in Ccc, but currently we 
        ! prefer to just compute some indicator variables for adaptive timestepping 
        ! and call this subroutine when evaluating the system source. 
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V
        REAL(8), DIMENSION(nccc),    INTENT(OUT) :: Ccc
        REAL(8), DIMENSION(ncc),     INTENT(OUT) :: Cc
        TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
        REAL(8)                                  :: Temp
        REAL(8)                                  :: tau1, tau2, Y, lam, mu, itaum, dEdc,Ydev,Yp
        REAL(8)                                  :: c1, c2, detA, alpha1, alpha0, alphar
        REAL(8), DIMENSION(3,3)                  :: A, G, devG
        REAL(8), PARAMETER                       :: floor_c2 = 1.0d-20
        REAL(8), PARAMETER                       :: detA_threshold = 2.0d0
        REAL(8), DIMENSION(3,3), PARAMETER       :: ID = reshape(&
            [1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0], [3,3])
        c1 = V(1)
        c2 = 1.0d0 - c1
        A = build_by_rows(V(2:10))
        G = matmul(transpose(A), A)
        devG  = G - trace(G)/3.0d0*ID
        detA = determinant_lu(A)
        CALL compute_l_m_dedc(lam, mu, dEdc, c1, devG, detA, ODE)
        CALL compute_total_stress(Y,Ydev,Yp, lam, mu, G, devG, detA, ODE)
        ! Compute the local relaxation times
        tau1 = ODE%tau0*exp(min(ODE%alpha1 - ODE%beta1*c2*Y, 690.0d0))  
        tau2 = ODE%tau0*exp(min(ODE%alpha2 - ODE%beta2*c1*Y, 690.0d0)) 
        ! Compute the relaxation time for the mixture; itaum = c2/mu1*tau1 + c1/(mu2*tau2)
        itaum = (c1*ODE%mu1*tau1 + c2*ODE%mu2*tau2)/(ODE%mu1*tau1*ODE%mu2*tau2 + 1.0d-200)
        itaum = min(itaum, 1.0d0/ODE%taumin)
        itaum = itaum*min(1.0d0, (detA_threshold/abs(detA))**(6.0d0/3.0d0))
        IF (c2 .ge. floor_c2) THEN
            alpha1 = -dEdc*ODE%KK*c2*(c2*abs(Y/ODE%Y0)**ODE%aexp + c1*abs(Y/ODE%Y1))
        ELSE
            alpha1 = 0.0d0
        END IF
        alpha0 = alpha1*ODE%xieps 
        alphar = alpha0*alpha1/(alpha1**2 + 1.0d-30)
        ! Ccc: 
        ! 1-10. strain-rate tensor components
        ! 11. first  coefficient of the solution of the damage parameter evolution (with square root scaling)
        ! 12. second coefficient of the solution of the damage parameter evolution (with square root scaling)
        ! 13. characteristic timescale for the strain rate relaxation
        ! 14. characteristic timescale for the strain rate relaxation (including A scaling)
        Ccc(1:10) = V(1:10)
        Ccc(11) = sign(1.0d0, alpha1)*abs(alpha1)**0.5d0
        Ccc(12) = sign(1.0d0, alphar)*abs(alphar)**0.5d0
        Ccc(13) = itaum**0.5d0
        Ccc(14) = (itaum*detA**(6.0d0/3.0d0))**0.5d0
        
#ifdef EQNTYPEC99
        ! Add the constcoeff for the thermal impulse
        !call compute_temperature(Temp, detA, lam, mu,ODE%rho0,ODE%gamma,ODE%cv, ODE%p0,ODE%s, ODE%EOS, ODE%T0) ! (!) rho0 here is the local rho (!) 
        !alpha1      = ODE%rho0*Temp/ODE%tau2
        !Ccc(15)     = alpha1 ! Coefficient of Jk
        !Ccc(16:18)  = ODE%F0J(:)*alpha1/(alpha1**2 + 1.0d-30) ! Coefficients given by the forcing terms         
#endif
        
        
        cc=0.0d0
        ! cc(1) = abs(Y/ODE%Y0)**ODE%aexp
        ! cc(2) = dedc
        ! for use with linearized_system_solution_auto put V in Ccc(1:nvarode) (or in whichever other 
        ! position you prefer but remember to specify the position in linearized_system_solution_auto)
    END SUBROUTINE system_linearization_gprdr

!    PURE SUBROUTINE compute_l_m_dedc(lam, mu, dEdc, c, devG, detA, ODE)
!        ! Elasticity constant homogeneization and computation of the derivative of 
!        ! energy with respect to the damage parameter.
!        REAL(8),                 INTENT(OUT) :: lam, mu, dEdc
!        REAL(8), DIMENSION(3,3), INTENT(IN)  :: devG
!        TYPE(T_ode_parameters),  INTENT(IN)  :: ODE
!        REAL(8),                 INTENT(IN)  :: c, detA
!        REAL(8)                              :: mus, Ks, K, dM, dK
!        REAL(8), PARAMETER                   :: zero = 1.0d-100
!        mus = c*ODE%mu1 + (1.0d0 - c)*ODE%mu2
!        Ks  = c*ODE%K1  + (1.0d0 - c)*ODE%K2
!        mu  = ODE%mu1*ODE%mu2/(mus + zero)
!        K   = ODE%K1*ODE%K2/(Ks + zero)
!        lam = K - 2.0d0/3.0d0*mu
!        dM  = (ODE%mu1 - ODE%mu2)*mu/mus
!        dK  = (ODE%K1 - ODE%K2)*K/Ks
!        dEdc = -(0.5d0*dK*(1.0d0 - detA)**2 + 0.25d0*dM*sum(devG**2))/ODE%rho0 
!    END SUBROUTINE compute_l_m_dedc


    
! ---- auxiliary subroutines for the evaluation of exact jacobian of S wrt A ----------------- !

    PURE SUBROUTINE evaluate_source_A(S, A, itaum, ODE)
        ! Source of the A tensor
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: A
        REAL(8),                 INTENT(IN)  :: itaum
        TYPE(T_ode_parameters),  INTENT(IN)  :: ODE
        REAL(8), DIMENSION(9),   INTENT(OUT) :: S
        REAL(8), DIMENSION(3,3)              :: devG
        REAL(8)                              :: detA
        REAL(8), DIMENSION(3,3), PARAMETER   :: ID = reshape(&
            [1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0], [3,3])
        devG = matmul(transpose(A), A)
        devG = devG - ID*((devG(1,1) + devG(2,2) + devG(3,3))/3.0d0)
        detA = abs(determinant_lu(A))
        S = -3.0d0*itaum*detA**(5.0d0/3.0d0)*flatten_rows(matmul(A, devG)) + ODE%F0
    END SUBROUTINE evaluate_source_A

    PURE SUBROUTINE compute_source_jacobian_A(Jm, A, itaum)
        ! Analitically computed jacobian matrix of the strain-rate subsystem (with fixed taum).
        REAL(8), DIMENSION(9,9), INTENT(OUT) :: Jm
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: A
        REAL(8),                 INTENT(IN)  :: itaum
        REAL(8), DIMENSION(3,3)              :: devG, AT, G, H, crA
        REAL(8), DIMENSION(3,3,3,3)          :: B
        REAL(8), DIMENSION(9)                :: U1, V1, U2, V2
        REAL(8)                              :: detA
        INTEGER                              :: i, j
        REAL(8), DIMENSION(3,3), PARAMETER   :: ID = reshape(&
            [1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0], [3,3])
        AT = transpose(A)
        G = matmul(AT, A)
        H = matmul(A, AT)
        devG = G - trace(G)/3.0d0*ID
        detA = abs(determinant_lu(A))
        crA(1,:) = cross_product(A(2,:), A(3,:))
        crA(2,:) = cross_product(A(3,:), A(1,:))
        crA(3,:) = cross_product(A(1,:), A(2,:))
        U1 = flatten_rows(matmul(A, devG))
        V1 = - 5.0d0*detA**(2.0d0/3.0d0)/3.0d0*flatten_rows(crA)
        U2 = flatten_rows(A)
        V2 = 2.0d0*detA**(5.0d0/3.0d0)/3.0d0*U2
        DO j = 1, 3
            DO i = 1, 3
                B(:,:,i,j) = diadic_product(AT(:,j), AT(:,i)) + ID(i,j)*devG + H(i,j)*ID
            END DO
        END DO
        B = -detA**(5.0d0/3.0d0)*B
        Jm = 3.0d0*itaum*(halven_rank42_rows(B) + diadic_product(U1, V1) + diadic_product(U2, V2))
    END SUBROUTINE compute_source_jacobian_A

    PURE SUBROUTINE compute_source_jacobian_A_block_diagonal(Jm, A, itaum)
        ! This is a simplified approximate computation of the same object given by 
        ! compute_source_jacobian_A, that is the jacobian matrix of the strain-rate 
        ! subsystem (with fixed taum).
        REAL(8), DIMENSION(3,3,3), INTENT(OUT) :: Jm
        REAL(8), DIMENSION(3,3),   INTENT(IN)  :: A
        REAL(8),                   INTENT(IN)  :: itaum
        REAL(8), DIMENSION(9,9)                :: Jmfull
        CALL compute_source_jacobian_A(Jmfull, A, itaum)
        Jm(1,1,1) = Jmfull(1,1)
        Jm(2,1,1) = Jmfull(5,1)
        Jm(3,1,1) = Jmfull(9,1)
        Jm(1,2,1) = Jmfull(1,5)
        Jm(2,2,1) = Jmfull(5,5)
        Jm(3,2,1) = Jmfull(9,5)
        Jm(1,3,1) = Jmfull(1,9)
        Jm(2,3,1) = Jmfull(5,9)
        Jm(3,3,1) = Jmfull(9,9)
        Jm(:,:,2) = Jmfull(2:4,2:4)
        Jm(:,:,3) = Jmfull(6:8,6:8)
        ! account for the two jacobian entries that have been suppressed
        ! the following two lines may be removed in case they cause issues
        Jm(2,2,2) = 2.0d0*Jm(2,2,2)
        Jm(2,2,3) = 2.0d0*Jm(2,2,3)
    END SUBROUTINE compute_source_jacobian_A_block_diagonal

    PURE SUBROUTINE compute_source_jacobian_A_approximate(Jm, A, itaum, ODE)
        ! Finite difference approximation for the jacobian of the strain-rate subsystem (with fixed taum).
        REAL(8), DIMENSION(9,9), INTENT(OUT) :: Jm
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: A
        REAL(8),                 INTENT(IN)  :: itaum
        TYPE(T_ode_parameters),  INTENT(IN)  :: ODE
        REAL(8), DIMENSION(3,3)              :: Al, Ar
        REAL(8), DIMENSION(9)                :: Sl, Sr
        REAL(8), PARAMETER                   :: eps = 1.0d-6
        INTEGER                              :: i, j, k
        k = 0
        DO i = 1, 3
            DO j = 1, 3
                k = k + 1
                Al = A
                Ar = A
                Al(i,j) = Al(i,j) - eps
                Ar(i,j) = Ar(i,j) + eps
                CALL evaluate_source_A(Sl, Al, itaum, ODE)
                CALL evaluate_source_A(Sr, Ar, itaum, ODE)
                Jm(:,k) = (Sr - Sl)/(2.0d0*eps)
            END DO
        END DO
    END SUBROUTINE compute_source_jacobian_A_approximate

    ! ---- for output ---------------------------------------------------------------------------- !

    !PURE SUBROUTINE compute_auxiliary_variables(V_aux, V, ODE)
    !    ! Output variables
    !    REAL(8), DIMENSION(naux),    INTENT(OUT) :: V_aux 
    !    REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V
    !    REAL(8), DIMENSION(nccc)                 :: Ccc
    !    REAL(8), DIMENSION(ncc)                  :: Cc
    !    TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
    !    CALL system_linearization(Ccc, Cc, V, ODE)
    !    V_aux(1) = V(1)   ! c
    !    V_aux(2) = V(2)   ! a11
    !    V_aux(3) = V(6)   ! a22
    !    V_aux(4) = V(10)  ! a33
    !    V_aux(5) = V(3)   ! a12
    !    V_aux(6) = Ccc(11) ! alpha1
    !    V_aux(7) = Ccc(12) ! alphar
    !    V_aux(8) = Ccc(13) ! itaum
    !    ! V_aux(7) = Cc(1) 
    !END SUBROUTINE compute_auxiliary_variables

END MODULE expintegrator_ode
#endif
#endif
