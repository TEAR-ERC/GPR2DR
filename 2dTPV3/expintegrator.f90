#ifndef HEADER_ODERSOLVER
#define HEADER_ODERSOLVER
#if defined(ODESOLVER)
!
!#include "MainVariables.f90"
!#include "expintegrator_type.f90"
!#include "expintegrator_ode.f90"
!#include "expintegrator_linalg.f90"

MODULE expintegrator
    use MainVariables, only : tEquation
    use expintegrator_type, only : T_odeexp_settings,nvarode, nccc, ncc, T_ode_parameters
    USE expintegrator_ode, only: convert_parameters, linearized_system_solution, timescale, system_linearization, test_positivity
    USE expintegrator_linalg, only:
    ! USE linalg, only:
    ! USE expintegrator_ode_tools, only:
    USE, INTRINSIC :: ieee_arithmetic

    IMPLICIT NONE

    
    PRIVATE
    PUBLIC :: expintegrator_adaptive
    PUBLIC :: expintegrator_adaptive_full

CONTAINS

    PURE FUNCTION tol_fpi_default(delta_max)
        REAL(8), INTENT(IN) :: delta_max
        REAL(8)             :: tol_fpi_default
        tol_fpi_default = min(1.0d-4, max(1.0d-9, 1.0d-2*delta_max**2))
    END FUNCTION tol_fpi_default

    PURE SUBROUTINE solve_nonlinear_system(dt, V0, Ccc_tn, Cc_tn, ODE, SETTINGS, tol_fpi, &
        V, Ccc_tnp12, Cc_tnp12, Ccc_tnp1, Cc_tnp1, success, dt_next, niter, statusflag)
        REAL(8),                     INTENT(IN)  :: dt
        REAL(8), DIMENSION(nvarode), INTENT(IN)  :: V0
        REAL(8), DIMENSION(nccc),    INTENT(IN)  :: Ccc_tn
        REAL(8), DIMENSION(ncc),     INTENT(IN)  :: Cc_tn
        TYPE(T_ode_parameters),      INTENT(IN)  :: ODE
        TYPE(T_odeexp_settings),     INTENT(IN)  :: SETTINGS
        REAL(8),                     INTENT(IN)  :: tol_fpi
        REAL(8), DIMENSION(nvarode), INTENT(OUT) :: V
        REAL(8), DIMENSION(nccc),    INTENT(OUT) :: Ccc_tnp12
        REAL(8), DIMENSION(ncc),     INTENT(OUT) :: Cc_tnp12
        REAL(8), DIMENSION(nccc),    INTENT(OUT) :: Ccc_tnp1
        REAL(8), DIMENSION(ncc),     INTENT(OUT) :: Cc_tnp1
        LOGICAL,                     INTENT(OUT) :: success
        REAL(8),                     INTENT(OUT) :: dt_next
        INTEGER,                     INTENT(OUT) :: niter
        INTEGER,                     INTENT(OUT) :: statusflag
        REAL(8), DIMENSION(nvarode)              :: Vx, Vm
        REAL(8), DIMENSION(nccc)                 :: eps_ccc
        REAL(8)                                  :: eps_fpi, increment_ccfl, delta_V
        REAL(8)                                  :: delta, delta_max, decrement_accuracy, decrement_positivity, delta0
        INTEGER                                  :: maxiter_accept, maxiter_reject, posstatus, i
        delta_max            = SETTINGS%delta_max
        maxiter_accept       = SETTINGS%maxiter_accept
        maxiter_reject       = SETTINGS%maxiter_reject
        eps_fpi              = SETTINGS%eps_fpi
        eps_ccc              = SETTINGS%eps_ccc
        increment_ccfl       = SETTINGS%increment_ccfl
        decrement_accuracy   = SETTINGS%decrement_accuracy
        decrement_positivity = SETTINGS%decrement_positivity
        V = V0
        Vx = V0
        Cc_tnp1 = Cc_tn
        Ccc_tnp1 = Ccc_tn
        delta = 1.0d10
        success = .false.
        niter = maxiter_accept
        ! if too many iterations for convergence, probably the linearization would have been 
        !     rejected at convergence, so reduce the timestep size
        dt_next = dt*0.5d0 ! 0.5 seems very good (this is used if past maxiter_reject)
        statusflag = 0
        DO i = 1, maxiter_reject
            ! 1. linearize the system at the average state for the timestep
            Vm = 0.5d0*(V0 + V)
            CALL system_linearization(Ccc_tnp12, Cc_tnp12, Vm, ODE)
            ! 2. save old iteration and solve linearized ODE for pressure and alpha
            Vx = V
            CALL linearized_system_solution(V, dt, V0, Ccc_tnp12, Cc_tnp12, ODE)
            ! 3. check positivity
            CALL test_positivity(V, ODE, posstatus) 
            IF (posstatus .gt. 0 .and. .not. any(ieee_is_nan(V)) .and. all(ieee_is_finite(V))) THEN
                ! 4. check if convergence criterion is satisfied
                delta_V = maxval(abs(V - Vx)/(abs(V) + abs(Vx) + eps_fpi))
                IF (delta_V .lt. tol_fpi .or. i .eq. maxiter_accept) THEN
                    ! 5. check the validity of the midpoint linearization (might in some cases not be strictly necessary)
                    delta0 = maxval(abs(Ccc_tnp12 - Ccc_tn)/(abs(Ccc_tnp12) + abs(Ccc_tn) + eps_ccc))
                    statusflag = maxloc(abs(Ccc_tnp12 - Ccc_tn)/(abs(Ccc_tnp12) + abs(Ccc_tn) + eps_ccc), 1)
                    ! 6. compute and validate linearisation for next timestep
                    CALL system_linearization(Ccc_tnp1, Cc_tnp1, V, ODE)
                    delta = max(delta0, maxval(abs(Ccc_tnp1 - Ccc_tn)/(abs(Ccc_tnp1) + abs(Ccc_tn) + eps_ccc)))
                    IF (delta .gt. delta0) THEN
                        statusflag = maxloc(abs(Ccc_tnp1 - Ccc_tn)/(abs(Ccc_tnp1) + abs(Ccc_tn) + eps_ccc), 1)
                    END IF
                    ! write(*,"(A, 1P 13 E27.16)") "", abs(Ccc_tnp12 - Ccc_tn)/(abs(Ccc_tnp12) + abs(Ccc_tn) + eps_ccc)
                    ! write(*,"(A, 1P 13 E27.16)") "", abs(Ccc_tnp1 - Ccc_tn)/(abs(Ccc_tnp1) + abs(Ccc_tn) + eps_ccc)
                    niter = i
                    IF (delta .lt. delta_max) THEN
                        ! if everything ok, choose next timestep size
                        dt_next = dt*(increment_ccfl*delta_max/(delta + 1.0d-14))
                        ! dt_next = min(dt_next, 0.01d0)
                        success = .true.
                        ! statusflag = 0
                    ELSE
                        ! if the linearization is not sufficiently accurate, reduce timestep size
                        dt_next = dt*decrement_accuracy ! choose, 1/3 maybe
                        ! statusflag = 1
                    END IF
                    EXIT
                END IF
            ELSE
                ! if positivity has been violated, reduce timestep size
                dt_next = dt*decrement_positivity ! 1/3 seems very good
                ! positivity is violated, then reduce timestep in the outer loop. I have never seen this happen.
                niter = i
                statusflag = -1
                EXIT
            END IF
        END DO
    END SUBROUTINE solve_nonlinear_system

    SUBROUTINE expintegrator_adaptive(V0, EQN, tend, SETTINGS, V, substep_count, total_iterations,ODE)
        REAL(8), DIMENSION(nvarode),    INTENT(IN)  :: V0
        TYPE(tEquation), INTENT(IN)  :: EQN
        REAL(8),                     INTENT(IN)  :: tend
        REAL(8), DIMENSION(nvarode),    INTENT(OUT) :: V
        TYPE(T_odeexp_settings),     INTENT(IN)  :: SETTINGS
        INTEGER,                     INTENT(OUT) :: substep_count
        INTEGER,                     INTENT(OUT) :: total_iterations
        TYPE(T_ode_parameters)                   :: ODE
        REAL(8), DIMENSION(nvarode)                 :: Vi, Vj
        REAL(8), DIMENSION(nccc)                 :: Ccc_tn, Ccc_tnp12, Ccc_tnp1
        REAL(8), DIMENSION(ncc)                  :: Cc_tn, Cc_tnp12, Cc_tnp1
        REAL(8)                                  :: tol_fpi, k_dt0, t, dt, dt_next
        INTEGER                                  :: max_timesteps, max_retry, i, j, niter, posstatus, statusflag
        LOGICAL :: last_step, try_last_step, success
        k_dt0 = 1.0d-9
        tol_fpi = tol_fpi_default(SETTINGS%delta_max)
        max_timesteps = 2000000
        max_retry = 20000
        ! initialize parameters and set initial conditions
        Vi = V0
        !CALL convert_parameters(V0, ODE, EQN)
        CALL test_positivity(V0, ODE, posstatus)
        IF (posstatus .lt. 0) THEN
            ! WRITE(*,"(A)") "invalid input data"
            ! DO i = 1, nvarode
            !     write(*,"(A, I3, 1PE27.16)") "", i, V0(i)
            ! END DO
            V = V0
            substep_count = -1
            total_iterations = -1
            RETURN
        END IF
        CALL system_linearization(Ccc_tn, Cc_tn, Vi, ODE)
        ! choose first timestep
        dt = k_dt0*timescale(V0, ODE)
        ! begin timestepping
        t = 0.0d0
        total_iterations = 0
        last_step = .false.
        try_last_step = .false.
        DO i = 1, max_timesteps
            IF (dt .lt. 1.0d-100) THEN
                ! WRITE(*,"(A)") "timestep too small in expintegrator!"
                V = Vi
                substep_count = -i
                total_iterations = -1 - abs(total_iterations)
                RETURN
            END IF
            IF (t + dt .ge. tend - 1.0d-14) THEN
                dt = tend - t
                try_last_step = .true.
            END IF
            ! retry step until successful
            DO j = 1, max_retry
                ! update from the fixed point iteration
                CALL solve_nonlinear_system(dt, Vi, Ccc_tn, Cc_tn, ODE, SETTINGS, tol_fpi, & 
                    Vj, Ccc_tnp12, Cc_tnp12, Ccc_tnp1, Cc_tnp1, success, dt_next, niter, statusflag)
                total_iterations = total_iterations + niter
                IF (success) THEN
                    ! advance time
                    t = t + dt
                    ! confirm solution
                    Vi = Vj
                    ! confirm next timestep value
                    ! dt = dt_next
                    IF (j .eq. 1) THEN
                        dt = dt_next
                    END IF
                    ! confirm constant coefficients
                    Ccc_tn = Ccc_tnp1
                    Cc_tn = Cc_tnp1
                    last_step = try_last_step
                    EXIT
                ELSE
                    ! retry with a smaller timestep
                    last_step = .false.
                    try_last_step = .false.
                    dt = dt_next
                END IF
            END DO
            IF (last_step) THEN
                EXIT
            END IF
        END DO
        IF (last_step) THEN
            substep_count = i
        ELSE
            substep_count = -i
            total_iterations = -1 - abs(total_iterations)
        END IF
        V = Vi
    END SUBROUTINE expintegrator_adaptive

    SUBROUTINE expintegrator_adaptive_full(V0, EQN, tend, SETTINGS, Nt, time, Vt, Ccc, Cc, itercount, statusflag)
        REAL(8), DIMENSION(nvarode),    INTENT(IN)  :: V0
        TYPE(tEquation), INTENT(IN)  :: EQN
        REAL(8),                     INTENT(IN)  :: tend
        TYPE(T_odeexp_settings),     INTENT(IN)  :: SETTINGS
        INTEGER,                     INTENT(OUT) :: nt
        TYPE(T_ode_parameters)                   :: ODE
        REAL(8), DIMENSION(nvarode)                 :: Vi, Vj
        REAL(8), DIMENSION(:,:), INTENT(OUT)     :: Vt, Cc, Ccc
        REAL(8), DIMENSION(:),   INTENT(OUT)     :: time
        INTEGER, DIMENSION(:),   INTENT(OUT)     :: itercount
        INTEGER, DIMENSION(:,:), INTENT(OUT)     :: statusflag
        REAL(8), DIMENSION(nccc)                 :: Ccc_tn, Ccc_tnp12, Ccc_tnp1
        REAL(8), DIMENSION(ncc)                  :: Cc_tn, Cc_tnp12, Cc_tnp1
        REAL(8)                                  :: tol_fpi, k_dt0, t, dt, dt_next
        INTEGER                                  :: max_timesteps, max_retry, i, j, niter, posstatus, statusflag_i
        LOGICAL :: last_step, try_last_step, success
        k_dt0 = 1.0d-9
        tol_fpi = tol_fpi_default(SETTINGS%delta_max)
        max_timesteps = size(time)
        max_retry = 20000
        ! initialize parameters and set initial conditions
        Vi = V0
        !CALL convert_parameters(V0, ODE, EQN)
        CALL test_positivity(V0, ODE, posstatus)
        IF (posstatus .lt. 0) THEN
            WRITE(*,"(A)") "invalid input data"
            DO i = 1, nvarode
                write(*,"(A, I3, 1PE27.16)") "", i, V0(i)
            END DO
            RETURN
        END IF
        CALL system_linearization(Ccc_tn, Cc_tn, Vi, ODE)
        ! choose first timestep
        dt = k_dt0*timescale(V0, ODE)
        ! begin timestepping
        t = 0.0d0
        try_last_step = .false.
        last_step = .false.
        time(1) = t
        Vt(:,1) = V0
        itercount = 0
        statusflag = 0
        DO i = 2, max_timesteps
            IF (dt .lt. 1.0d-100) THEN
                WRITE(*,"(A)") "timestep too small in expintegrator!"
                RETURN
            END IF
            IF (t + dt .ge. tend - 1.0d-14) THEN
                dt = tend - t
                try_last_step = .true.
            END IF
            ! retry step until successful
            DO j = 1, max_retry
                write(*,"(A, I6, 1P 2 E27.16)") "i, t, dt", i, t, dt
                ! update from the fixed point iteration
                CALL solve_nonlinear_system(&
                    dt, Vi, Ccc_tn, Cc_tn, ODE, SETTINGS, tol_fpi, Vj, Ccc_tnp12, Cc_tnp12, &
                    Ccc_tnp1, Cc_tnp1, success, dt_next, niter, statusflag_i)
                ! write(*,"(A, 1P 10 E27.16)") "", Vi
                itercount(i) = itercount(i) + niter
                statusflag(3,i) = statusflag_i
                IF (success) THEN
                    ! advance time
                    t = t + dt
                    ! confirm solution
                    Vi = Vj
                    time(i) = t
                    Vt(:,i) = Vi
                    Ccc(:,i) = Ccc_tnp12
                    Cc(:,i) = Cc_tnp12
                    ! confirm next timestep value
                    IF (j .eq. 1) THEN
                        dt = dt_next
                    END IF 
                    ! confirm constant coefficients
                    Ccc_tn = Ccc_tnp1
                    Cc_tn = Cc_tnp1
                    last_step = try_last_step
                    EXIT
                ELSE
                    ! retry with a smaller timestep
                    last_step = .false.
                    try_last_step = .false.
                    dt = dt_next
                    statusflag(1,i) = statusflag(1,i) + max(sign(1, statusflag_i), 0) ! accumulate rejections due to precision requirements
                    statusflag(2,i) = statusflag(2,i) + min(sign(1, statusflag_i), 0) ! accumulate rejections due to positivity violation
                END IF
            END DO
            IF (last_step) THEN
                EXIT
            END IF
        END DO
        Nt = min(i, max_timesteps)
    END SUBROUTINE expintegrator_adaptive_full

END MODULE expintegrator
    
#endif
#endif
