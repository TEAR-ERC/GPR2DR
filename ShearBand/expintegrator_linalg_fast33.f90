#ifndef HEADER_EXPINTEGRATOR_LINALG_FAST
#define HEADER_EXPINTEGRATOR_LINALG_FAST
#if defined(ODESOLVER) || defined(EQNTYPEC99) || defined(EQNTYPED99)
    
MODULE expintegrator_linalg_fast33

    IMPLICIT NONE

    PRIVATE
    
    PUBLIC  :: fast33_compute_matrix_exponential
    PUBLIC  :: fast33_build_jacobian_matrix_mask
    PUBLIC  :: fast33_build_reduced_jacobian_matrix
    PUBLIC  :: fast33_compute_first_order_linear_ode_solution
    PUBLIC  :: fast33_select_solution_from_mask

CONTAINS


! robust matrix exponential

    ! PURE FUNCTION matrix_norm1(A)
    !     REAL(8), DIMENSION(:,:), INTENT(IN) :: A
    !     REAL(8)                             :: matrix_norm1
    !     ! take the largest l1 norm of all columns (sum over column)
    !     matrix_norm1 = maxval(sum(abs(A), dim=1)) 
    ! END FUNCTION matrix_norm1

    ! PURE FUNCTION fast33_outerprod(a, b)
        ! REAL(8), DIMENSION(3), INTENT(IN)   :: a, b
        ! REAL(8), DIMENSION(3,3)             :: fast33_outerprod
        ! fast33_outerprod = spread(a, dim=2, ncopies=3)*spread(b, dim=1, ncopies=3)
    ! END FUNCTION fast33_outerprod

    ! PURE FUNCTION fast22_outerprod(a, b)
        ! REAL(8), DIMENSION(2), INTENT(IN)   :: a, b
        ! REAL(8), DIMENSION(2,2)             :: fast22_outerprod
        ! fast22_outerprod = spread(a, dim=2, ncopies=2)*spread(b, dim=1, ncopies=2)
    ! END FUNCTION fast22_outerprod

    PURE FUNCTION fast33_matrix_norm1_pow(A, p)
        REAL(8), DIMENSION(:,:), INTENT(IN)      :: A
        INTEGER,                 INTENT(IN)      :: p
        REAL(8), DIMENSION(3)                    :: v
        REAL(8), DIMENSION(3,3)                  :: AT
        REAL(8)                                  :: fast33_matrix_norm1_pow
        INTEGER                                  :: i
        v = 1.0d0
        AT = transpose(A)
        DO i = 1, p
            v = matmul(AT, v)
        END DO
        fast33_matrix_norm1_pow = maxval(v)
    END FUNCTION fast33_matrix_norm1_pow

    PURE FUNCTION fast22_matrix_norm1_pow(A, p)
        REAL(8), DIMENSION(:,:), INTENT(IN)      :: A
        INTEGER,                 INTENT(IN)      :: p
        REAL(8), DIMENSION(2)                    :: v
        REAL(8), DIMENSION(2,2)                  :: AT
        REAL(8)                                  :: fast22_matrix_norm1_pow
        INTEGER                                  :: i
        v = 1.0d0
        AT = transpose(A)
        DO i = 1, p
            v = matmul(AT, v)
        END DO
        fast22_matrix_norm1_pow = maxval(v)
    END FUNCTION fast22_matrix_norm1_pow

    PURE FUNCTION abs_c_inv(p)
        INTEGER, INTENT(IN) :: p
        INTEGER             :: r, i
        REAL(8)             :: abs_c_inv
        r = 1
        DO i = 1, p
            r = r*(p + i)
        END DO
        r = r**2*(2*p + 1)
        abs_c_inv = real(r, 8)
    END FUNCTION abs_c_inv

    ! ---------------------------------------------------------------------------------------------
    ! fast matrix exponential for 3 by 3 systems (with hardcoded matrix size)

    PURE SUBROUTINE fast33_compute_matrix_exponential(EA, A0)
        ! as the above, but with hardcoded 3 by 3 matrix size
        REAL(8), DIMENSION(3,3), INTENT(OUT) :: EA
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: A0
        REAL(8), DIMENSION(5),   PARAMETER   :: theta = [1.495585217958292d-2, &
            2.539398330063230d-1, 9.504178996162932d-1, 2.097847961257068d0, 4.25d0]
        INTEGER :: i, s, ludcmp_d
        INTEGER, PARAMETER :: N = 3
        REAL(8), DIMENSION(3,3) :: ID, A, U, V, Q
        REAL(8), DIMENSION(3,3,4) :: A_pow
        INTEGER, DIMENSION(3) :: ludcmp_indx
        REAL(8) :: d4, d6, d8, d10, eta1, eta2, eta3, eta4, eta5
        s = 0
        ID = 0.0d0
        DO i = 1, N
            ID(i,i) = 1.0d0
        END DO
        A = A0
        A_pow(:,:,1) = matmul(A, A)
        A_pow(:,:,2) = matmul(A_pow(:,:,1), A_pow(:,:,1))
        A_pow(:,:,3) = matmul(A_pow(:,:,2), A_pow(:,:,1))
        d4 = fast33_matrix_norm1(A_pow(:,:,2))**(1.0d0/4.0d0)
        d6 = fast33_matrix_norm1(A_pow(:,:,3))**(1.0d0/6.0d0)
        eta1 = max(d4, d6)
        eta2 = eta1
        IF (eta1 .le. theta(1) .and. fast33_ell(A, 3) .eq. 0) THEN
            CALL fast33_pade3(U, V, ID, A, A_pow)
        ELSE IF (eta2 .le. theta(2) .and. fast33_ell(A, 5) .eq. 0) THEN
            CALL fast33_pade5(U, V, ID, A, A_pow)
        ELSE
            A_pow(:,:,4) = matmul(A_pow(:,:,3), A_pow(:,:,1))
            d8 = fast33_matrix_norm1(A_pow(:,:,4))**(1.0d0/8.0d0)
            eta3 = max(d6, d8)
            IF (eta3 .le. theta(3) .and. fast33_ell(A, 7) .eq. 0) THEN
                CALL fast33_pade7(U, V, ID, A, A_pow)
            ELSE IF (eta3 .le. theta(4) .and. fast33_ell(A, 9) .eq. 0) THEN
                CALL fast33_pade9(U, V, ID, A, A_pow)
            ELSE
                ! since the degree 13 pade approximant does not use A^8, we can overwrite A^10 on A^8
                A_pow(:,:,4) = matmul(A_pow(:,:,2), A_pow(:,:,3))
                d10 = fast33_matrix_norm1(A_pow(:,:,4))**(1.0d0/10.0d0)
                eta4 = max(d8, d10)
                eta5 = min(eta3, eta4)
                s = max(0, ceiling(log((eta5 + 1.0d-300)/theta(5))/log(2.0d0)))
                s = s + fast33_ell(2.0d0**(-s)*A, 13)
                A = A*2.0d0**(-s)
                A_pow(:,:,1) = A_pow(:,:,1)*2.0d0**(-2*s)
                A_pow(:,:,2) = A_pow(:,:,2)*2.0d0**(-4*s)
                A_pow(:,:,3) = A_pow(:,:,3)*2.0d0**(-6*s)
                CALL fast33_pade13(U, V, ID, A, A_pow)
            END IF
        END IF
        Q  = V - U
        EA = V + U ! EA ( = R ) = P directly, since we ovwewrite the rhs in lu_back_substitution
        CALL fast33_ludcmp(Q, ludcmp_indx, ludcmp_d)
        DO i = 1, N
            CALL fast33_lubksb(Q, ludcmp_indx, EA(:,i))
        END DO
        DO i = 1, s
            EA = matmul(EA, EA)
        END DO
    END SUBROUTINE fast33_compute_matrix_exponential

    PURE FUNCTION fast33_ell(A, m)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        INTEGER,                 INTENT(IN) :: m
        INTEGER                             :: fast33_ell
        REAL(8), PARAMETER                  :: eps = 2.0d0**(-53)
        INTEGER                             :: p
        REAL(8)                             :: aci
        p = 2*m + 1
        aci = abs_c_inv(p)
        fast33_ell = max(0, ceiling(log((fast33_matrix_norm1_pow(abs(A), p) + 1.0d-300)/(&
            eps*fast33_matrix_norm1(A)*aci))/(log(2.0d0)*2.0d0*real(m, 8))))
    END FUNCTION fast33_ell

    PURE FUNCTION fast33_matrix_norm1(A)
        REAL(8), DIMENSION(3,3), INTENT(IN) :: A
        REAL(8)                             :: fast33_matrix_norm1
        ! take the largest l1 norm of all columns (sum over column)
        fast33_matrix_norm1 = maxval(sum(abs(A), dim=1)) 
    END FUNCTION fast33_matrix_norm1

    PURE SUBROUTINE fast33_pade3(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(3,3),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(3,3),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(3,3,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(4),     PARAMETER   :: b = [120.0d0, 60.0d0, 12.0d0, 1.0d0]
        U = matmul(A, b(4)*A_pow(:,:,1) + b(2)*ID)
        V = b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE fast33_pade3

    PURE SUBROUTINE fast33_pade5(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(3,3),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(3,3),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(3,3,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(6),     PARAMETER   :: b = [30240.0d0, 15120.0d0, 3360.0d0, &
            420.0d0, 30.0d0, 1.0d0]
        U = matmul(A, b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID)
        V = b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE fast33_pade5

    PURE SUBROUTINE fast33_pade7(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(3,3),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(3,3),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(3,3,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(8),     PARAMETER   :: b = [17297280.0d0, 8648640.0d0, 1995840.0d0, &
            277200.0d0, 25200.0d0, 1512.0d0, 56.0d0, 1.0d0]
        U = matmul(A, b(8)*A_pow(:,:,3) + b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID)
        V = b(7)*A_pow(:,:,3) + b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE fast33_pade7

    PURE SUBROUTINE fast33_pade9(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(3,3),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(3,3),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(3,3,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(10),    PARAMETER   :: b = real([17643225600_8, 8821612800_8, &
            2075673600_8, 302702400_8, 30270240_8, 2162160_8, 110880_8, 3960_8, 90_8, 1_8], 8)
        U = matmul(A, (b(10)*A_pow(:,:,4) + b(8)*A_pow(:,:,3) + b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID))
        V = (b(9)*A_pow(:,:,4) + b(7)*A_pow(:,:,3) + b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID)
    END SUBROUTINE fast33_pade9

    PURE SUBROUTINE fast33_pade13(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(3,3),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(3,3),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(3,3,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(14),    PARAMETER   :: b = real([64764752532480000_8, &
            32382376266240000_8, 7771770303897600_8, 1187353796428800_8, 129060195264000_8, &
            10559470521600_8, 670442572800_8, 33522128640_8, 1323241920_8, 40840800_8, 960960_8, &
            16380_8, 182_8, 1_8], 8)
        U  = matmul(A, matmul(A_pow(:,:,3), b(14)*A_pow(:,:,3) + b(12)*A_pow(:,:,2) + b(10)*A_pow(:,:,1)) + &
            b(8)*A_pow(:,:,3) + b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID)
        V  = matmul(A_pow(:,:,3), b(13)*A_pow(:,:,3) + b(11)*A_pow(:,:,2) + b(9)*A_pow(:,:,1)) + &
            b(7)*A_pow(:,:,3) + b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE fast33_pade13

    PURE SUBROUTINE fast33_ludcmp(a, indx, d)
        REAL(8), DIMENSION(3,3), INTENT(INOUT) :: a
        INTEGER, DIMENSION(3),   INTENT(OUT)   :: indx
        INTEGER,                 INTENT(OUT)   :: d
        REAL(8), DIMENSION(3)                  :: vv
        REAL(8)                                :: temp
        REAL(8), PARAMETER                     :: TINY = 1.0d-20
        INTEGER, PARAMETER                     :: n = 3
        ! REAL(8), PARAMETER                     :: VERYTINY = 1.0d-200
        INTEGER                                :: i, j, imax, ii, jj
        d = 1
        vv = maxval(abs(a), dim=2)
        ! IF (any(abs(vv) < VERYTINY)) THEN
        !     WRITE(*,*) "singular matrix in ludcmp"
        ! END IF
        vv = 1.0d0/vv 
        DO j = 1, n
            imax = (j - 1) + maxloc(vv(j:n)*abs(a(j:n,j)), 1) 
            IF (j .ne. imax) THEN
                ! swap 
                DO i = 1, n
                    temp = a(imax,i)
                    a(imax,i) = a(j,i)
                    a(j,i) = temp
                END DO
                d        = -d
                vv(imax) = vv(j)
            END IF
            indx(j) = imax
            IF (abs(a(j,j)) .lt. TINY) THEN
                ! WRITE(*,*) "TINY IN LUDCMP"
                a(j,j) = TINY
            END IF
            a(j+1:n,j)     = a(j+1:n,j)/a(j,j)
            ! a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - fast33_outerprod(a(j+1:n,j), a(j,j+1:n))
            ! explicit outerprod
            do jj = j+1,n
                do ii = j+1,n
                    a(ii,jj) = a(ii,jj) - a(ii,j)*a(j,jj)
                end do
            end do
        END DO
    END SUBROUTINE fast33_ludcmp

    PURE SUBROUTINE fast33_lubksb(a, indx, b)
        REAL(8), DIMENSION(3,3), INTENT(IN)    :: a
        INTEGER, DIMENSION(3),   INTENT(IN)    :: indx
        REAL(8), DIMENSION(3),   INTENT(INOUT) :: b
        INTEGER                                :: i, ii, ll
        INTEGER, PARAMETER                     :: n = 3
        REAL(8)                                :: summ
        ii = 0 
        DO i = 1, n
            ll = indx(i)
            summ = b(ll)
            b(ll) = b(i)
            IF (ii .ne. 0) THEN
                summ = summ - dot_product(a(i,ii:i-1), b(ii:i-1))
            ELSE IF (abs(summ) .gt. 1.0d-200) THEN
                ii = i 
            END IF
            b(i) = summ
        END DO
        DO i = n, 1, -1 
            b(i) = (b(i) - dot_product(a(i,i+1:n), b(i+1:n)))/a(i,i)
        END DO
    END SUBROUTINE fast33_lubksb

    PURE SUBROUTINE fast33_solve(a, x, b)
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: a
        REAL(8), DIMENSION(3),   INTENT(IN)  :: b
        REAL(8), DIMENSION(3),   INTENT(OUT) :: x
        INTEGER, DIMENSION(3)                :: ludcmp_indx
        REAL(8), DIMENSION(3,3)              :: ludcmp_a
        INTEGER                              :: ludcmp_d
        ludcmp_a = a
        x = b
        CALL fast33_ludcmp(ludcmp_a, ludcmp_indx, ludcmp_d)
        CALL fast33_lubksb(ludcmp_a, ludcmp_indx, x)
    END SUBROUTINE fast33_solve

    ! ---- fast 3by3 size helper functions for exact solution of first order linear ODEs --------- !

    PURE SUBROUTINE fast33_build_jacobian_matrix_mask(mask, nv_reduced, Jm, eps)
        ! builds a mask for empty rows of the jacobian matrix
        INTEGER, DIMENSION(3),   INTENT(OUT) :: mask
        INTEGER,                 INTENT(OUT) :: nv_reduced
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: Jm
        REAL(8),                 INTENT(IN)  :: eps
        REAL(8)                              :: scale
        INTEGER                              :: i, k
        scale = eps*(1.0d0 + maxval(abs(Jm)))
        mask = 1
        k = 0
        DO i = 1, 3 ! == size(mask)
            IF (maxval(abs(Jm(i,:))) .lt. scale) THEN
                mask(i) = 0
            ELSE 
                k = k + 1
                mask(i) = k
            END IF
        END DO
        nv_reduced = k
    END SUBROUTINE fast33_build_jacobian_matrix_mask

    PURE SUBROUTINE fast33_build_reduced_jacobian_matrix(Jm_reduced, S0_reduced, V0_reduced, Jm, S0, V0, mask)
        REAL(8), DIMENSION(3,3), INTENT(OUT) :: Jm_reduced
        REAL(8), DIMENSION(3),   INTENT(OUT) :: S0_reduced, V0_reduced
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: Jm
        REAL(8), DIMENSION(3),   INTENT(IN)  :: S0, V0
        INTEGER, DIMENSION(3),   INTENT(IN)  :: mask
        INTEGER                              :: i, j
        DO i = 1, 3  ! == size(mask)
            IF (mask(i) .gt. 0) THEN
                DO j = 1, 3
                    IF (mask(j) .gt. 0) THEN
                        Jm_reduced(mask(i),mask(j)) = Jm(i,j) 
                    END IF
                END DO
                S0_reduced(mask(i)) = S0(i)
                V0_reduced(mask(i)) = V0(i)
            END IF
        END DO    
    END SUBROUTINE fast33_build_reduced_jacobian_matrix

    PURE SUBROUTINE fast33_compute_first_order_linear_ode_solution(Ve, Jm_reduced, V0_reduced, S0_reduced, nv2, t)
        REAL(8), DIMENSION(3),   INTENT(OUT) :: Ve
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: Jm_reduced
        REAL(8), DIMENSION(3),   INTENT(IN)  :: V0_reduced
        REAL(8), DIMENSION(3),   INTENT(IN)  :: S0_reduced
        INTEGER,                 INTENT(IN)  :: nv2
        REAL(8),                 INTENT(IN)  :: t
        REAL(8), DIMENSION(nv2,nv2)          :: Jmt
        REAL(8), DIMENSION(nv2,nv2)          :: EJmt
        REAL(8), DIMENSION(nv2)              :: Vf_reduced
        IF (nv2 .eq. 3) THEN
            Jmt = Jm_reduced
            CALL fast33_solve(Jmt, Vf_reduced, S0_reduced)
            Jmt = Jmt*t
            CALL fast33_compute_matrix_exponential(EJmt, Jmt)
            ! compute the matrix exponential and the solution at the new time level
            Ve = matmul(EJmt, V0_reduced + Vf_reduced) - Vf_reduced
        ELSE IF (nv2 .eq. 2) THEN
            Jmt = Jm_reduced(:nv2,:nv2)
            CALL fast22_solve(Jmt, Vf_reduced, S0_reduced(:nv2))
            Jmt = Jmt*t
            CALL fast22_compute_matrix_exponential(EJmt, Jmt)
            ! compute the matrix exponential and the solution at the new time level
            Ve(:nv2) = matmul(EJmt, V0_reduced(:nv2) + Vf_reduced) - Vf_reduced
        ELSE IF (nv2 .eq. 1) THEN
            Vf_reduced(1) = S0_reduced(1)/Jm_reduced(1,1)
            EJmt(1,1) = exp(Jm_reduced(1,1)*t)
            Ve(1) = EJmt(1,1)*(V0_reduced(1) + Vf_reduced(1)) - Vf_reduced(1)
        ELSE
            Ve = 0.0d0
        END IF
    END SUBROUTINE fast33_compute_first_order_linear_ode_solution

    PURE SUBROUTINE fast33_select_solution_from_mask(V, V0, S0, Ve, mask, t)
        REAL(8), DIMENSION(3), INTENT(OUT) :: V
        REAL(8), DIMENSION(3), INTENT(IN)  :: V0, S0, Ve
        INTEGER, DIMENSION(3), INTENT(IN)  :: mask
        REAL(8),               INTENT(IN)  :: t
        INTEGER                            :: i
        DO i = 1, size(V0)
            IF (mask(i) .eq. 0) THEN
                ! for the rows where the jacobian had to be masked, plug in the trivial solution
                V(i) = V0(i) + t*S0(i)
            ELSE
                V(i) = Ve(mask(i))
            END IF
        END DO
    END SUBROUTINE fast33_select_solution_from_mask

    ! ---- also the two by two case is needed for 3 by 3 system solution

    PURE SUBROUTINE fast22_ludcmp(a, indx, d)
        REAL(8), DIMENSION(2,2), INTENT(INOUT) :: a
        INTEGER, DIMENSION(2),   INTENT(OUT)   :: indx
        INTEGER,                 INTENT(OUT)   :: d
        REAL(8), DIMENSION(2)                  :: vv
        REAL(8)                                :: temp
        REAL(8), PARAMETER                     :: TINY = 1.0d-20
        INTEGER, PARAMETER                     :: n = 2
        ! REAL(8), PARAMETER                     :: VERYTINY = 1.0d-200
        INTEGER                                :: i, j, imax, jj, ii
        d = 1
        vv = maxval(abs(a), dim=2)
        ! IF (any(abs(vv) < VERYTINY)) THEN
        !     WRITE(*,*) "singular matrix in ludcmp"
        ! END IF
        vv = 1.0d0/vv 
        DO j = 1, n
            imax = (j - 1) + maxloc(vv(j:n)*abs(a(j:n,j)), 1) 
            IF (j .ne. imax) THEN
                ! swap 
                DO i = 1, n
                    temp = a(imax,i)
                    a(imax,i) = a(j,i)
                    a(j,i) = temp
                END DO
                d        = -d
                vv(imax) = vv(j)
            END IF
            indx(j) = imax
            IF (abs(a(j,j)) .lt. TINY) THEN
                ! WRITE(*,*) "TINY IN LUDCMP"
                a(j,j) = TINY
            END IF
            a(j+1:n,j)     = a(j+1:n,j)/a(j,j)
            ! a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - fast22_outerprod(a(j+1:n,j), a(j,j+1:n))
            ! explicit outerprod
            do jj = j+1,n
                do ii = j+1,n
                    a(ii,jj) = a(ii,jj) - a(ii,j)*a(j,jj)
                end do
            end do
        END DO
    END SUBROUTINE fast22_ludcmp

    PURE SUBROUTINE fast22_lubksb(a, indx, b)
        REAL(8), DIMENSION(2,2), INTENT(IN)    :: a
        INTEGER, DIMENSION(2),   INTENT(IN)    :: indx
        REAL(8), DIMENSION(2),   INTENT(INOUT) :: b
        INTEGER                                :: i, ii, ll
        INTEGER, PARAMETER                     :: n = 3
        REAL(8)                                :: summ
        ii = 0 
        DO i = 1, n
            ll = indx(i)
            summ = b(ll)
            b(ll) = b(i)
            IF (ii .ne. 0) THEN
                summ = summ - dot_product(a(i,ii:i-1), b(ii:i-1))
            ELSE IF (abs(summ) .gt. 1.0d-200) THEN
                ii = i 
            END IF
            b(i) = summ
        END DO
        DO i = n, 1, -1 
            b(i) = (b(i) - dot_product(a(i,i+1:n), b(i+1:n)))/a(i,i)
        END DO
    END SUBROUTINE fast22_lubksb

    PURE SUBROUTINE fast22_solve(a, x, b)
        REAL(8), DIMENSION(2,2), INTENT(IN)  :: a
        REAL(8), DIMENSION(2),   INTENT(IN)  :: b
        REAL(8), DIMENSION(2),   INTENT(OUT) :: x
        INTEGER, DIMENSION(2)                :: ludcmp_indx
        REAL(8), DIMENSION(2,2)              :: ludcmp_a
        INTEGER                              :: ludcmp_d
        ludcmp_a = a
        x = b
        CALL fast22_ludcmp(ludcmp_a, ludcmp_indx, ludcmp_d)
        CALL fast22_lubksb(ludcmp_a, ludcmp_indx, x)
    END SUBROUTINE fast22_solve

    PURE SUBROUTINE fast22_compute_matrix_exponential(EA, A0)
        ! as the above, but with hardcoded 3 by 3 matrix size
        REAL(8), DIMENSION(2,2), INTENT(OUT) :: EA
        REAL(8), DIMENSION(2,2), INTENT(IN)  :: A0
        REAL(8), DIMENSION(5),   PARAMETER   :: theta = [1.495585217958292d-2, &
            2.539398330063230d-1, 9.504178996162932d-1, 2.097847961257068d0, 4.25d0]
        INTEGER :: i, s, ludcmp_d
        INTEGER, PARAMETER :: N = 2
        REAL(8), DIMENSION(2,2) :: ID, A, U, V, Q
        REAL(8), DIMENSION(2,2,4) :: A_pow
        INTEGER, DIMENSION(2) :: ludcmp_indx
        REAL(8) :: d4, d6, d8, d10, eta1, eta2, eta3, eta4, eta5
        s = 0
        ID = 0.0d0
        DO i = 1, N
            ID(i,i) = 1.0d0
        END DO
        A = A0
        A_pow(:,:,1) = matmul(A, A)
        A_pow(:,:,2) = matmul(A_pow(:,:,1), A_pow(:,:,1))
        A_pow(:,:,3) = matmul(A_pow(:,:,2), A_pow(:,:,1))
        d4 = fast22_matrix_norm1(A_pow(:,:,2))**(1.0d0/4.0d0)
        d6 = fast22_matrix_norm1(A_pow(:,:,3))**(1.0d0/6.0d0)
        eta1 = max(d4, d6)
        eta2 = eta1
        IF (eta1 .le. theta(1) .and. fast22_ell(A, 3) .eq. 0) THEN
            CALL fast22_pade3(U, V, ID, A, A_pow)
        ELSE IF (eta2 .le. theta(2) .and. fast22_ell(A, 5) .eq. 0) THEN
            CALL fast22_pade5(U, V, ID, A, A_pow)
        ELSE
            A_pow(:,:,4) = matmul(A_pow(:,:,3), A_pow(:,:,1))
            d8 = fast22_matrix_norm1(A_pow(:,:,4))**(1.0d0/8.0d0)
            eta3 = max(d6, d8)
            IF (eta3 .le. theta(3) .and. fast22_ell(A, 7) .eq. 0) THEN
                CALL fast22_pade7(U, V, ID, A, A_pow)
            ELSE IF (eta3 .le. theta(4) .and. fast22_ell(A, 9) .eq. 0) THEN
                CALL fast22_pade9(U, V, ID, A, A_pow)
            ELSE
                ! since the degree 13 pade approximant does not use A^8, we can overwrite A^10 on A^8
                A_pow(:,:,4) = matmul(A_pow(:,:,2), A_pow(:,:,3))
                d10 = fast22_matrix_norm1(A_pow(:,:,4))**(1.0d0/10.0d0)
                eta4 = max(d8, d10)
                eta5 = min(eta3, eta4)
                s = max(0, ceiling(log((eta5 + 1.0d-300)/theta(5))/log(2.0d0)))
                s = s + fast22_ell(2.0d0**(-s)*A, 13)
                A = A*2.0d0**(-s)
                A_pow(:,:,1) = A_pow(:,:,1)*2.0d0**(-2*s)
                A_pow(:,:,2) = A_pow(:,:,2)*2.0d0**(-4*s)
                A_pow(:,:,3) = A_pow(:,:,3)*2.0d0**(-6*s)
                CALL fast22_pade13(U, V, ID, A, A_pow)
            END IF
        END IF
        Q  = V - U
        EA = V + U ! EA ( = R ) = P directly, since we ovwewrite the rhs in lu_back_substitution
        CALL fast22_ludcmp(Q, ludcmp_indx, ludcmp_d)
        DO i = 1, N
            CALL fast22_lubksb(Q, ludcmp_indx, EA(:,i))
        END DO
        DO i = 1, s
            EA = matmul(EA, EA)
        END DO
    END SUBROUTINE fast22_compute_matrix_exponential

    PURE FUNCTION fast22_ell(A, m)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        INTEGER,                 INTENT(IN) :: m
        INTEGER                             :: fast22_ell
        REAL(8), PARAMETER                  :: eps = 2.0d0**(-53)
        INTEGER                             :: p
        REAL(8)                             :: aci
        p = 2*m + 1
        aci = abs_c_inv(p)
        fast22_ell = max(0, ceiling(log((fast22_matrix_norm1_pow(abs(A), p) + 1.0d-300)/(&
            eps*fast22_matrix_norm1(A)*aci))/(log(2.0d0)*2.0d0*real(m, 8))))
    END FUNCTION fast22_ell

    PURE FUNCTION fast22_matrix_norm1(A)
        REAL(8), DIMENSION(2,2), INTENT(IN) :: A
        REAL(8)                             :: fast22_matrix_norm1
        ! take the largest l1 norm of all columns (sum over column)
        fast22_matrix_norm1 = maxval(sum(abs(A), dim=1)) 
    END FUNCTION fast22_matrix_norm1

    PURE SUBROUTINE fast22_pade3(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(2,2),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(2,2),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(2,2,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(4),     PARAMETER   :: b = [120.0d0, 60.0d0, 12.0d0, 1.0d0]
        U = matmul(A, b(4)*A_pow(:,:,1) + b(2)*ID)
        V = b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE fast22_pade3

    PURE SUBROUTINE fast22_pade5(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(2,2),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(2,2),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(2,2,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(6),     PARAMETER   :: b = [30240.0d0, 15120.0d0, 3360.0d0, &
            420.0d0, 30.0d0, 1.0d0]
        U = matmul(A, b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID)
        V = b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE fast22_pade5

    PURE SUBROUTINE fast22_pade7(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(2,2),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(2,2),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(2,2,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(8),     PARAMETER   :: b = [17297280.0d0, 8648640.0d0, 1995840.0d0, &
            277200.0d0, 25200.0d0, 1512.0d0, 56.0d0, 1.0d0]
        U = matmul(A, b(8)*A_pow(:,:,3) + b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID)
        V = b(7)*A_pow(:,:,3) + b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE fast22_pade7

    PURE SUBROUTINE fast22_pade9(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(2,2),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(2,2),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(2,2,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(10),    PARAMETER   :: b = real([17643225600_8, 8821612800_8, &
            2075673600_8, 302702400_8, 30270240_8, 2162160_8, 110880_8, 3960_8, 90_8, 1_8], 8)
        U = matmul(A, (b(10)*A_pow(:,:,4) + b(8)*A_pow(:,:,3) + b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID))
        V = (b(9)*A_pow(:,:,4) + b(7)*A_pow(:,:,3) + b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID)
    END SUBROUTINE fast22_pade9

    PURE SUBROUTINE fast22_pade13(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(2,2),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(2,2),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(2,2,4), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(14),    PARAMETER   :: b = real([64764752532480000_8, &
            32382376266240000_8, 7771770303897600_8, 1187353796428800_8, 129060195264000_8, &
            10559470521600_8, 670442572800_8, 33522128640_8, 1323241920_8, 40840800_8, 960960_8, &
            16380_8, 182_8, 1_8], 8)
        U  = matmul(A, matmul(A_pow(:,:,3), b(14)*A_pow(:,:,3) + b(12)*A_pow(:,:,2) + b(10)*A_pow(:,:,1)) + &
            b(8)*A_pow(:,:,3) + b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID)
        V  = matmul(A_pow(:,:,3), b(13)*A_pow(:,:,3) + b(11)*A_pow(:,:,2) + b(9)*A_pow(:,:,1)) + &
            b(7)*A_pow(:,:,3) + b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE fast22_pade13

END MODULE expintegrator_linalg_fast33

#endif
#endif
