#ifndef HEADER_EXPINTEGRATOR_LINALG
#define HEADER_EXPINTEGRATOR_LINALG
!
#if defined(ODESOLVER) || defined(EQNTYPEC99) || defined(EQNTYPED99)
    
MODULE expintegrator_linalg

    IMPLICIT NONE

    PRIVATE
    
    ! gpr
    PUBLIC :: build_by_rows
    PUBLIC :: flatten_rows
    PUBLIC :: halven_rank42_rows
    PUBLIC :: diadic_product
    PUBLIC :: cross_product
    PUBLIC :: trace
    
    ! lu decomposition and matrix inversion
    PUBLIC  :: compute_matrix_inverse
    PRIVATE :: ludcmp_preallocated
    PUBLIC  :: solve_large
    PUBLIC  :: solve
    PUBLIC  :: compute_matrix_inverse_large
    PRIVATE :: ludcmp
    PRIVATE :: lubksb
    PRIVATE :: outerprod
    PUBLIC  :: determinant_lu
    
    ! matrix exponential
    PRIVATE :: matrix_norm1
    PRIVATE :: pade3
    PRIVATE :: pade5
    PRIVATE :: pade7
    PRIVATE :: pade9
    PRIVATE :: pade13
    PRIVATE :: matrix_norm1_pow
    PRIVATE :: abs_c_inv
    PRIVATE :: ell
    PUBLIC  :: compute_matrix_exponential
    PUBLIC  :: compute_matrix_exponential_simple

    ! ode
    PUBLIC  :: build_jacobian_matrix_mask
    PUBLIC  :: build_reduced_jacobian_matrix
    PUBLIC  :: compute_first_order_linear_ode_solution
    PUBLIC  :: select_solution_from_mask

CONTAINS

    ! subroutines for computation of source of gpr with rupture

    PURE FUNCTION trace(A)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8)                             :: trace
        INTEGER                             :: i
        trace = 0.0d0
        DO i = 1, size(A, 1)
            trace = trace + A(i,i)
        END DO
    END FUNCTION trace

    PURE FUNCTION diadic_product(a, b)
        REAL(8), DIMENSION(:), INTENT(IN)    :: a, b
        REAL(8), DIMENSION(size(a), size(b)) :: diadic_product
        diadic_product = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a))
    END FUNCTION diadic_product

    PURE FUNCTION cross_product(u, v)
        REAL(8), DIMENSION(3), INTENT(IN) :: u, v
        REAL(8), DIMENSION(3)             :: cross_product
        cross_product(1) = u(2)*v(3) - u(3)*v(2)
        cross_product(2) = u(3)*v(1) - u(1)*v(3)
        cross_product(3) = u(1)*v(2) - u(2)*v(1)
    END FUNCTION cross_product

    PURE FUNCTION build_by_rows(v)
        REAL(8), DIMENSION(9), INTENT(IN) :: v
        REAL(8), DIMENSION(3,3) :: build_by_rows
        build_by_rows(1,:) = v(1:3)
        build_by_rows(2,:) = v(4:6)
        build_by_rows(3,:) = v(7:9)
    END FUNCTION build_by_rows

    PURE FUNCTION flatten_rows(M)
        REAL(8), DIMENSION(3,3), INTENT(IN) :: M
        REAL(8), DIMENSION(9)               :: flatten_rows
        flatten_rows(1:3) = M(1,:)
        flatten_rows(4:6) = M(2,:)
        flatten_rows(7:9) = M(3,:)
    END FUNCTION flatten_rows

    PURE FUNCTION halven_rank42_rows(M)
        REAL(8), DIMENSION(3,3,3,3), INTENT(IN) :: M
        REAL(8), DIMENSION(9,9)                 :: halven_rank42_rows
        halven_rank42_rows(1:3,1:3) = M(:,:,1,1)
        halven_rank42_rows(1:3,4:6) = M(:,:,1,2)
        halven_rank42_rows(1:3,7:9) = M(:,:,1,3)
        halven_rank42_rows(4:6,1:3) = M(:,:,2,1)
        halven_rank42_rows(4:6,4:6) = M(:,:,2,2)
        halven_rank42_rows(4:6,7:9) = M(:,:,2,3)
        halven_rank42_rows(7:9,1:3) = M(:,:,3,1)
        halven_rank42_rows(7:9,4:6) = M(:,:,3,2)
        halven_rank42_rows(7:9,7:9) = M(:,:,3,3)
    END FUNCTION halven_rank42_rows

    ! LU

    PURE SUBROUTINE compute_matrix_inverse(ainv, a)
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: a
        REAL(8), DIMENSION(:,:), INTENT(OUT) :: ainv
        INTEGER, DIMENSION(size(a,1))           :: ludcmp_indx
        REAL(8), DIMENSION(size(a,1))           :: ludcmp_v
        REAL(8), DIMENSION(size(a,1),size(a,1)) :: LU
        INTEGER                              :: ludcmp_d
        INTEGER                              :: i, j, k, imax, ii, n
        REAL(8)                              :: temp, summ
        REAL(8), PARAMETER                   :: ludcmp_epsilon = 1.0d-20
        n = size(a,1)
        LU = a
        ! perform LU decomposition
        ludcmp_d = 1
        ludcmp_v = maxval(abs(LU), dim=2)
        ludcmp_v = 1.0d0/ludcmp_v 
        DO j = 1, n
            imax = (j - 1) + maxloc(ludcmp_v(j:n)*abs(LU(j:n,j)), 1) 
            IF (j .ne. imax) THEN
                ! swap 
                DO i = 1, n
                    temp = LU(imax,i)
                    LU(imax,i) = LU(j,i)
                    LU(j,i) = temp
                END DO
                ludcmp_d = -ludcmp_d
                ludcmp_v(imax) = ludcmp_v(j)
            END IF
            ludcmp_indx(j) = imax
            IF (abs(LU(j,j)) .lt. ludcmp_epsilon) THEN
                LU(j,j) = ludcmp_epsilon
            END IF
            LU(j+1:n,j)     = LU(j+1:n,j)/LU(j,j)
            LU(j+1:n,j+1:n) = LU(j+1:n,j+1:n) - outerprod(LU(j+1:n,j), LU(j,j+1:n))
        END DO
        ! perform LU decomposition - end
        ! LU = a
        ! CALL ludcmp_preallocated(LU, ludcmp_indx, ludcmp_d, ludcmp_v) !!!!!! for debugging
        ainv = 0.0d0
        DO j = 1, n
            ! perform LU back substitution
            ainv(j,j) = 1.0d0
            ii = 0 
            DO i = 1, n
                k = ludcmp_indx(i)
                summ = ainv(k,j)
                ainv(k,j) = ainv(i,j)
                IF (ii .ne. 0) THEN
                    summ = summ - dot_product(LU(i,ii:i-1), ainv(ii:i-1,j))
                ELSE IF (abs(summ) .gt. 1.0d-200) THEN
                    ii = i 
                END IF
                ainv(i,j) = summ
            END DO
            DO i = n, 1, -1 
                ainv(i,j) = (ainv(i,j) - dot_product(LU(i,i+1:n), &
                    ainv(i+1:n,j)))/LU(i,i)
            END DO
            ! perform LU back substitution - end
        END DO
    END SUBROUTINE compute_matrix_inverse

    PURE SUBROUTINE ludcmp_preallocated(a, indx, d, vv)
        REAL(8), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER, DIMENSION(:),   INTENT(OUT)   :: indx
        INTEGER,                 INTENT(OUT)   :: d
        REAL(8), DIMENSION(:),   INTENT(OUT)   :: vv
        REAL(8)                                :: temp
        REAL(8), PARAMETER                     :: TINY = 1.0d-20
        INTEGER                                :: i, j, n, imax
        n = size(indx)
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
            a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - outerprod(a(j+1:n,j), a(j,j+1:n))
        END DO
    END SUBROUTINE ludcmp_preallocated 

    PURE SUBROUTINE solve_large(a, x, b)
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: a
        REAL(8), DIMENSION(:),   INTENT(IN)  :: b
        REAL(8), DIMENSION(:),   INTENT(OUT) :: x
        INTEGER, DIMENSION(:),   ALLOCATABLE :: ludcmp_indx
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: LU
        REAL(8), DIMENSION(:),   ALLOCATABLE :: ludcmp_v
        INTEGER                              :: ludcmp_d
        ALLOCATE(LU(size(b), size(b)), ludcmp_indx(size(b)), ludcmp_v(size(b)))
        LU = a
        x = b
        CALL ludcmp_preallocated(LU, ludcmp_indx, ludcmp_d, ludcmp_v)
        CALL lubksb(LU, ludcmp_indx, x)
    END SUBROUTINE solve_large

    PURE SUBROUTINE solve(a, x, b)
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: a
        REAL(8), DIMENSION(:),   INTENT(IN)  :: b
        REAL(8), DIMENSION(:),   INTENT(OUT) :: x
        INTEGER, DIMENSION(size(b))          :: ludcmp_indx
        REAL(8), DIMENSION(size(b),size(b))  :: ludcmp_a
        INTEGER                              :: ludcmp_d
        ludcmp_a = a
        x = b
        CALL ludcmp(ludcmp_a, ludcmp_indx, ludcmp_d)
        CALL lubksb(ludcmp_a, ludcmp_indx, x)
    END SUBROUTINE solve

    PURE SUBROUTINE compute_matrix_inverse_large(ainv, a)
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: a
        REAL(8), DIMENSION(:,:), INTENT(OUT) :: ainv
        INTEGER, DIMENSION(:),   ALLOCATABLE :: ludcmp_indx
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: LU
        INTEGER                              :: ludcmp_d
        INTEGER                              :: i
        ALLOCATE(LU(size(a,1), size(a,1)), ludcmp_indx(size(a,1)))
        LU = a
        CALL ludcmp(LU, ludcmp_indx, ludcmp_d)
        ainv = 0.0d0
        DO i = 1, size(a,1)
            ainv(i,i) = 1.0d0
            CALL lubksb(LU, ludcmp_indx, ainv(:,i))
        END DO
        DEALLOCATE(LU, ludcmp_indx)
    END SUBROUTINE compute_matrix_inverse_large

    PURE SUBROUTINE ludcmp(a, indx, d)
        REAL(8), DIMENSION(:,:), INTENT(INOUT) :: a
        INTEGER, DIMENSION(:),   INTENT(OUT)   :: indx
        INTEGER,                 INTENT(OUT)   :: d
        REAL(8), DIMENSION(size(A,1))          :: vv
        REAL(8)                                :: temp
        REAL(8), PARAMETER                     :: TINY = 1.0d-20
        ! REAL(8), PARAMETER                     :: VERYTINY = 1.0d-200
        INTEGER                                :: i, j, n, imax
        n = size(indx)
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
            a(j+1:n,j+1:n) = a(j+1:n,j+1:n) - outerprod(a(j+1:n,j), a(j,j+1:n))
        END DO
    END SUBROUTINE ludcmp  

    PURE SUBROUTINE lubksb(a, indx, b)
        REAL(8), DIMENSION(:,:), INTENT(IN)    :: a
        INTEGER, DIMENSION(:),   INTENT(IN)    :: indx
        REAL(8), DIMENSION(:),   INTENT(INOUT) :: b
        INTEGER                                :: i, n, ii, ll
        REAL(8)                                :: summ
        n = size(indx)
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
    END SUBROUTINE lubksb

    PURE FUNCTION outerprod(a, b)
        REAL(8), DIMENSION(:), INTENT(IN)   :: a, b
        REAL(8), DIMENSION(size(a),size(b)) :: outerprod
        outerprod = spread(a, dim=2, ncopies=size(b))*spread(b, dim=1, ncopies=size(a))
    END FUNCTION outerprod

    PURE FUNCTION determinant_lu(A)
        REAL(8), DIMENSION(:,:), INTENT(IN)     :: A
        REAL(8), DIMENSION(size(A,1),size(A,1)) :: Alu
        INTEGER, DIMENSION(size(A,1))           :: indx
        REAL(8)                                 :: determinant_lu
        INTEGER                                 :: dd, i
        Alu = A
        CALL ludcmp(Alu, indx, dd)
        determinant_lu = real(dd, 8)
        DO i = 1, size(A, 1)
            determinant_lu = determinant_lu*Alu(i,i)
        END DO
    END FUNCTION determinant_lu

! robust matrix exponential

    PURE SUBROUTINE pade3(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(:,:),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(:,:),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(:,:,:), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(4),     PARAMETER   :: b = [120.0d0, 60.0d0, 12.0d0, 1.0d0]
        U = matmul(A, b(4)*A_pow(:,:,1) + b(2)*ID)
        V = b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE pade3

    PURE SUBROUTINE pade5(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(:,:),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(:,:),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(:,:,:), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(6),     PARAMETER   :: b = [30240.0d0, 15120.0d0, 3360.0d0, &
            420.0d0, 30.0d0, 1.0d0]
        U = matmul(A, b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID)
        V = b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE pade5

    PURE SUBROUTINE pade7(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(:,:),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(:,:),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(:,:,:), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(8),     PARAMETER   :: b = [17297280.0d0, 8648640.0d0, 1995840.0d0, &
            277200.0d0, 25200.0d0, 1512.0d0, 56.0d0, 1.0d0]
        U = matmul(A, b(8)*A_pow(:,:,3) + b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID)
        V = b(7)*A_pow(:,:,3) + b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE pade7

    PURE SUBROUTINE pade9(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(:,:),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(:,:),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(:,:,:), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(10),    PARAMETER   :: b = real([17643225600_8, 8821612800_8, &
            2075673600_8, 302702400_8, 30270240_8, 2162160_8, 110880_8, 3960_8, 90_8, 1_8], 8)
        U = matmul(A, (b(10)*A_pow(:,:,4) + b(8)*A_pow(:,:,3) + b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID))
        V = (b(9)*A_pow(:,:,4) + b(7)*A_pow(:,:,3) + b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID)
    END SUBROUTINE pade9

    PURE SUBROUTINE pade13(U, V, ID, A, A_pow)
        REAL(8), DIMENSION(:,:),   INTENT(OUT) :: U, V
        REAL(8), DIMENSION(:,:),   INTENT(IN)  :: ID, A
        REAL(8), DIMENSION(:,:,:), INTENT(IN)  :: A_pow
        REAL(8), DIMENSION(14),    PARAMETER   :: b = real([64764752532480000_8, &
            32382376266240000_8, 7771770303897600_8, 1187353796428800_8, 129060195264000_8, &
            10559470521600_8, 670442572800_8, 33522128640_8, 1323241920_8, 40840800_8, 960960_8, &
            16380_8, 182_8, 1_8], 8)
        U  = matmul(A, matmul(A_pow(:,:,3), b(14)*A_pow(:,:,3) + b(12)*A_pow(:,:,2) + b(10)*A_pow(:,:,1)) + &
            b(8)*A_pow(:,:,3) + b(6)*A_pow(:,:,2) + b(4)*A_pow(:,:,1) + b(2)*ID)
        V  = matmul(A_pow(:,:,3), b(13)*A_pow(:,:,3) + b(11)*A_pow(:,:,2) + b(9)*A_pow(:,:,1)) + &
            b(7)*A_pow(:,:,3) + b(5)*A_pow(:,:,2) + b(3)*A_pow(:,:,1) + b(1)*ID
    END SUBROUTINE pade13

    PURE FUNCTION matrix_norm1(A)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        REAL(8)                             :: matrix_norm1
        ! take the largest l1 norm of all columns (sum over column)
        matrix_norm1 = maxval(sum(abs(A), dim=1)) 
    END FUNCTION matrix_norm1

    PURE FUNCTION matrix_norm1_pow(A, p)
        REAL(8), DIMENSION(:,:), INTENT(IN)      :: A
        INTEGER,                 INTENT(IN)      :: p
        REAL(8), DIMENSION(size(A,1))            :: v
        REAL(8), DIMENSION(size(A,1), size(A,2)) :: AT
        REAL(8)                                  :: matrix_norm1_pow
        INTEGER                                  :: i
        v = 1.0d0
        AT = transpose(A)
        DO i = 1, p
            v = matmul(AT, v)
        END DO
        matrix_norm1_pow = maxval(v)
    END FUNCTION matrix_norm1_pow

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

    PURE FUNCTION ell(A, m)
        REAL(8), DIMENSION(:,:), INTENT(IN) :: A
        INTEGER,                 INTENT(IN) :: m
        INTEGER                             :: ell
        REAL(8), PARAMETER                  :: eps = 2.0d0**(-53)
        INTEGER                             :: p
        REAL(8)                             :: aci
        p = 2*m + 1
        aci = abs_c_inv(p)
        ell = max(0, ceiling(log((matrix_norm1_pow(abs(A), p) + 1.0d-300)/(&
            eps*matrix_norm1(A)*aci))/(log(2.0d0)*2.0d0*real(m, 8))))
    END FUNCTION ell

    PURE SUBROUTINE compute_matrix_exponential(EA, A0)
        REAL(8), DIMENSION(:,:), INTENT(OUT) :: EA
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: A0
        REAL(8), DIMENSION(5),   PARAMETER   :: theta = [1.495585217958292d-2, &
            2.539398330063230d-1, 9.504178996162932d-1, 2.097847961257068d0, 4.25d0]
        INTEGER :: N, i, s, ludcmp_d
        REAL(8), DIMENSION(size(A0,1),size(A0,1)) :: ID, A, U, V, Q
        REAL(8), DIMENSION(size(A0,1),size(A0,1),4) :: A_pow
        INTEGER, DIMENSION(size(A0,1)) :: ludcmp_indx
        REAL(8) :: d4, d6, d8, d10, eta1, eta2, eta3, eta4, eta5
        N = size(A0,1)
        s = 0
        ID = 0.0d0
        DO i = 1, N
            ID(i,i) = 1.0d0
        END DO
        A = A0
        A_pow(:,:,1) = matmul(A, A)
        A_pow(:,:,2) = matmul(A_pow(:,:,1), A_pow(:,:,1))
        A_pow(:,:,3) = matmul(A_pow(:,:,2), A_pow(:,:,1))
        d4 = matrix_norm1(A_pow(:,:,2))**(1.0d0/4.0d0)
        d6 = matrix_norm1(A_pow(:,:,3))**(1.0d0/6.0d0)
        eta1 = max(d4, d6)
        eta2 = eta1
        IF (eta1 .le. theta(1) .and. ell(A, 3) .eq. 0) THEN
            CALL pade3(U, V, ID, A, A_pow)
        ELSE IF (eta2 .le. theta(2) .and. ell(A, 5) .eq. 0) THEN
            CALL pade5(U, V, ID, A, A_pow)
        ELSE
            A_pow(:,:,4) = matmul(A_pow(:,:,3), A_pow(:,:,1))
            d8 = matrix_norm1(A_pow(:,:,4))**(1.0d0/8.0d0)
            eta3 = max(d6, d8)
            IF (eta3 .le. theta(3) .and. ell(A, 7) .eq. 0) THEN
                CALL pade7(U, V, ID, A, A_pow)
            ELSE IF (eta3 .le. theta(4) .and. ell(A, 9) .eq. 0) THEN
                CALL pade9(U, V, ID, A, A_pow)
            ELSE
                ! since the degree 13 pade approximant does not use A^8, we can overwrite A^10 on A^8
                A_pow(:,:,4) = matmul(A_pow(:,:,2), A_pow(:,:,3))
                d10 = matrix_norm1(A_pow(:,:,4))**(1.0d0/10.0d0)
                eta4 = max(d8, d10)
                eta5 = min(eta3, eta4)
                s = max(0, ceiling(log((eta5 + 1.0d-300)/theta(5))/log(2.0d0)))
                s = s + ell(2.0d0**(-s)*A, 13)
                A = A*2.0d0**(-s)
                A_pow(:,:,1) = A_pow(:,:,1)*2.0d0**(-2*s)
                A_pow(:,:,2) = A_pow(:,:,2)*2.0d0**(-4*s)
                A_pow(:,:,3) = A_pow(:,:,3)*2.0d0**(-6*s)
                CALL pade13(U, V, ID, A, A_pow)
            END IF
        END IF
        Q  = V - U
        EA = V + U ! EA ( = R ) = P directly, since we ovwewrite the rhs in lu_back_substitution
        CALL ludcmp(Q, ludcmp_indx, ludcmp_d)
        DO i = 1, N
            CALL lubksb(Q, ludcmp_indx, EA(:,i))
        END DO
        DO i = 1, s
            EA = matmul(EA, EA)
        END DO
    END SUBROUTINE compute_matrix_exponential

    PURE SUBROUTINE compute_matrix_exponential_simple(EA, A0)
        REAL(8), DIMENSION(:,:), INTENT(OUT) :: EA
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: A0
        REAL(8), DIMENSION(13),  PARAMETER   :: b = real([32382376266240000_8, 7771770303897600_8, &
            1187353796428800_8, 129060195264000_8, 10559470521600_8, 670442572800_8, 33522128640_8, &
            1323241920_8, 40840800_8, 960960_8, 16380_8, 182_8, 1_8], 8)
        REAL(8), DIMENSION(5),   PARAMETER  :: theta = [1.495585217958292d-2, 2.539398330063230d-1, &
            9.504178996162932d-1, 2.097847961257068d0, 5.371920351148152d0]
        REAL(8), PARAMETER :: b0 = real(64764752532480000_8, 8)
        INTEGER :: N, i, imax, s, ludcmp_d
        REAL(8) :: mu, Anorm
        REAL(8), DIMENSION(size(A0,1),size(A0,1)) :: ID, A, U, V, R, Q, Uk ! maybe can eliminate Uk
        REAL(8), DIMENSION(size(A0,1),size(A0,1),4) :: A_pow
        INTEGER, DIMENSION(size(A0,1)) :: ludcmp_indx
        N = size(A0,1)
        ! mu = trace(A0)/real(N, 8)
        mu = 0.0d0 ! dont extract the trace for now
        ID = 0.0d0
        DO i = 1, N
            ID(i,i) = 1.0d0
        END DO
        A = A0 - mu*ID
        ! Anorm = matrix_norm1_estimate_simple(A)
        Anorm = matrix_norm1(A)
        imax = 5 ! for 13 and above 13
        DO i = 1, 5
            IF (Anorm .le. theta(i)) THEN
                imax = i
                EXIT
            END IF
        END DO
        s = ceiling(log((Anorm + 1.0d-300)/theta(5))/log(2.0d0)) ! note: theta(5) = theta13
        ! evaluate pade approximant
        IF (s .lt. 1) THEN
            A_pow(:,:,1) = matmul(A, A)
            U = b(1)*ID + b(3)*A_pow(:,:,1)
            V = b0*ID   + b(2)*A_pow(:,:,1)
            DO i = 2, min(imax, 3)
                A_pow(:,:,i) = matmul(A_pow(:,:,1), A_pow(:,:,i-1))
                U = U + b(2*i+1)*A_pow(:,:,i)
                V = V + b(2*i  )*A_pow(:,:,i)
            END DO
            IF (imax .eq. 4) THEN ! m .eq. 9 (m = 2*imax + 1)
                A_pow(:,:,4) = matmul(A_pow(:,:,1), A_pow(:,:,i-1))
                U = U + b(9)*A_pow(:,:,4)
                V = V + b(8)*A_pow(:,:,4)
            END IF
            IF (imax .le. 4) THEN ! m .le. 9 (m = 2*imax + 1)
                Uk = matmul(A, U)
            ELSE         
                Uk = matmul(A, U + matmul(A_pow(:,:,3), b(9)*A_pow(:,:,1) + b(11)*A_pow(:,:,2) + b(13)*A_pow(:,:,3)))
                V = V + matmul(A_pow(:,:,3), b(8)*A_pow(:,:,1) + b(10)*A_pow(:,:,2) + b(12)*A_pow(:,:,3))
            END IF
        ELSE
            A = A/real(2**s, 8)
            A_pow(:,:,1) = matmul(A, A)
            U = b(1)*ID + b(3)*A_pow(:,:,1)
            V = b0*ID   + b(2)*A_pow(:,:,1)
            DO i = 2, 3
                A_pow(:,:,i) = matmul(A_pow(:,:,1), A_pow(:,:,i-1))
                U = U + b(2*i+1)*A_pow(:,:,i)
                V = V + b(2*i  )*A_pow(:,:,i)
            END DO
            Uk = matmul(A, U + matmul(A_pow(:,:,3), b(9)*A_pow(:,:,1) + b(11)*A_pow(:,:,2) + b(13)*A_pow(:,:,3)))
            V = V + matmul(A_pow(:,:,3), b(8)*A_pow(:,:,1) + b(10)*A_pow(:,:,2) + b(12)*A_pow(:,:,3))
        END IF
        Q = V - Uk
        R = V + Uk ! R = P directly, since we ovwewrite the rhs in lu_back_substitution
        CALL ludcmp(Q, ludcmp_indx, ludcmp_d)
        DO i = 1, N
            CALL lubksb(Q, ludcmp_indx, R(:,i))
        END DO
        IF (s .ge. 1) THEN
            R = matmul(R, R)
            DO i = 2, s
                R = matmul(R, R)
            END DO
        END IF
        EA = exp(mu)*R
    END SUBROUTINE compute_matrix_exponential_simple

    ! ---- helper functions for exact solution of first order linear ODEs ------------------------ !

    PURE SUBROUTINE build_jacobian_matrix_mask(mask, nv_reduced, Jm, eps)
        ! builds a mask for empty rows of the jacobian matrix
        INTEGER, DIMENSION(:),   INTENT(OUT) :: mask
        INTEGER,                 INTENT(OUT) :: nv_reduced
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: Jm
        REAL(8),                 INTENT(IN)  :: eps
        REAL(8)                              :: scale
        INTEGER                              :: i, k
        scale = eps*(1.0d0 + maxval(abs(Jm)))
        mask = 1
        k = 0
        DO i = 1, size(Jm, 1) ! == size(mask)
            IF (maxval(abs(Jm(i,:))) .lt. scale) THEN
                mask(i) = 0
            ELSE 
                k = k + 1
                mask(i) = k
            END IF
        END DO
        nv_reduced = k
    END SUBROUTINE build_jacobian_matrix_mask

    PURE SUBROUTINE build_reduced_jacobian_matrix(Jm_reduced, S0_reduced, V0_reduced, Jm, S0, V0, mask)
        REAL(8), DIMENSION(:,:), INTENT(OUT) :: Jm_reduced
        REAL(8), DIMENSION(:),   INTENT(OUT) :: S0_reduced, V0_reduced
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: Jm
        REAL(8), DIMENSION(:),   INTENT(IN)  :: S0, V0
        INTEGER, DIMENSION(:),   INTENT(IN)  :: mask
        INTEGER                              :: i, j
        DO i = 1, size(Jm, 1)  ! == size(mask)
            IF (mask(i) .gt. 0) THEN
                DO j = 1, size(Jm, 1)
                    IF (mask(j) .gt. 0) THEN
                        Jm_reduced(mask(i),mask(j)) = Jm(i,j) 
                    END IF
                END DO
                S0_reduced(mask(i)) = S0(i)
                V0_reduced(mask(i)) = V0(i)
            END IF
        END DO    
    END SUBROUTINE build_reduced_jacobian_matrix

    PURE SUBROUTINE compute_first_order_linear_ode_solution(Ve, Jm_reduced, V0_reduced, S0_reduced, nv2, t)
        REAL(8), DIMENSION(:),   INTENT(OUT) :: Ve
        REAL(8), DIMENSION(:,:), INTENT(IN)  :: Jm_reduced
        REAL(8), DIMENSION(:),   INTENT(IN)  :: V0_reduced
        REAL(8), DIMENSION(:),   INTENT(IN)  :: S0_reduced
        INTEGER,                 INTENT(IN)  :: nv2
        REAL(8),                 INTENT(IN)  :: t
        REAL(8), DIMENSION(nv2,nv2)          :: Jmt
        REAL(8), DIMENSION(nv2,nv2)          :: EJmt
        REAL(8), DIMENSION(nv2)              :: Vf_reduced
        IF (nv2 .gt. 0) THEN
            Jmt = Jm_reduced(:nv2,:nv2)
            CALL solve(Jmt, Vf_reduced, S0_reduced(:nv2))
            Jmt = Jmt*t
            CALL compute_matrix_exponential(EJmt, Jmt)
            ! compute the matrix exponential and the solution at the new time level
            Ve(:nv2) = matmul(EJmt, V0_reduced(:nv2) + Vf_reduced) - Vf_reduced
        ELSE
            Ve = 0.0d0
        END IF
    END SUBROUTINE compute_first_order_linear_ode_solution

    PURE SUBROUTINE select_solution_from_mask(V, V0, S0, Ve, mask, t)
        REAL(8), DIMENSION(:), INTENT(OUT) :: V
        REAL(8), DIMENSION(:), INTENT(IN)  :: V0, S0, Ve
        INTEGER, DIMENSION(:), INTENT(IN)  :: mask
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
    END SUBROUTINE select_solution_from_mask

END MODULE expintegrator_linalg
#endif
    
#endif
