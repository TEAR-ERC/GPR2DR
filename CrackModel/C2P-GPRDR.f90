#ifndef HEADER_CPGPRDR
#define HEADER_CPGPRDR

! The Con2Prim and Prim2Con routines for MHD.
! Should be merged with SRHD's.

 ! #include "MainVariables.f90"
 ! #include "SpecificVarEqn99.f90"
 ! #include "expintegrator_ode.f90"
 ! #include "expintegrator_type.f90"


RECURSIVE SUBROUTINE PDEPrim2Cons(Q,V)
  USE MainVariables, ONLY: nVar, nDim, EQN
#if defined(EQNTYPEC99) || defined(EQNTYPED99)
    USE SpecificVarEqn99, only : compute_l_m_mix
#endif
#if defined(EQNTYPED99)
    USE expintegrator_ode, only : convert_parameters
    use expintegrator_type, only : t_ode_parameters
#endif
  IMPLICIT NONE
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: V
  INTENT(OUT) :: Q 
  ! Local variable declaration
  REAL :: LL_gpr, MM_gpr, dMM_gpr, dKK_gpr, Kelast, rhoeh,rhoevv
  REAL :: A(3,3),detA, S, G(3,3), trG, devG(3,3), Id(3,3),temp(3,3)
#ifdef EQNTYPED99
    type(t_ode_parameters) :: ODE
#endif
  !
  Q=V
    A(1,:) = (/ V( 9), V(10), V(11) /) 
    A(2,:) = (/ V(12), V(13), V(14) /) 
    A(3,:) = (/ V(15), V(16), V(17) /) 
    detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
    G      = MATMUL( TRANSPOSE(A), A ) 
    trG    = G(1,1) + G(2,2) + G(3,3) 
    Id     = 0. 
    Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0     
    devG   = G - 1./3.*trG*Id
    temp   = detA*MATMUL( TRANSPOSE(devG), devG )    ! Changed
    
    S      = V(5)
    !
    !call ComputeLMc(LL_gpr,MM_gpr,dMM_gpr,dKK_gpr,V(20),V(23),V(19),V(24),V(21))
#ifdef EQNTYPED99
    call convert_parameters(V, ODE, EQN)
    call compute_l_m_mix(LL_GPR, MM_GPR,V(19),V(20),ODE%lam1lam2ratio,ODE%mu1mu2ratio,V(21))
#else
    call compute_l_m_mix(LL_GPR, MM_GPR,V(19),V(20),V(23),V(22),V(21))
#endif
    !call compute_l_m_mix(LL_GPR, MM_GPR,V(19),V(20),V(23),V(22),V(21))
    ! Compute the pressure as p=rho^2*E_rho -------------------------------------------------------------------------------
    Kelast  = LL_gpr+2.0/3.0*MM_gpr
    select case(EQN%EOS)
        case(1)
            rhoeh   = Kelast/(EQN%gamma*(EQN%gamma-1))*(detA)**(EQN%gamma)*EXP(S/EQN%cv) + (Kelast-EQN%gamma*EQN%p0)/(EQN%gamma)  ! internal energy as a function of density and entropy     
        case(2)
            rhoeh   = 0.5*detA*Kelast*(1-detA)**2+V(1)*EQN%cv*detA**EQN%gamma*EQN%T0*(EXP(S/EQN%cv)-1.)
    end select
    rhoevv  = MM_gpr/4.0*(temp(1,1)+temp(2,2)+temp(3,3))
    
    !
    Q=V
    Q(1)=V(1)
    Q(2:4)=V(2:4)*V(18)*V(1)
    Q(25:27)=V(25:27)
    Q(5)   = rhoeh + 0.5*V(1)*( V(2)**2+V(3)**2+V(4)**2 ) + 0.5*V(1)*EQN%ch**2*( V(6)**2 + V(7)**2 + V(8)**2 ) + rhoevv
END SUBROUTINE PDEPrim2Cons

RECURSIVE SUBROUTINE PDECons2Prim(V,Q)
  USE MainVariables, ONLY: nVar, nDim, EQN
#if defined(EQNTYPEC99) || defined(EQNTYPED99)
    USE SpecificVarEqn99, only : compute_l_m_mix
#endif
#if defined(EQNTYPED99)
    USE expintegrator_ode, only : convert_parameters
    use expintegrator_type      , only : t_ode_parameters
#endif
  IMPLICIT NONE
  !--------------------------------------------!
  ! Argument list declaration
  REAL :: Q(nVar), V(nVar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: V 
  ! Local variables
  REAL :: U(3), dKK_gpr, dMM_gpr, ec, evvc,ffc,Kelast, KK_gpr, LL_gpr, MM_gpr,ejc
  REAL :: alpha,ialpha, A(3,3),detA, G(3,3), trG, Id(3,3), devG(3,3),tempp(3,3)
  REAL :: rhoE3,rhoE2, CIE,KIC
#ifdef EQNTYPED99
  type(t_ode_parameters) :: ODE
#endif
  V=Q
  
	alpha=Q(18)
    if(EQN%epsilon1<0) then
        ialpha=1.0/alpha
    else
        ialpha=alpha/(alpha**2+EQN%epsilon1*(1-alpha))
    end if
    
    A(1,:)  = (/ Q( 9), Q(10), Q(11) /) 
    A(2,:)  = (/ Q(12), Q(13), Q(14) /) 
    A(3,:)  = (/ Q(15), Q(16), Q(17) /) 
    detA    = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
    G       = MATMUL( TRANSPOSE(A), A ) 
    trG     = G(1,1) + G(2,2) + G(3,3) 
    Id      = 0. 
    Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0     
    devG    = G - 1./3.*trG*Id
    tempp    = detA*MATMUL( TRANSPOSE(devG), devG ) ! changed
#ifdef EQNTYPED99
    call convert_parameters(Q, ODE, EQN)
    call compute_l_m_mix(LL_GPR, MM_GPR,Q(19),Q(20),ODE%lam1lam2ratio,ODE%mu1mu2ratio,Q(21))
#else
    call compute_l_m_mix(LL_GPR, MM_GPR,Q(19),Q(20),Q(23),Q(22),Q(21))
#endif
    !call compute_l_m_mix(LL_GPR, MM_GPR,Q(19),Q(20),Q(23),Q(22),Q(21))
    ! Compute the pressure as p=rho^2*E_rho -------------------------------------------------------------------------------
    Kelast  = LL_gpr+2.0/3.0*MM_gpr
    u       =Q(2:4)/Q(1)*ialpha
    V(25:27)=Q(25:27)
    select case(EQN%EOS)
        case(1)
            ec      = 1./detA/Kelast*(Q(5)-0.5*Q(1)*sum(u(1:3)**2))
            evvc    = 1./4.*MM_gpr/Kelast*(tempp(1,1)+tempp(2,2)+tempp(3,3))
            ffc     = 1./detA*(1./EQN%gamma+EQN%p0/Kelast)
            ejc     = 0.5*Q(1)/detA*EQN%ch**2*( Q(6)**2 + Q(7)**2 + Q(8)**2 )
    
            V       = Q
            V(2:4)  = u(1:3)
            V(5)    = EQN%cv*LOG(         (   ec-evvc-ejc-ffc   )*EQN%gamma*(EQN%gamma-1.0)*detA**(1-EQN%gamma)            )
        case(2)
            
            rhoE3   = 0.5*Q(1)*sum(u(1:3)**2)
            rhoE2   = 0.25*MM_gpr*(tempp(1,1)+tempp(2,2)+tempp(3,3))+0.5*Q(1)*EQN%ch**2*( Q(6)**2 + Q(7)**2 + Q(8)**2 )
            CIE     = 0.5*detA*Kelast*(1-detA)**2
            KIC     = Q(1)*EQN%cv*EQN%T0*detA**EQN%gamma
            ! Cons2prim assignaation
            V       = Q
            V(2:4)  = u(1:3)
            if(EQN%T0 .eq. 0) then
                V(5)    = 0.   
            else
                V(5)    = EQN%cv*LOG(1./KIC*(Q(5)-rhoE2-rhoE3+KIC-CIE))
            end if
    end select
END SUBROUTINE PDECons2Prim

#endif
