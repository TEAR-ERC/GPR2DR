#ifndef HEADER_PDE
#define HEADER_PDE

!#include "MainVariables.f90"
!#include "SpecificVarEqn99.f90"
!#include "expintegrator_ode.f90"
!#include "expintegrator_type.f90"
!#include "SpecificVarEqn99.f90"

! GRMHD PDE.f90
! Trento (EQNTYPE4)

RECURSIVE SUBROUTINE DynamicRupture(x_in, t, Q,slp)
  USE, INTRINSIC :: ISO_C_BINDING
  USE MainVariables, ONLY :  nVar, nDim, EQN
  USE SpecificVarEqn99, only : compute_total_stress,compute_hydrodynamic_pressure, compute_temperature, computeLMc, compute_l_m_mix, computeGPRLEstressGen
  USE expintegrator_ode , only : convert_parameters
  use expintegrator_type, only : t_ode_parameters

        IMPLICIT NONE
        type(t_ode_parameters) :: ODE

        ! Argument list 
        REAL, INTENT(IN)               :: x_in(nDim), t,slp   ! 
        REAL, INTENT(OUT)              :: Q(nVar)        ! 
        ! Local variables
	INTEGER :: iErr
        REAL :: stressnorm, sigma_vec(6),ee,u(3),xGP(3), V(nVar), LL_GPR, MM_GPR
	REAL :: mus,mud,dc, mu_t
    REAL :: syy,sxy,pf,TT(3,3),A(3,3), detA,Id(3,3),devG(3,3), GT(3,3), Yeq,Ydev,Yp
        ! Compute the normal stress using the Formula by Barton (IJNMF 2009)
        xGP=0.
        xGP(1:nDim)=x_in(1:nDim)

      !call computeGPRLEstressGen(aux(7),aux(1:6),Q(9:17), LL_gpr,MM_gpr, 0.,0.,.true.)
      CALL PDECons2Prim(V,Q,xGP,t,iErr)
      call convert_parameters(V, ODE, EQN)
      call compute_l_m_mix(LL_GPR, MM_GPR,V(19),V(20),ODE%lam1lam2ratio,ODE%mu1mu2ratio,V(21))
      A(1,:) = (/ Q( 9), Q(10), Q(11) /) 
      A(2,:) = (/ Q(12), Q(13), Q(14) /) 
      A(3,:) = (/ Q(15), Q(16), Q(17) /) 
      detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
      Id     = 0. 
      Id(1,1) = 1.0;    Id(2,2) = 1.0;  Id(3,3) = 1.0 
      GT     = MATMUL( TRANSPOSE(A), A )              ! Compute G=At*A
      devG   = GT - (GT(1,1)+GT(2,2)+GT(3,3))/3.*Id   ! Compute devG
      CALL compute_total_stress(Yeq,Ydev,Yp, LL_gpr, MM_gpr, GT, devG, detA, ODE)

      ! Compute the stress tensor
       TT      = -detA*MM_gpr*MATMUL(GT,devG) ! changed
       pf=  -(LL_gpr+2.0/3.0*MM_gpr)*(detA)**2*(1-detA)

       syy  = TT(2,2)-pf
       sxy  = TT(1,2)
      ! static friction =0.677*syy 
      ! slp = abs(V(25))
       mus= 0.5
       mud = 0.25
       dc = 250
       mu_t = max(mud,mus-(mus-mud)/dc*(2000.0*t-abs(xGP(1))))

!      if(slp>0.1) print *, 'slp=',slp

       if(sxy.ge.mu_t*abs(syy)) then
          Q(24) = 1000.0*V(21)
        else 
     	  Q(24) = 0.0
        end if

end SUBROUTINE DynamicRupture


RECURSIVE SUBROUTINE PDEFlux(f,g,hz,Q)
  USE MainVariables, ONLY : nVar, nDim, EQN
#if defined(EQNTYPEC99) || defined(EQNTYPED99)
    USE SpecificVarEqn99, only : compute_hydrodynamic_pressure, compute_temperature, computeLMc, compute_l_m_mix!, computeGPRLEpressure
#endif
#if defined(EQNTYPED99)
    USE expintegrator_ode , only : convert_parameters
    use expintegrator_type, only : t_ode_parameters
#endif
#if defined(EQNTYPEB99)
    USE SpecificVarEqn99, only : computeLMc, computeGPRLEpressure, GPRM, computeGPRLEstressGen, compute_hydrodynamic_pressure, compute_l_m_mix
#endif
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: f(nVar), g(nVar), h(nVar), hz(nVar), Q(nVar), V(nVar)
  REAL, PARAMETER :: epsilon = 1e-14 
  INTENT(IN)  :: Q
  INTENT(OUT) :: f, g, hz
  ! Local varialbes
  REAL :: p, A(3,3), AU(3), detA, GT(3,3), devG(3,3), TT(3,3), Id(3,3), T,falpha
  INTEGER :: iErr
  REAL :: ialpha,LEalpha,u(3)
  real :: LL_gpr,MM_gpr,dMM_gpr,dKK_gpr,YY,etaloc
  REAL :: x(3),time,rho,Jx,Jy,Jz, Jv
#ifdef EQNTYPED99
    type(t_ode_parameters) :: ODE
#endif

 if(abs(Q(1)).lt.1)then
     f=0.0
     g=0.0
     hz=0.0
     print *,'flux zero'
     return
  end if

    f=0
    g=0
    h=0
	CALL PDECons2Prim(V,Q,x,time,iErr)

! Extract J and A from Q
    rho  = Q(1)
    Jx   = Q(6)
    Jy   = Q(7)
    Jz   = Q(8)
    A(1,:) = (/ Q( 9), Q(10), Q(11) /) 
    A(2,:) = (/ Q(12), Q(13), Q(14) /) 
    A(3,:) = (/ Q(15), Q(16), Q(17) /) 
    detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 

    LEalpha=Q(18)                                   ! Define diffuse interface alpha from Q

    if(EQN%epsilon1<0) then
        ialpha=1/LEalpha                            ! Take the inverse if EQN%epsilon1<0, this can be used with alpha0>0
    else
        ialpha=LEalpha/(LEalpha**2+EQN%epsilon1*(1-LEalpha))    ! Smooth the inversion process, this allows to take alpha0=0
    end if
    
    u=Q(2:4)/Q(1)*ialpha                            ! Local velocity vector
    AU = MATMUL( A, u(1:3) )                        ! Compute the product of A*u
    ! Initialize the identity matrix needed for the computation of devG
    Id     = 0.     
    Id(1,1) = 1.0;    Id(2,2) = 1.0;  Id(3,3) = 1.0 
    GT     = MATMUL( TRANSPOSE(A), A )              ! Compute G=At*A
    devG   = GT - (GT(1,1)+GT(2,2)+GT(3,3))/3.*Id   ! Compute devG

#ifdef EQNTYPED99
    call convert_parameters(Q, ODE, EQN)
    call compute_l_m_mix(LL_GPR, MM_GPR,Q(19),Q(20),ODE%lam1lam2ratio,ODE%mu1mu2ratio,Q(21))
#else
    call compute_l_m_mix(LL_GPR, MM_GPR,Q(19),Q(20),Q(23),Q(22),Q(21))
#endif

    call compute_hydrodynamic_pressure(p, detA, LL_gpr, MM_gpr,V(5),EQN%gamma,EQN%cv, EQN%p0, V(1), EQN%EOS, EQN%T0)
    call compute_temperature(T, detA, LL_gpr, MM_gpr,rho,EQN%gamma,EQN%cv, EQN%p0,V(5), EQN%EOS, EQN%T0)
 
    TT      = -detA*MM_gpr*MATMUL(GT,devG)       ! Compute the stress matrix TT
    Jv   = u(1)*Jx + u(2)*Jy + u(3)*Jz            ! Compute u*J
    
	f(1)=Q(2)*ialpha
    f(2)=Q(2)*u(1)      +LEalpha*p      -LEalpha*TT(1,1)    + LEalpha*rho*EQN%ch**2*Q(6)*Q(6)
    f(3)=Q(3)*u(1)                      -LEalpha*TT(2,1)    + LEalpha*rho*EQN%ch**2*Q(6)*Q(7)
    f(4)=Q(4)*u(1)                      -LEalpha*TT(3,1)    + LEalpha*rho*EQN%ch**2*Q(6)*Q(8)
    f(5) = u(1)*(Q(5) + p ) + rho*EQN%ch**2*(u(1)*Q(6)*Q(6)+u(2)*Q(6)*Q(7)+u(3)*Q(6)*Q(8)) - u(1)*TT(1,1) - u(2)*TT(2,1) - u(3)*TT(3,1) + rho*T*EQN%ch**2*Jx 
    f(6) = Jv + T 
    f(7) = 0.0 
    f(8) = 0.0 
    f(9 )=AU(1)
    f(12)=AU(2)
    f(15)=AU(3)
    
    g(1)=Q(3)*ialpha
    g(2)=Q(2)*u(2)                      -LEalpha*TT(1,2)   + LEalpha*rho*EQN%ch**2*Q(7)*Q(6)
    g(3)=Q(3)*u(2)      +LEalpha*p      -LEalpha*TT(2,2)   + LEalpha*rho*EQN%ch**2*Q(7)*Q(7)
    g(4)=Q(4)*u(2)                      -LEalpha*TT(3,2)   + LEalpha*rho*EQN%ch**2*Q(7)*Q(8)
    g(5) = u(2)*(Q(5) + p ) + rho*EQN%ch**2*(u(1)*Q(7)*Q(6)+u(2)*Q(7)*Q(7)+u(3)*Q(7)*Q(8)) - u(1)*TT(1,2) - u(2)*TT(2,2) - u(3)*TT(3,2) + rho*T*EQN%ch**2*Jy  
    g(6) = 0.0 
    g(7) = Jv + T 
    g(8) = 0.0     
    g(10)=AU(1)
    g(13)=AU(2)
    g(16)=AU(3)
    g(28)=-V(2) ! D.Li 7.May
    
    h(1)=Q(4)*ialpha
    h(2)=Q(2)*u(3)                      -LEalpha*TT(1,3)    + LEalpha*rho*EQN%ch**2*Q(8)*Q(6)
    h(3)=Q(3)*u(3)                      -LEalpha*TT(2,3)    + LEalpha*rho*EQN%ch**2*Q(8)*Q(7)
    h(4)=Q(4)*u(3)      +LEalpha*p      -LEalpha*TT(3,3)    + LEalpha*rho*EQN%ch**2*Q(8)*Q(8)
    h(5) = u(3)*(Q(5) + p ) + rho*EQN%ch**2*(u(1)*Q(7)*Q(6)+u(2)*Q(7)*Q(7)+u(3)*Q(7)*Q(8)) - u(1)*TT(1,3) - u(2)*TT(2,3) - u(3)*TT(3,3) + rho*T*EQN%ch**2*Jz 
    h(6) = 0.0 
    h(7) = 0.0 
    h(8) = Jv + T  
    h(11)=AU(1)
    h(14)=AU(2)
    h(17)=AU(3)	

  IF(nDim==3) THEN
    hz=h
  END IF  
  END SUBROUTINE PDEFlux


RECURSIVE SUBROUTINE PDENCP(BgradQ,Q,gradQ) 
   USE MainVariables, ONLY :  nVar, nDim, EQN
   IMPLICIT NONE
   ! 11. Oct 21:40: This was a matrix BGradQ(nVar, nDim) but is a vector in spaceTimePredictorNonlinear
   REAL, INTENT(OUT) :: BgradQ(nVar)
   REAL, INTENT(IN)  :: gradQ(nVar, 3)
   REAL, INTENT(IN)  :: Q(nVar)
  ! Linear elasticity variables
   REAL :: u(3),VP(nVar) 
   REAL :: AQx(nVar), BQy(nVar), CQz(nVar) , Qx(nVar), Qy(nVar), Qz(nVar)
	REAL :: alpha, ialpha
 if(abs(Q(1)).lt.1)then
     BgradQ=0.0
     print *, 'NCP zero'
     return
  end if
	
   
	Qx = gradQ(:,1)
	Qy = gradQ(:,2)
	IF(nDim==3) THEN
		Qz = gradQ(:,3)
	ELSE
		Qz = 0.0 
	ENDIF 
     ! Initialize NCP ---!
    AQx=0.              !
    BQy=0.              !
    CQz=0.              !
    ! ------------------!
	alpha=Q(18) ! GPR+LE+DI
    if(EQN%epsilon1<0) then
        ialpha=1/alpha
    else
        ialpha=alpha/(alpha**2+EQN%epsilon1*(1-alpha))
    end if
    u=Q(2:4)/Q(1)*ialpha
 
 ! NCP for the J vector
    AQx(6)  =  -u(2)*Qx(7) - u(3)*Qx(8)
    AQx(7)  =  +u(1)*Qx(7)
    AQx(8)  =  +u(1)*Qx(8)
    ! NCP for the A matrix
	AQx(9)  =  -u(2)*Qx(10) - u(3)*Qx(11)
	AQx(10) =  +u(1)*Qx(10)
	AQx(11) =  +u(1)*Qx(11)
	AQx(12) =  -u(2)*Qx(13) - u(3)*Qx(14)
	AQx(13) =  +u(1)*Qx(13)
	AQx(14) = +u(1)*Qx(14)
	AQx(15) = -u(2)*Qx(16) - u(3)*Qx(17)
	AQx(16) = +u(1)*Qx(16)
    AQx(17) = +u(1)*Qx(17)
    
    ! NCP for the J vector
    BQy(6 ) =  +u(2)*Qy(6)
    BQy(7 ) =  -u(1)*Qy(6) - u(3)*Qy(8)
    BQy(8 ) =  +u(2)*Qy(8)
    ! NCP for the A matrix
    BQy(9 ) =  +u(2)*Qy(9)
    BQy(10) =  -u(1)*Qy(9) - u(3)*Qy(11)
    BQy(11) =  +u(2)*Qy(11)
    BQy(12) =  +u(2)*Qy(12)
    BQy(13) = -u(1)*Qy(12) - u(3)*Qy(14)
    BQy(14) = +u(2)*Qy(14)
    BQy(15) = +u(2)*Qy(15)
    BQy(16) = -u(1)*Qy(15) - u(3)*Qy(17)
    BQy(17) = +u(2)*Qy(17)
    
    ! NCP for the J vector
    CQz(6 ) =  +u(3)*Qz(6)
    CQz(7 ) =  +u(3)*Qz(7)
    CQz(8 ) =  -u(1)*Qz(6) - u(2)*Qz(7)
    ! NCP for the A matrix
    CQz(9 ) =  +u(3)*Qz(9)
    CQz(10) =  +u(3)*Qz(10)
    CQz(11) =  -u(1)*Qz(9) - u(2)*Qz(10)
    CQz(12) =  +u(3)*Qz(12)
    CQz(13) = +u(3)*Qz(13) 
    CQz(14) = -u(1)*Qz(12) - u(2)*Qz(13)
    CQz(15) = +u(3)*Qz(15) 
    CQz(16) = +u(3)*Qz(16)
    CQz(17) = -u(1)*Qz(15) - u(2)*Qz(16)
   
    ! Transport equation
    AQx(18:24-1) = +u(1)*Qx(18:24-1)
    BQy(18:24-1) = +u(2)*Qy(18:24-1)
    CQz(18:24-1) = +u(3)*Qz(18:24-1)
    !AQx(25:27)   = +u(1)*Qx(25:27)
    !BQy(25:27)   = +u(2)*Qy(25:27)
    !CQz(25:27)   = +u(3)*Qz(25:27)
 
    if( nDim .eq. 2) then
        BgradQ = AQx + BQy         
    else
        BgradQ = AQx + BQy + CQz     
    end if
END SUBROUTINE PDENCP


RECURSIVE SUBROUTINE PDEEigenvalues(L,Q,n)
  USE MainVariables, ONLY :  nVar, nDim, EQN
  USE iso_c_binding
  IMPLICIT NONE
  REAL :: L(nVar), n(3), Q(nVar), Vp(nVar)
  INTENT(IN)  :: Q,n
  INTENT(OUT) :: L 
  ! Local Variables 
  REAL :: lam, mu, irho, VPR(3), cs,c0,uu
   
  if(abs(Q(1)).lt.1)then 
     L=1
     print *, 'eigenvalue Zero'
     return
  end if

   lam  = Q(19)   
   mu   = Q(20)
   irho = 1./Q(1)
  
   L    = 0 
   VPR=Q(2:4)/Q(1)  ! Skip alpha for now!
   uu = SQRT(sum(VPR(:)**2) ) 
   cs=sqrt(mu*irho)
   c0 = SQRT((lam+2*mu)*irho + EQN%ch**2 )
   L(1)=uu+c0
   L(2)=uu-c0
   L(3)=uu+cs
   L(4)=uu-cs
   L(5)=uu+cs
   L(6)=uu-cs
   
   L(18:24-1)=uu
   L(25:27)=0.0 !  D.Li slip(U,V,W)

END SUBROUTINE PDEEigenvalues

RECURSIVE SUBROUTINE PDESource(S,Q) 
  USE MainVariables, ONLY:  nVar, nDim,EQN
#if defined(EQNTYPEB99) || defined(EQNTYPEA99) || defined(EQNTYPEC99) || defined(EQNTYPED99)
  USE SpecificVarEqn99
#endif
#if defined(EQNTYPED99)
    USE expintegrator_ode , only : convert_parameters
    use expintegrator_type, only : t_ode_parameters
#endif
  USE iso_c_binding
  IMPLICIT NONE
  ! --------------------------------------------
  ! Argument list declaration
  REAL :: S(nvar), Q(nvar)
  INTENT(IN)  :: Q 
  INTENT(OUT) :: S
  ! --------------------------------------------
#if defined(EQNTYPEB99) || defined(EQNTYPEC99)  || defined(EQNTYPED99)
    real :: LLc,MMc,dMM,dKK,Pressure,YY_sh,YY_p, YY, tau1,tau2,taum,dEc, Kelast
#endif
#ifdef EQNTYPED99
    type(t_ode_parameters) :: ODE
#endif
  ! Local variable declaration 
  REAL :: AM(3,3), Id(3,3), G(3,3) ,devG(3,3), detA, detA2,psiM(3,3), temp2, T
  REAL :: V(nvar)
  REAL :: X(3),time, theta2
  INTEGER :: iErr

 if(abs(Q(1)).lt.1)then
     S=0.0
     print *, 'source zero'
     return
  end if

  S = 0.
! --------------- SOURCE CONTRIBUTION TO A ----------------------------------
#ifndef ODESOLVER 
    ! Recontruct the matrix A
    AM(1,:) = (/ Q( 9), Q(10), Q(11) /) 
    AM(2,:) = (/ Q(12), Q(13), Q(14) /)
    AM(3,:) = (/ Q(15), Q(16), Q(17) /)
    G      = MATMUL( TRANSPOSE(AM), AM )
    Id     = 0. 
    Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0 
    devG   = G - (G(1,1)+G(2,2)+G(3,3))/3.*Id
    detA   = AM(1,1)*AM(2,2)*AM(3,3)-AM(1,1)*AM(2,3)*AM(3,2)-AM(2,1)*AM(1,2)*AM(3,3)+AM(2,1)*AM(1,3)*AM(3,2)+AM(3,1)*AM(1,2)*AM(2,3)-AM(3,1)*AM(1,3)*AM(2,2)
    psiM   = 3./(detA)*MATMUL(AM,devG)
    
    tau =EQN%tau0
    
    S(9:17) = -(/ psiM(1,1), psiM(1,2), psiM(1,3), & 
                  psiM(2,1), psiM(2,2), psiM(2,3), & 
                  psiM(3,1), psiM(3,2), psiM(3,3) /) /tau*detA**(8./3.)*Q(18)
#endif
! ----------------------------------------------------------------------------

   ! ------------------------------------------------------------------------------------------------------
   ! Do not use ODESolver for J! -----------------------------------------------------------------------
   ! For debug remove the source for J
   CALL PDECons2Prim(V,Q,x,time,iErr)
   AM(1,:) = (/ Q( 9), Q(10), Q(11) /) 
   AM(2,:) = (/ Q(12), Q(13), Q(14) /) 
   AM(3,:) = (/ Q(15), Q(16), Q(17) /) 
   detA   = AM(1,1)*AM(2,2)*AM(3,3)-AM(1,1)*AM(2,3)*AM(3,2)-AM(2,1)*AM(1,2)*AM(3,3)+AM(2,1)*AM(1,3)*AM(3,2)+AM(3,1)*AM(1,2)*AM(2,3)-AM(3,1)*AM(1,3)*AM(2,2)     ! this is the determinant we have     
   
   G      = MATMUL( TRANSPOSE(AM), AM ) 
   Id     = 0. 
   Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0 
   devG   = G - (G(1,1)+G(2,2)+G(3,3))/3.*Id   
   
   !call ComputeLMc(LLc,MMc,dMM,dKK,Q(20),Q(23),Q(19),Q(24),Q(21))
#ifdef EQNTYPED99
    call convert_parameters(Q, ODE, EQN)
    call compute_l_m_mix(LLc, MMc,Q(19),Q(20),ODE%lam1lam2ratio,ODE%mu1mu2ratio,Q(21))
#else
    call compute_l_m_mix(LLc, MMc,Q(19),Q(20),Q(23),Q(22),Q(21))
#endif
   ! Compute the pressure as p=rho^2*E_rho -------------------------------------------------------------------------------
   Kelast=LLc+2.0/3.0*MMc
   !p = Kelast/EQN%gamma*(detA)**EQN%gamma*EXP(V(5)/EQN%cv) - Kelast/EQN%gamma + EQN%p0
   !T = 1.0/EQN%cv/(EQN%gamma-1.0)/EQN%gamma*Kelast/Q(1)*detA**EQN%gamma
   !call compute_hydrodynamic_pressure(p, detA, LLc, MMc,V(5),EQN%gamma,EQN%cv, EQN%p0, V(1), EQN%EOS, EQN%T0)
   call compute_temperature(T, detA, LLc, MMc,Q(1),EQN%gamma,EQN%cv, EQN%p0,V(5), EQN%EOS, EQN%T0)
   ! Other way to compute the pressure -------------------------------------------------------------------------------------------------------
   !tempp  = MATMUL( TRANSPOSE(devG), devG )
   !p    = (EQN%gamma-1.0)*( Q(5) - 0.5*rho*sum(u**2) - 0.5*rho*EQN%ch**2*(Jx**2+Jy**2+Jz**2) - MMc/4.*(tempp(1,1)+tempp(2,2)+tempp(3,3)) )
   !T    = 1.0/EQN%cv/(EQN%gamma-1.0) * p/rho
   ! -----------------------------------------------------------------------------------------------------------------------------------------
   
   theta2 = EQN%tau2 / ( Q(1) * T ) 
   S(6:8) = -Q(6:8)/theta2          
   
   ! ------------------------------------------------------------------------------------------------------
   ! Computer other sources
   S(2:4)=Q(1)*Q(18)*EQN%gvec(1:3)
   S(5)=Q(1)*dot_product(V(2:4),EQN%gvec(1:3))

  ! Duo April 9 - linear or exponential; S(23)=0 -> heviside func. 
   S(23) = -Q(24)!*Q(23)
   S(25:27) = V(2:4)  ! D.Li slip source term

END SUBROUTINE PDESource

RECURSIVE SUBROUTINE PDEVarName(MyNameOUT,ind) 
  USE MainVariables, ONLY: nVar  
  IMPLICIT NONE     
  CHARACTER(LEN=10):: MyName(nVar),MyNameOUT
  INTEGER			:: ind

    MyName(1 ) = 'rho' 
    MyName(2 ) = 'u' 
    MyName(3 ) = 'v'
    MyName(4 ) = 'w' 
      
    MyName(5 ) = 's'
    MyName(6 ) = 'Jx'
    MyName(7 ) = 'Jy'
    MyName(8 ) = 'Jz'
          
    MyName(9 ) = 'A11' 
    MyName(10) = 'A12' 
    MyName(11) = 'A13' 
    MyName(12) = 'A21' 
    MyName(13) = 'A22' 
    MyName(14) = 'A23'
    MyName(15) = 'A31'
    MyName(16) = 'A32' 
    MyName(17) = 'A33'
    
    MyName(18) = 'alpha'
    MyName(19) = 'lambda'
    MyName(20) = 'mu'
    MyName(21) = 'xi'
    MyName(22) = 'xi1'
    MyName(23) = 'sliding'
    MyName(24) = 'sliding_xi'
    MyName(25) = 'U'
    MyName(26) = 'V'
    MyName(27) = 'W'
    MyName(28) = 'derivY'
	
	MyNameOUT=MyName(ind+1)
    END SUBROUTINE PDEVarName

RECURSIVE SUBROUTINE PDEAuxName(MyNameOUT,ind) 
	USE MainVariables, ONLY: nVar,nAux  
	IMPLICIT NONE     
	CHARACTER(LEN=10):: AuxName(nAux),MyNameOUT
	INTEGER			:: ind


	  ! EQNTYPE99
	  AuxName(1) = 'sxx'
	  AuxName(2) = 'syy'
	  AuxName(3) = 'szz'
	  AuxName(4) = 'sxy'
	  AuxName(5) = 'syz'
	  AuxName(6) = 'sxz'
	  AuxName(7) = 'YY'
	  AuxName(8) = 'YYsh'
	  AuxName(9) = 'YYp'
	  AuxName(10) = 'tau1'
	  AuxName(11) = 'tau2'
	  AuxName(12) = 'LL'
	  AuxName(13) = 'MM'
	  AuxName(14) = 'T'
	  AuxName(15) = 'diff|ADG|'
	  AuxName(16) = 'equiMu'

	MyNameOUT=AuxName(ind+1)
END SUBROUTINE PDEAuxName
	
RECURSIVE SUBROUTINE PDEAuxVar(aux,Q,x,time)
    USE MainVariables, ONLY : nVar,nAux,EQN
#if defined(EQNTYPEC99) || defined(EQNTYPED99)
    use SpecificVarEqn99, only : Compute_l_m_mix, compute_total_stress, compute_temperature, compute_hydrodynamic_pressure, computeGPRLEstressGen
    use expintegrator_type, only : t_ode_parameters
    use expintegrator_ode, only : convert_parameters
#endif
	implicit none
#if defined(ODESOLVER) || defined(EQNTYPEC99) || defined(EQNTYPEB99) || defined(EQNTYPED99)
  real ::V0(10),G(3,3),Ydev,Yp, tau1,tau2
  type(t_ode_parameters) :: ODE
#endif
	real :: aux(nAux),Q(nVar),x(3),time
	REAL :: LL_gpr, MM_gpr
	real :: detA,A(3,3), Id(3,3),GT(3,3),devG(3,3),TT(3,3)
	integer :: iErr
	REAL :: V(nVar)
	! print *, "ENTERED HERE ------------------------------------------------------------------------------"
	aux=0.
	!return
	CALL PDECons2Prim(V,Q,x,0.,iErr)
	
#ifdef EQNTYPED99
if(V(21)>1.e-2) then
    continue
    end if
    call convert_parameters(V, ODE, EQN)
    call compute_l_m_mix(LL_GPR, MM_GPR,V(19),V(20),ODE%lam1lam2ratio,ODE%mu1mu2ratio,V(21))
#else
    call compute_l_m_mix(LL_GPR, MM_GPR,Q(19),Q(20),Q(23),Q(22),Q(21))
#endif
    !call compute_l_m_mix(LL_GPR, MM_GPR,Q(19),Q(20),Q(23),Q(22),Q(21))
    aux(12)=LL_gpr      ! Mixture parameter lambda
    aux(13)=MM_gpr      ! Mixture parameter mu
      
    A(1,:) = (/ Q( 9), Q(10), Q(11) /) 
    A(2,:) = (/ Q(12), Q(13), Q(14) /) 
    A(3,:) = (/ Q(15), Q(16), Q(17) /) 
    detA   = A(1,1)*A(2,2)*A(3,3)-A(1,1)*A(2,3)*A(3,2)-A(2,1)*A(1,2)*A(3,3)+A(2,1)*A(1,3)*A(3,2)+A(3,1)*A(1,2)*A(2,3)-A(3,1)*A(1,3)*A(2,2) 
    Id     = 0.     
    Id(1,1) = 1.0;    Id(2,2) = 1.0;  Id(3,3) = 1.0 
    GT     = MATMUL( TRANSPOSE(A), A )              ! Compute G=At*A
    devG   = GT - (GT(1,1)+GT(2,2)+GT(3,3))/3.*Id   ! Compute devG

! Generalized version
    call compute_hydrodynamic_pressure(aux(9), detA, LL_gpr, MM_gpr,V(5),EQN%gamma,EQN%cv, EQN%p0, Q(1), EQN%EOS, EQN%T0)
    call compute_temperature(aux(14), detA, LL_gpr, MM_gpr,Q(1),EQN%gamma,EQN%cv, EQN%p0,V(5), EQN%EOS, EQN%T0)

    CALL convert_parameters(V, ODE, EQN)
    CALL compute_total_stress(aux(7),Ydev,Yp, LL_gpr, MM_gpr, GT, devG, detA, ODE)
    aux(9)=-(LL_gpr+2.0/3.0*MM_gpr)*(detA)**2*(1-detA)
    ! Compute the stress tensor
    TT      = -detA*MM_gpr*MATMUL(GT,devG) ! changed
    aux(1)  = TT(1,1)-aux(9)
    aux(2)  = TT(2,2)-aux(9)
    aux(3)  = TT(3,3)-aux(9)
    aux(4)  = TT(1,2)
    aux(5)  = TT(2,3)
    aux(6)  = TT(1,3)
    call computeGPRLEstressGen(aux(7),aux(1:6),Q(9:17), LL_gpr,MM_gpr, 0.,0.,.true.)
    
    aux(9)=Yp
    aux(8)=Ydev
    
    tau1 = ODE%tau0*exp(min(ODE%alpha1 - ODE%beta1*(1.-V(21))*aux(7), 690.0d0))  
    tau2 = ODE%tau0*exp(min(ODE%alpha2 - ODE%beta2*V(21)*aux(7), 690.0d0)) 
    ! Compute the relaxation time for the mixture; itaum = c2/mu1*tau1 + c1/(mu2*tau2)
    aux(10) = min(1.e+16,(ODE%mu1*tau1*ODE%mu2*tau2)/(V(21)*ODE%mu1*tau1 + (1.-V(21))*ODE%mu2*tau2))
    aux(11) =EQN%tau2
    aux(15) = maxval(abs(matmul(A,devG)-matmul(devG,A)))
 
   ! equivalent mu
    aux(16) = max(-1e+37,min(1e+37,aux(4)/abs(aux(2))) ) 

END SUBROUTINE PDEAuxVar
	
	
RECURSIVE SUBROUTINE PDECritialStress(CS,Q)
    USE MainVariables, ONLY : nVar,nAux,EQN
#if defined(EQNTYPEC99) || defined(EQNTYPED99)
    use SpecificVarEqn99        , only : StaticLimiterEQ99,compute_l_m_mix, compute_total_stress, compute_hydrodynamic_pressure
    use expintegrator_linalg    , only : determinant_lu,build_by_rows,trace 
    use expintegrator_type      , only : t_ode_parameters
    use expintegrator_ode       , only : convert_parameters
#endif
	implicit none
#if defined(EQNTYPEC99) || defined(EQNTYPED99)
  real      :: stressnorm,sigma_vec(6),dMM_gpr,dKK_gpr, Ydev,Yp
  real      :: V0(nVar),G(3,3)
  type(t_ode_parameters) :: ODE
#endif
	real :: CS,Q(nVar),x(3),time
	REAL :: LL_gpr, MM_gpr
	real :: detA,A(3,3), Id(3,3),GT(3,3),devG(3,3),TT(3,3)
	integer :: iErr
	REAL :: V(nVar)
	! print *, "ENTERED HERE ------------------------------------------------------------------------------"
	CS=0.
	x=0.
	!return
	ID=0.
	ID(1,1)=1.
	ID(2,2)=1.
	ID(3,3)=1.
	A = build_by_rows(Q(9:17))
	G = matmul(transpose(A), A)
	devG  = G - trace(G)/3.0d0*ID
	detA = determinant_lu(A)
#ifdef EQNTYPED99
        CALL PDECons2Prim(V0,Q,x,0.,iErr)
        call convert_parameters(V0, ODE, EQN)
        call compute_l_m_mix(LL_GPR, MM_GPR,V0(19),V0(20),ODE%lam1lam2ratio,ODE%mu1mu2ratio,V0(21))
        CALL compute_total_stress(stressnorm,Ydev,Yp, LL_gpr, MM_gpr, G, devG, detA, ODE)
        IF(stressnorm>ODE%Y0*EQN%Y0lim .and. V0(18)>1.e-4) then      ! stress>Y0/0.9 V0(25)*1.0 ! .or. V0(5)<-1.e-2
            CS=1.
        end if
#else
        CALL PDECons2Prim(V0,Q,lxb,0.,iErr)
        call convert_parameters(V0, ODE, EQN)
        call compute_l_m_mix(LL_GPR, MM_GPR,V0(19),V0(20),V0(23),V0(22),V0(21))
        CALL compute_total_stress(stressnorm,Ydev,Yp, LL_gpr, MM_gpr, G, devG, detA, ODE)
        IF(stressnorm>V0(25)*EQN%Y0lim) then      ! stress>Y0/0.9 V0(25)*1.0 ! .or. V0(5)<-1.e-2
            CS=1.
        end if
#endif
    
END SUBROUTINE PDECritialStress
	
RECURSIVE SUBROUTINE PDEIntermediateFields(RL,LL,iRL,Q,nv) 
  USE MainVariables, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  IMPLICIT NONE
  ! Argument list 
  REAL :: RL(nVar,nLin), LL(nLin,nLin), iRL(nLin,nVar)
  REAL :: Q(nVar), nv(3)
  REAL :: x(nDim), time
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: RL,LL,iRL 
  ! Local variables
  INTEGER :: i,j,k, zero,iErr, minl(1), maxl(1) 
  REAL :: R(nVar,nVar), L(nVar,nVar), Lv(nVar), iR(nVar,nVar)
  !  
  RL=0.
  LL=0.
  iRL=0.
  CALL PDEEigenvectors(R,L,iR,Q,nv)
#if defined(EQNTYPEC99)
    CALL PDEEigenvectors(R,L,iR,Q,nv)
    RL(:,1:17)  = R(:,18:34)
    iRL(1:17,:) = iR(18:34,:)
    LL = 0.
    do i=1,17
        LL(i,i) = L(17+i,17+i)
    end do
#elif defined(EQNTYPED99)
    CALL PDEEigenvectors(R,L,iR,Q,nv)
	do i=1,nVar
		RL(i,1:7)  = R(i,18:24)
		iRL(1:7,i) = iR(18:24,i)
	end do
    LL = 0.
    do i=1,7
        LL(i,i) = L(17+i,17+i)
    end do
  !
#endif
END SUBROUTINE PDEIntermediateFields
    
    
RECURSIVE SUBROUTINE PDEEigenvectors(R,L,iR,Q,nv) 
	USE MainVariables, ONLY : nVar, nDim, EQN
	USE iso_c_binding  
	IMPLICIT NONE
	! Argument list 
	REAL :: R(nVar,nVar), L(nVar,nVar), iR(nVar,nVar)
	REAL :: Q(nVar), nv(3)
	REAL :: x(nDim), time
	INTENT(IN)  :: Q,nv
	INTENT(OUT) :: R,L,iR 
	! Local variables
	REAL    :: Lambda(nVar),sv(3),tv(3),TM(3,3),iTM(3,3),Qrot(3),ialpha,uv(3)
    INTEGER :: j,i    

	Lambda=0.
	R=0.
	L = 0.
	iR=0.
	do j=1,nVar ! D.Li April 27
		R(j,j)=1.
		iR(j,j)=1.
	end do

!   Add velocity for eigenvalues
    ! 
    IF(ABS(ABS(nv(1))-1.0).LE.1e-14) THEN 
        sv = (/ 0., 1., 0. /) 
        tv = (/ 0., 0., 1. /) 
    ENDIF 
    IF(ABS(ABS(nv(2))-1.0).LE.1e-14) THEN 
        sv = (/ 0., 0., 1. /) 
        tv = (/ 1., 0., 0. /) 
    ENDIF 
    IF(ABS(ABS(nv(3))-1.0).LE.1e-14) THEN 
        sv = (/ 1., 0., 0. /) 
        tv = (/ 0., 1., 0. /) 
    ENDIF
    TM(1,1:3) = (/ nv(1), sv(1), tv(1) /)     
    TM(2,1:3) = (/ nv(2), sv(2), tv(2) /)     
    TM(3,1:3) = (/ nv(3), sv(3), tv(3) /)  
    iTM(1,1:3) = (/ nv(1), nv(2), nv(3) /)     
    iTM(2,1:3) = (/ sv(1), sv(2), sv(3) /)     
    iTM(3,1:3) = (/ tv(1), tv(2), tv(3) /) 
        
    Qrot(2:4) = MATMUL( iTM(1:3,1:3), Q(2:4))
    ialpha=Q(18)/(Q(18)**2+EQN%epsilon1*(1-Q(18)))
    uv=Qrot(2:4)*1.0/Q(1)  
    Lambda(18:24-1)=uv(1)
   ! Lambda(25:27)=uv(1)
    do i=18,24-1
        L(i,i)=Lambda(i)
    end do
    !do i = 25,27
    !   L(i,i)=Lambda(i)
    !end do
    END SUBROUTINE PDEEigenvectors 
	
RECURSIVE SUBROUTINE PDEJacobian(An,Q,gradQ,nv)
  USE MainVariables, ONLY : nVar, nDim
  USE iso_c_binding
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), gradQ(nVar,nDim), nv(3) 
  INTENT(IN)  :: Q, gradQ, nv
  INTENT(OUT) :: An 
  
  CALL PDEMatrixB(An,Q,nv) 
  
END SUBROUTINE PDEJacobian
	
RECURSIVE SUBROUTINE PDEMatrixB(An,Q,nv) 
  USE MainVariables, ONLY : nVar, nDim, EQN
  USE iso_c_binding
  IMPLICIT NONE
  ! Argument list 
  REAL :: An(nVar,nVar)
  REAL :: Q(nVar), nv(3) 
  INTENT(IN)  :: Q,nv
  INTENT(OUT) :: An  
  ! Local variables
  ! Linear elasticity variables
  REAL :: A(nVar,nVar), B(nVar,nVar), C(nVar,nVar), Vp(nVar)
  REAL :: alpha , ialpha, uv(3)
  INTEGER :: i
   
    A=0.
    B=0.
    C=0.
    alpha=Q(18) ! GPR+LE+DI
    if(EQN%epsilon1<0) then
        ialpha=1/alpha
    else
        ialpha=alpha/(alpha**2+EQN%epsilon1*(1-alpha))
    end if  
	uv=Q(2:4)/Q(1)*ialpha
	
! J contribution
    A(6,7)  =   -uv(2)
    A(6,8)  =   -uv(3)
    A(7,7)  =   +uv(1)
    A(8,8)  =   +uv(1)
    ! A contribution
    A(9,10)  =  -uv(2)
    A(9,11)  =  -uv(3)
	A(10,10) =  +uv(1)
	A(11,11) =  +uv(1)
	A(12,13) =  -uv(2)
    A(12,14) =  -uv(3)
	A(13,13) =  +uv(1)
	A(14,14) =  +uv(1)
	A(15,16) =  -uv(2)
    A(15,17) =  -uv(3)
	A(16,16) =  +uv(1)
    A(17,17) =  +uv(1)
    ! Transport contribution of the parameters
    do i=18,24-1
        A(i,i) =  +uv(1)    
    end do
    
    ! J contribution
    B(6,6)   =  +uv(2)
    B(7,6)   =  -uv(1)
    B(7,8)   =  -uv(3)
    B(8,8)   =  +uv(2)
    ! A contribution
    B(9 ,9)  =  +uv(2)
    B(10 ,9) =  -uv(1)
    B(10 ,11)=  -uv(3)
    B(11 ,11)=  +uv(2)
    B(12 ,12)=  +uv(2)
    B(13 ,12)=  -uv(1)
    B(13,14) =  -uv(3)
    B(14,14) =  +uv(2)
    B(15,15) =  +uv(2)
    B(16,15) =  -uv(1)
    B(16,17) =  -uv(3)
    B(17,17) =  +uv(2)  
    ! Transport contribution of the parameters
    do i=18,24-1
        B(i,i) =  +uv(2)    
    end do
    
    ! J contribution
    C(6 ,6)  =  +uv(3)
    C(7 ,7)  =  +uv(3)
    C(8 ,6)  =  -uv(1)
    C(8 ,7)  =  -uv(2)
    ! A contribution
    C(9 ,9)  =  +uv(3)
    C(10 ,10)=  +uv(3)
    C(11 ,9) =  -uv(1)
    C(11 ,10)=  -uv(2)
    C(12 ,12)=  +uv(3)
    C(13 ,13)=  +uv(3)
    C(14,12) =  -uv(1)
    C(14,13) =  -uv(2)
    C(15,15) =  +uv(3) 
    C(16,16) =  +uv(3)
    C(17,15) =  -uv(1) 
    C(17,16) =  -uv(2)
    ! Transport contribution of the parameters
    do i=18,24-1
        C(i,i) =  +uv(3)    
    end do	

   !do i = 25,27
   ! A(i,i)=uv(1)
   ! B(i,i)=uv(2)
   ! C(i,i)=uv(3)
   ! end do	
	
	
    if( nDim .eq. 2) then
        An = A*nv(1) + B*nv(2)         
    else
        An = A*nv(1) + B*nv(2) + C*nv(3)     
    end if

END SUBROUTINE PDEMatrixB

RECURSIVE SUBROUTINE RoeMatrix(ARoe,QL,QR,nv) 
  USE MainVariables, ONLY : nVar, nDim,sGP3,wGP3,nGP3
  IMPLICIT NONE
  ! Argument list 
  REAL        :: ARoe(nVar,nVar), QL(nVar), QR(nVar), nv(3) 
  INTENT(IN)  :: QL, QR, nv  
  INTENT(OUT) :: ARoe 
  ! Local variables 
  INTEGER             :: iGP  
  REAL                :: psi(nVar) 
  REAL                :: A(nVar,nVar)
  REAL                :: gradQ(nVar,nDim) 
  ! Midpoint Rule 
  !REAL, PARAMETER     :: sGP(1) = (/ 0.5 /) 
  !REAL, PARAMETER     :: wGP(1) = (/ 1.0 /) 
  !INTEGER, PARAMETER :: nGP = 1 
  ! Trapezoidal Rule 
  !REAL, PARAMETER     :: sGP(2) = (/ 0.0, 1.0 /) 
  !REAL, PARAMETER     :: wGP(2) = (/ 0.5, 0.5 /) 
  !INTEGER, PARAMETER :: nGP = 2 
  ! 4-point Gaussian quadrature 
  !REAL, PARAMETER     :: sGP(4) = (/ 1.D0/2.D0-sqrt(525+70*sqrt(30.D0))/70,  1.D0/2.D0-sqrt(525-70*sqrt(30.D0))/70,  1.D0/2.D0+sqrt(525-70*sqrt(30.D0))/70, 1.D0/2.D0+sqrt(525+70*sqrt(30.D0))/70 /) 
  !REAL, PARAMETER     :: wGP(4) = (/ 1.D0/4.D0-sqrt(30.D0)/72, 1.D0/4.D0+sqrt(30.D0)/72,  1.D0/4.D0+sqrt(30.D0)/72,  1.D0/4.D0-sqrt(30.D0)/72 /) 
  !INTEGER, PARAMETER  :: nGP = 4  
  !
  gradQ = 0.0 
  !
  ARoe = 0. 
  DO iGP = 1, nGP3  
     psi = QL + sGP3(iGP)*(QR-QL)  
     CALL PDEJacobian(A,psi,gradQ,nv) 
     ARoe = ARoe + wGP3(iGP)*A   ! Numerical integration of the Roe-matrix  
  ENDDO
  !
END SUBROUTINE RoeMatrix 

!!!RECURSIVE SUBROUTINE HLLEMFluxFV(FL,FR,QL,QR,QavL,QavR,NormalNonZero) 
!!!  USE MainVariables, ONLY : nVar, nDim, nLin
!!!  USE iso_c_binding 
!!!  ! Local variables
!!!  INTEGER, INTENT(IN)   :: NormalNonZero
!!!  REAL, INTENT(IN)     :: QL(nVar)
!!!  REAL, INTENT(IN)     :: QR(nVar)
!!!  REAL, INTENT(INOUT)  :: FL(nVar)
!!!  REAL, INTENT(INOUT)  :: FR(nVar)
!!!  REAL    :: QavL(nVar), QavR(nVar)  
!!!    ! Local variables 
!!!  INTEGER           :: i,j,k,l, ml(1)  ,iErr
!!!  REAL              :: smax, Qav(nVar)
!!!  REAL              ::  nv(nDim), flattener(nLin)
!!!  REAL    :: absA(nVar,nVar), amax  
!!!  REAL    :: QM(nVar),LL(nVar),LR(nVar),LM(nVar)
!!!  REAL    :: deltaL(nLin,nLin),Lam(nLin,nLin),Lap(nLin,nLin) 
!!!  REAL    :: RL(nVar,nLin),iRL(nLin,nVar),LLin(nLin,nLin) 
!!!  REAL    :: Aroe(nVar,nVar),Aroep(nVar,nVar), Aroem(nVar,nVar), Dm(nVar), Dp(nVar), dQ(nVar)
!!!  REAL :: f1R(nVar), g1R(nVar), h1R(nVar) 
!!!  REAL :: f1L(nVar), g1L(nVar), h1L(nVar) , VM(nVar) 
!!!  
!!!  REAL :: XX0(3),TIME0
!!!  XX0=0.
!!!  TIME0=0.
!!!  !  
!!!  nv(:)=0.
!!!  nv(NormalNonZero+1)=1.
!!!  !
!!!  flattener=1.0 !0.8
!!!  !
!!!  CALL PDEFlux(f1L,g1L,h1L,QL)
!!!  CALL PDEFlux(f1R,g1R,h1R,QR)
!!!  !
!!!  fR = f1R*nv(1)+g1R*nv(2)+h1R*nv(3)
!!!  fL = f1L*nv(1)+g1L*nv(2)+h1L*nv(3)
!!!  !
!!!!IF(ANY(fR(6:8).NE.0)) THEN
!!!!     PRINT *,"f1R",f1R
!!!!     PRINT *,"g1R",g1R
!!!!     PRINT *,"h1R",h1R
!!!!     PRINT *,"f1L",f1L
!!!!     PRINT *,"g1L",g1L
!!!!     PRINT *,"h1L",h1L
!!!!     PRINT *,"dQ",dQ
!!!!     PRINT *,"QR",QR
!!!!     PRINT *,"QL",QL
!!!!     PRINT *,"dQ",dQ
!!!!     PRINT *,"BEFORE: fR(6:8).NE.0",fR(6:8)
!!!!    STOP
!!!!ENDIF
!!!!IF(ANY(fl(6:8).NE.0)) THEN
!!!!     PRINT *,"QR",QR
!!!!     PRINT *,"QL",QL
!!!!     PRINT *,"dQ",dQ
!!!!     PRINT *,"BEFORE: fl(6:8).NE.0",fl(6:8)
!!!!    STOP
!!!!ENDIF
!!!  !USE Parameters, ONLY : d,nVar,ndim 
!!!  QM = 0.5*(QL+QR) 
!!!  CALL PDECons2Prim(VM,QM,XX0,TIME0,iErr)
!!!  !CALL PDECons2PrimGRMHD(VM,QM,iErr)
!!!  CALL PDEEigenvalues(LL,QL,nv)  
!!!  CALL PDEEigenvalues(LR,QR,nv)  
!!!!IF(ANY(QM(6:8).NE.0)) THEN
!!!!    PRINT *, "HLLEMFluxFV QM(6:8)",QM(6:8)
!!!!    STOP
!!!!ENDIF
!!!  CALL PDEEigenvalues(LM,QM,nv)  
!!!  sL = MIN( 0., MINVAL(LL(:)), MINVAL(LM(:)) ) 
!!!  sR = MAX( 0., MAXVAL(LR(:)), MAXVAL(LM(:)) ) 
!!! ! PRINT *, "PDEIntermediateFields"
!!!  !DO i=1,nVar
!!!  !  WRITE(*,'(E16.6)'), QM(i)
!!!  !ENDDO
!!!  !print *,"*********************************************"
!!!  !WRITE(*,'(a,f18.10,f18.10,f18.10)')    "***** nv:",nv(1),nv(2),nv(3)
!!!  CALL PDEIntermediateFields(RL,LLin,iRL,QM,nv) 
!!!  !PRINT *, "PDEIntermediateFields finished"
!!!  Lam = 0.5*(LLin-ABS(LLin))
!!!  Lap = 0.5*(LLin+ABS(LLin)) 
!!!  deltaL = 0.0
!!!  !print *,"*********************************************"
!!!  !!print *,"*****LLin, QR(1),nv",QM(1),NormalNonZero
!!!  !WRITE(*,'(a,E16.6,i9)')    "***** LLin, QR(1),nv",QM(1),NormalNonZero
!!!  !print *,"**********************************************"
!!!
!!!  DO i = 1, nLin
!!!      deltaL(i,i) = (1. - Lam(i,i)/(sL-1e-14) - Lap(i,i)/(sR+1e-14) )*flattener(i)  
!!!      !print *,"i,DeltaL(i,i):",i,DeltaL(i,i)
!!!      !WRITE(*,'(a,i9,E16.6,E16.6,E16.6,E16.6,E16.6)') "i,QM(i),VM(i),DeltaL(i,i),LM(i):",i,QM(i),VM(i),DeltaL(i,i),LLin(i,i),LM(i)
!!!      !WRITE(*,'(a, i9, f18.10, f16.5, f16.5, i9)') 
!!!  ENDDO    
!!!  !print *,"**********************************************"
!!!  !STOP
!!!#ifdef VISCOUS
!!!  CALL PDEViscEigenvalues(LL,QL,nv)  
!!!  CALL PDEViscEigenvalues(LR,QR,nv)
!!!  amax = 2.0*MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) )/dist 
!!!#else
!!!  amax = 0. 
!!!#endif 
!!!  absA = 0. 
!!!  DO i = 1, nVar
!!!      absA(i,i) = sR*sL/(sR-sL)  - 0.5*amax ! regular HLL diffusion, only on the diagonal 
!!!  ENDDO  
!!!  !
!!!  IF(QR(1).LT.1e-9.OR.QL(1).LT.1e-9) THEN
!!!      deltaL = 0.
!!!  ENDIF
!!!  !
!!!  !deltaL = 0.
!!!  absA = absA - sR*sL/(sR-sL)*MATMUL( RL, MATMUL(deltaL, iRL) )  ! HLLEM anti-diffusion  
!!!  !    
!!!  !PRINT *, "RoeMatrix"
!!!  CALL RoeMatrix(ARoe,QL,QR,nv)
!!!  !PRINT *, "RoeMatrix done!"
!!!  !
!!!  ARoem = -sL*ARoe/(sR-sL)
!!!  ARoep = +sR*ARoe/(sR-sL)
!!!  ! 
!!!  !DR = ARoep
!!!  !DL = ARoem
!!!  !!!!FL(:) = 0.5*( FR(:) + FL(:) ) + MATMUL(absA, QR(:) - QL(:) )    ! purely conservative flux 
!!!  !!!!FR(:) = FL(:) - 0.5*ncp(:)                                                              ! subtract the non-conservative product 
!!!  !!!!FL(:) = FL(:) + 0.5*ncp(:)
!!!  !
!!!  dQ = QR - QL
!!!  fL = (sR*fL - sL*fR)/(sR-sL) + MATMUL( absA, dQ ) 
!!!  !
!!!  !Dp = -MATMUL(Aroep,dQ)
!!!  !Dm = -MATMUL(Aroem,dQ)        ! these are the path integral of the MatrixB from QL to QR. (the NCP as a first approximation)
!!!  !
!!!  Dp = MATMUL(Aroep,dQ)
!!!  Dm = MATMUL(Aroem,dQ)        ! these are the path integral of the MatrixB from QL to QR. (the NCP as a first approximation)
!!!  !
!!!  fR = fL - Dp
!!!  fL = fL + Dm
!!!  ! 
!!!  
!!!  ! if(QavL(21)>1.e-3) then
!!!  !   print *,'---------------------'
!!!  !   print *,QavL
!!!  !   print *, '====================='
!!!  !   print *,QavR
!!!  !   print *, '====================='
!!!  !   print *, fR
!!!  !   print *, '====================='
!!!  !   print *, fL
!!!  !   print *,'---------------------'
!!!  ! end if
!!!!IF(ANY(Dp(6:8).NE.0)) THEN
!!!!     PRINT *,"QR",QR
!!!!     PRINT *,"QL",QL
!!!!     PRINT *,"dQ",dQ
!!!!     PRINT *,"Dp(6:8).NE.0",Dp(6:8)
!!!!    STOP
!!!!ENDIF
!!!!IF(ANY(Dm(6:8).NE.0)) THEN
!!!!     PRINT *,"QR",QR
!!!!     PRINT *,"QL",QL
!!!!     PRINT *,"dQ",dQ
!!!!     PRINT *,"Dm(6:8).NE.0",Dm(6:8)
!!!!    STOP
!!!!ENDIF
!!!!IF(ANY(fR(6:8).NE.0)) THEN
!!!!     PRINT *,"QR",QR
!!!!     PRINT *,"QL",QL
!!!!     PRINT *,"dQ",dQ
!!!!     PRINT *,"fR(6:8).NE.0",fR(6:8)
!!!!    STOP
!!!!ENDIF
!!!!IF(ANY(fl(6:8).NE.0)) THEN
!!!!     PRINT *,"QR",QR
!!!!     PRINT *,"QL",QL
!!!!     PRINT *,"dQ",dQ
!!!!     PRINT *,"fl(6:8).NE.0",fl(6:8)
!!!!    STOP
!!!!ENDIF
!!!!  ! REMEMBER THE FOLLOWING: we are recursively updating qh as
!!!!  ! q_i^{n+1} = q_i^n - FL              .... i.e. F_=F_{i+1/2}_ right flux
!!!!  ! q_{i+1}^{n+1} = q_i^n + FR             .... i.e. FR=F_{i+1/2} left flux
!!!!  ! see musclhancock.cpph after "// 4. Solve Riemann problems"
!!!  ! 
!!!    END SUBROUTINE HLLEMFluxFV

RECURSIVE SUBROUTINE HLLEMFluxFV(FL,FR,QL,QR,NormalNonZero) 
  USE MainVariables, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  ! Local variables
  INTEGER, INTENT(IN)   :: NormalNonZero
  REAL, INTENT(IN)     :: QL(nVar)
  REAL, INTENT(IN)     :: QR(nVar)
  REAL, INTENT(INOUT)  :: FL(nVar)
  REAL, INTENT(INOUT)  :: FR(nVar)
  REAL    :: QavL(nVar), QavR(nVar)  
    ! Local variables 
  INTEGER           :: i,j,k,l, ml(1)  ,iErr
  REAL              :: smax, Qav(nVar)
  REAL              ::  nv(3), flattener(nLin)
  REAL    :: absA(nVar,nVar), amax  ,gradQ(nVar,3), ncp(nVar)
  REAL    :: QM(nVar),LL(nVar),LR(nVar),LM(nVar)
  REAL    :: deltaL(nLin,nLin),Lam(nLin,nLin),Lap(nLin,nLin) 
  REAL    :: RL(nVar,nLin),iRL(nLin,nVar),LLin(nLin,nLin) , TMPM(nLin, nVar),TMPM2(nVar,nVar)
  REAL    :: Aroe(nVar,nVar),Aroep(nVar,nVar), Aroem(nVar,nVar), Dm(nVar), Dp(nVar), dQ(nVar)
  REAL :: f1R(nVar), g1R(nVar), h1R(nVar) ,flux(nVar)
  REAL :: f1L(nVar), g1L(nVar), h1L(nVar) , VM(nVar) 
  
  REAL :: XX0(3),TIME0
  if(any(isnan(QL)).or.any(isnan(QR))) then
    print *, 'QL or OR NaN'
  end if
     
 
  XX0=0.
  TIME0=0.
  !  
  nv(:)=0.
  nv(NormalNonZero+1)=1.
  !
  flattener=0.8 !Duo: original=0.8
  !
CALL PDEFlux(f1L,g1L,h1L,QL)
CALL PDEFlux(f1R,g1R,h1R,QR)
!
fR = f1R*nv(1)+g1R*nv(2)+h1R*nv(3)
fL = f1L*nv(1)+g1L*nv(2)+h1L*nv(3)
  
  QM=0.5*(QL+QR)
  
  CALL PDEEigenvalues(LL,QL,nv,xg)  
  CALL PDEEigenvalues(LR,QR,nv,xg)  
  CALL PDEEigenvalues(LM,QM,nv,xg)  
  sL = MIN( 0., MINVAL(LL(:)), MINVAL(LM(:)) ) 
  sR = MAX( 0., MAXVAL(LR(:)), MAXVAL(LM(:)) )  
  CALL PDEIntermediateFields(RL,LLin,iRL,QM,nv) 
  Lam = 0.5*(LLin-ABS(LLin))
  Lap = 0.5*(LLin+ABS(LLin)) 

  deltaL = 0.0
  DO i = 1, nLin
      deltaL(i,i) = ( 1. - Lam(i,i)/(sL-1e-14) - Lap(i,i)/(sR+1e-14) )*flattener(i)  
  ENDDO  
  amax = 0. 

  absA = 0. 
  DO i = 1, nVar
      absA(i,i) = sR*sL/(sR-sL)  - 0.5*amax ! regular HLL diffusion, only on the diagonal 
  ENDDO  
  TMPM=MATMUL(deltaL, iRL)
  TMPM2=MATMUL( RL,TMPM)
  absA = absA - sR*sL/(sR-sL)*TMPM2 ! HLLEM anti-diffusion  
  !    
  !PRINT *, "RoeMatrix"
  CALL RoeMatrix(ARoe,QL,QR,nv)
  !PRINT *, "RoeMatrix done!"
  !
  ARoem = -sL*ARoe/(sR-sL)
  ARoep = +sR*ARoe/(sR-sL)
  ! 

  dQ = QR - QL
  flux = (sR*fL - sL*fR)/(sR-sL) + MATMUL( absA, QR - QL ) 
  !
  !Dp = -MATMUL(Aroep,dQ)
  !Dm = -MATMUL(Aroem,dQ)        ! these are the path integral of the MatrixB from QL to QR. (the NCP as a first approximation)
  !
  gradQ=0.
  gradQ(:,NormalNonZero+1) = QR(:) - QL(:) 
  CALL PDENCP(ncp,QM,gradQ)
  Dp = MATMUL(Aroep,dQ)
  Dm = MATMUL(Aroem,dQ)        ! these are the path integral of the MatrixB from QL to QR. (the NCP as a first approximation)
  !
  !if(abs(flux(21))>1.e-0) then
  !   print *, '*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
	!	print *,flux
	!	print *,'---------------------'
	!	print *,fl
	!	print *,'---------------------'
	!	print *,fr
	!	print *,'---------------------'
	!	print *,Ql
	!	print *,'---------------------'
	!	print *,qr
  !     print *, '====================='
  !     print *, sR
  !     print *, '====================='
  !     print *, sL
	! print *, '====================='
	!	print *,'---------------------'
	!	print *,deltaL
	!	print *,'---------------------'
	!	print *, MATMUL( absA, QR - QL ) 
	!	print *,'---------------------'
	!	print *, sR*sL/(sR-sL)
	!	print *,'---------------------'
	!	print *, absA(18:nVar,:)
	!	print *,'---------------------'
	!	print *, absA(:,:)
	!	print *, '*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'
	!	pause
  !end if
  !flux = 0.5*( fR + fL ) + MATMUL(absA, QR - QL)
  fR = flux - Dp
  fL = flux + Dm
  

  
 !fR(18:nVar) = 0.- 0.5*ncp(18:nVar)
 !fL(18:nVar) = 0.+ 0.5*ncp(18:nVar) 
  if(any(isnan(fR))) then
	!print *, 'Issue here'
    !print *,'---------------------'
    !print *,QL
    !print *, '====================='
    !print *,QR
    !print *, '====================='
    !print *, fR
    !print *, '====================='
    !print *, fL
    !print *, '====================='
    !print *, sR
    !print *, '====================='
    !print *, sL
	! print *, '====================='
	!print *, deltaL
    !print *, '***********************'
    !print *, Lam
    !print *, '***********************'
    !print *, ncp
	!print *, '***********************'
    !!print *, Lam
	!!    print *, '***********************'
    !!print *, Lam
    !print *,'---------------------'	
	pause
  end if
    END SUBROUTINE HLLEMFluxFV	

RECURSIVE SUBROUTINE HLLEMRiemannSolver(basisSize,NormalNonZero,lFbndLio,lFbndRio,lQbndL,lQbndR,QavL,QavR) 
  USE MainVariables, ONLY : nVar, nDim, nLin
  USE iso_c_binding 
  ! Local variables
  INTEGER, INTENT(IN)   :: NormalNonZero, basisSize
  REAL, INTENT(IN)     :: lQbndL(nVar,basisSize,basisSize)
  REAL, INTENT(IN)     :: lQbndR(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndLio(nVar,basisSize,basisSize)
  REAL, INTENT(INOUT)  :: lFbndRio(nVar,basisSize,basisSize)
  
  REAL				   :: lFbndL(nVar,basisSize,basisSize)
  REAL				   :: lFbndR(nVar,basisSize,basisSize)
    ! Local variables 
	REAL :: f(nVar), g(nVar), h(nVar)
INTEGER           :: i,j,k,l, ml(1)  
REAL              :: aux(nDim), Id(nVar,nVar), smax, Qav(nVar),QavL(nVar), QavR(nVar) 
REAL              ::  xGP, yGP, xdeb, ydeb  
REAL              :: Bn(nVar,nVar), DM(nVar,nVar), ncp(nVar), nv(3)
REAL              :: gradQ(nVar,3), src(nVar),flattener(nLin)
  REAL    :: absA(nVar,nVar), amax  
  REAL    :: QM(nVar),LL(nVar),LR(nVar),LM(nVar)
  REAL    :: deltaL(nLin,nLin),Lam(nLin,nLin),Lap(nLin,nLin) 
  REAL    :: RL(nVar,nLin),iRL(nLin,nVar),LLin(nLin,nLin) 
  REAL :: f1R(nVar), g1R(nVar), h1R(nVar)
  REAL :: f1L(nVar), g1L(nVar), h1L(nVar)
  
  lFbndL=lFbndLio
  lFbndR=lFbndRio
  nv(:)=0.
  nv(NormalNonZero+1)=1.
  !print *, "Normal non zero in fortran=" NormalNonZero
  !print *, "basisSize=", basisSize
  !print *, "NormalNonZero=", NormalNonZero
  !print *, "QavR=",QavR(1)
  !return
  !nv(NormalNonZero)=1.;
  !FL=0.
  !FR=0.
  ! CALL PDEFlux(f1L,g1L,h1L,QL)
  ! CALL PDEFlux(f1R,g1R,h1R,QR)
  !!
  !lFbndR = f1R*nv(1)+g1R*nv(2)+h1R*nv(3)
  !lFbndL = f1L*nv(1)+g1L*nv(2)+h1L*nv(3)
  
	flattener=1.
	!lFbndL=0.
	!lFbndR=0.
    CALL PDEEigenvalues(LL,QavL,nv) 
    CALL PDEEigenvalues(LR,QavR,nv) 
    smax = MAX( MAXVAL(ABS(LL)), MAXVAL(ABS(LR)) )
    ! Identity matrix 
    Id = 0.0 
    DO i=1,nVar
        Id(i,i)=1.0
    ENDDO
    gradQ = 0.0      
    ! HLLEM
    Qav = 0.5*(QavL+QavR) 
    CALL PDEIntermediateFields(RL,LLin,iRL,Qav,nv) 
    Lam = 0.5*(LLin-ABS(LLin))
    Lap = 0.5*(LLin+ABS(LLin)) 
    deltaL = 0.0
    DO i = 1, nLin
        deltaL(i,i) = ( 1. - Lam(i,i)/(-smax-1e-14) - Lap(i,i)/(smax+1e-14) )*flattener(i)  
    ENDDO    
    absA = 0. 
    DO i = 1, nVar
        absA(i,i) = -0.5*smax  ! regular Rusanov diffusion, only on the diagonal 
    ENDDO  
    absA = absA + 0.5*smax*MATMUL( RL, MATMUL(deltaL, iRL) )  ! HLLEM anti-diffusion

        DO k = 1, basisSize
          DO j = 1, basisSize
                Qav = 0.5*(lQbndR(:,j,k)+lQbndL(:,j,k)) 
                gradQ(:,NormalNonZero+1) = lQbndR(:,j,k) - lQbndL(:,j,k) 
                CALL PDENCP(ncp,Qav,gradQ)
				!lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) + MATMUL(absA, lQbndR(:,j,k) - lQbndL(:,j,k) )    ! purely conservative flux 
				lFbndL(:,j,k) = 0.5*( lFbndR(:,j,k) + lFbndL(:,j,k) ) - 0.5*smax*( lQbndR(:,j,k) - lQbndL(:,j,k) )      ! purely conservative flux 
				lFbndR(:,j,k) = lFbndL(:,j,k) - 0.5*ncp(:)                                                              ! subtract the non-conservative product 
                lFbndL(:,j,k) = lFbndL(:,j,k) + 0.5*ncp(:)
            ENDDO
        ENDDO			
	lFbndLio=	lFbndL
    lFbndRio=	lFbndR
	
  !if(QavL(21)>1.e-3) then
	!  print *,'---------------------'
	!  print *,QavL
	!  print *, '====================='
	!  print *,QavR
	!  print *, '====================='
	!  print *, lFbndR(:,1,1)
	!  print *, '====================='
	!  print *, lFbndL(:,1,1)
	!  print *,'---------------------'
  !end if
END SUBROUTINE HLLEMRiemannSolver

!RECURSIVE SUBROUTINE InitTECPLOT(N_in,M_in)
!	USE TECPLOTPLOTTERmod
!	implicit none
!	INTEGER :: N_in,M_in
!	CALL SetMainParameters(N_in,M_in)
!END SUBROUTINE InitTECPLOT
!
!RECURSIVE SUBROUTINE getNumericalSolution(V,Q) 
!  USE MainVariables, ONLY: nVar  
!  IMPLICIT NONE     
!  REAL				:: V(nVar), Q(nVar)
!  CALL PDECons2Prim(V,Q)
!END SUBROUTINE getNumericalSolution
!
!RECURSIVE SUBROUTINE getExactSolution(V,pos, timeStamp) 
!  USE MainVariables, ONLY: nVar , nDim  
!  IMPLICIT NONE     
!  REAL				:: V(nVar), Q(nVar), pos(nDim), timeStamp
!  call InitialData(pos, timeStamp, Q)
!  CALL PDECons2Prim(V,Q)
!  !V(1)=pos(1)**2!*pos(2)**2
!END SUBROUTINE getExactSolution


#endif
