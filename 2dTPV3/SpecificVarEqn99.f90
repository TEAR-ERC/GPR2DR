#ifndef __HEADER_SPECIFICVAREN__
#define __HEADER_SPECIFICVAREN__
!#include "expintegrator_type.f90"
!#include "MainVariables.f90"
!#include "expintegrator_linalg.f90"

    module SpecificVarEqn99

        IMPLICIT NONE
        private
        public :: ComputeLMc
        public :: GPRsigma2ASGEOS
        public :: compute_l_m_dedc
        public :: compute_hydrodynamic_pressure
        public :: compute_total_stress
        public :: compute_coulomb
        public :: compute_temperature
        public :: compute_s_of_t
        public :: StaticLimiterEQ99
        public :: SmoothInterface
        public :: compute_l_m_mix
        public :: computeGPRLEstressGen 
    contains

subroutine computeGPRLEstressGen(stressnorm,sigma_vec,Avec, LL,MM,rho,ss, addp)
    implicit none
    real, intent(in) :: Avec(9),LL,MM,rho,ss        ! For a general EOS we need the A components, LL,MM, density rho and entropy ss
    real, intent(out) :: sigma_vec(6)
    real :: detA,stressnorm,pressure
    logical, intent(in) :: addp
    real :: AM(3,3), G(3,3),Id(3,3),devG(3,3),TT(3,3), GT(3,3),YY,YY_sh,YY_p
    
    AM(1,:) = (/ Avec(1), Avec(2), Avec(3) /) 
    AM(2,:) = (/ Avec(4), Avec(5), Avec(6) /)
    AM(3,:) = (/ Avec(7), Avec(8), Avec(9)/) 
    
    G      = MATMUL( TRANSPOSE(AM), AM )
    Id     = 0. 
    Id(1,1) = 1.0; Id(2,2) = 1.0; Id(3,3) = 1.0 
    GT     = MATMUL( TRANSPOSE(AM), AM )
    devG   = GT - (GT(1,1)+GT(2,2)+GT(3,3))/3.*Id
    TT     = -MM*MATMUL(GT,devG) 
    sigma_vec(1)=TT(1,1)
    sigma_vec(2)=TT(2,2)
    sigma_vec(3)=TT(3,3)
    sigma_vec(4)=TT(1,2)
    sigma_vec(5)=TT(2,3)
    sigma_vec(6)=TT(1,3)

    if(addp) then ! Add the pressure contribution diff_k p delta_ik
        !pressure=(LL+2.0/3.0*MM)*(detA)**2*(1-detA)
        detA = Avec(1) * Avec(5) * Avec(9) - Avec(1) * Avec(6) * Avec(8) - Avec(2) * Avec(4) * Avec(9) + Avec(2) * Avec(6) * Avec(7) + Avec(3) * Avec(4) * Avec(8) - Avec(3) * Avec(5) * Avec(7)
                pressure=-(LL+2.0/3.0*MM)*(detA)**2*(1-detA)
        !call computeGPRLEpressure(pressure,Avec, LL,MM,rho,ss)
        sigma_vec(1)=sigma_vec(1)-pressure       ! A
        sigma_vec(2)=sigma_vec(2)-pressure
        sigma_vec(3)=sigma_vec(3)-pressure 
    end if

!    stressnorm  = SQRT( 0.5 * ( (sigma_vec(1)-sigma_vec(2))**2 + (sigma_vec(2)-sigma_vec(3))**2 + (sigma_vec(3)-sigma_vec(1))**2 + 6.*(sigma_vec(4)**2+sigma_vec(6)**2+sigma_vec(5)**2) ) )  
    ! stressnorm = 0.25*abs(sigma_vec(1)-sigma_vec(2))+0.8*abs(sigma_vec(1)+2.0*sigma_vec(2))/3.0
    !YY_sh=SQRT( 0.5 * ( (TT(1,1)-TT(2,2))**2 + (TT(2,2)-TT(3,3))**2 + (TT(3,3)-TT(1,1))**2 + 6.*(TT(1,2)**2+TT(1,3)**2+TT(2,3)**2) ) )  
    !YY_p=abs(pressure)
    !YY=0.8*YY_p
    !!if(YY_sh .ge. 2.0e+9) then
    !    YY=YY+0.25*YY_sh*(0.5+0.5*erf(  (YY_sh-2.e+9)/5.e+8   ))
    !!end if
    !stressnorm=YY
    !call computeEquivalentStress(stressnorm,TT,pressure)
end subroutine computeGPRLEstressGen    
    
    ! --------------------------------------------------------------------------------------------------------------------
    pure subroutine ComputeLMc(LL,MM,dMM,dKK,mu1,mu2,lam1,lam2,c)
        ! Compute LL,MM ,dMM/dc, dKK/dc according to the damage parameter c
        implicit none
        real, intent(out) :: LL, MM,dMM,dKK
        real, intent(in)  ::  lam1,lam2,mu1,mu2,c
        real :: mu_tilda,K1,K2,KK_tilda

        mu_tilda=c*mu1+(1-c)*mu2
        K1=lam1+2.0/3.0*mu1
        K2=lam2+2.0/3.0*mu2
        KK_tilda=c*K1+(1-c)*K2
    
        MM=mu1*mu2/mu_tilda
        LL=K1*K2/KK_tilda-2.0/3.0*MM

        dMM=(mu1-mu2)*mu1*mu2/mu_tilda**2
        dKK=(K1-K2)*K1*K2/KK_tilda**2

    end subroutine ComputeLMc
    ! --------------------------------------------------------------------------------------------------------------------
    function GPRsigma2ASGEOS(sigmaK,theta,lambda1,mu1,rho0,gamma,s,cv,p0,eps,EOS_mode,T0)
        implicit none
        real, intent(in)    :: sigmaK(6),theta,lambda1,mu1,rho0, T0
        integer, intent(in) :: EOS_mode
        real                :: GPRsigma2ASGEOS(9),Aloc(9),yy
        ! Local variables
        real                :: eps,unkvarsig(6),fdf(6),fp(6),JJ(6,6),ff(6),fpdeb(6),gamma,s,cv,p0
        integer             :: iNewton,maxNewton,ierr
        
        maxNewton=50
        unkvarsig=0.1
        unkvarsig(1)=1
        unkvarsig(3)=1
        unkvarsig(6)=1
        !unkvarsig=sigmaK

        do iNewton=1,maxNewton
            call A2sigmaJacobianSGEOS(JJ,unkvarsig,theta,lambda1,mu1,rho0,gamma,s,cv,p0,EOS_mode, T0)
            
            call A2sigmaComputeFSGEOS(fp,unkvarsig,theta,lambda1,mu1,rho0,gamma,s,cv,p0,EOS_mode, T0)
            
            ff=fp-sigmaK
            call LinSolve(6,JJ,ff,fdf)
            unkvarsig=unkvarsig-fdf
            if(sqrt(abs(sum(ff**2)))<eps) then
                goto 146    
            end if
        end do
        print *, ff
        print *, sum(ff(1:6))
        print *, 'Newton in sigma2A does not converge, err=',sqrt(abs(sum(ff(1:6)**2))),">eps=",eps
		print *, 'Failure:', lambda1,mu1,rho0,s,cv,p0,EOS_mode,T0
146     continue
        call A2sigmaAssignQA(GPRsigma2ASGEOS,unkvarsig, theta)
    end function GPRsigma2ASGEOS
    ! --------------------------------------------------------------------------------------------------------------------
    subroutine A2sigmaAssignQA(QA,sig, theta)
        implicit none
        real,intent(out) :: QA(9)
        real             :: Aij(3,3)
        real, intent(in) :: sig(6), theta
        integer          :: k,i,j
        Aij=0.
        Aij(1,1) = sig(1) * cos(theta)
        Aij(1,2) = -sig(1) * sin(theta)
        Aij(2,1) = sig(2) * cos(theta) + sig(3) * sin(theta)
        Aij(2,2) = -sig(2) * sin(theta) + sig(3) * cos(theta)
        Aij(3,1) = sig(4) * cos(theta) + sig(5) * sin(theta)
        Aij(3,2) = -sig(4) * sin(theta) + sig(5) * cos(theta)
        Aij(3,3) = sig(6)
        k=0
        do i=1,3
            do j=1,3
                k=k+1
                QA(k)=Aij(i,j)
            end do
        end do
        
        
    end subroutine A2sigmaAssignQA
    ! --------------------------------------------------------------------------------------------------------------------
subroutine A2sigmaJacobianSGEOS(J,sig,theta,lambda1,mu1,rho0,g,s,cv,p0,EOS_mode, T0)
        implicit none
        real :: s1,s2,s3,s4
        real,intent(out)    :: J(6,6)
        real, intent(in)    :: sig(6),theta,lambda1,mu1,rho0,g,s,cv,p0, T0
        integer, intent(in) :: EOS_mode
        
        select case(EOS_mode)
        case(1)
    J(1,1) = -rho0 * (cos(theta) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.2D&
     &1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(th&
     &eta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) *&
     &* 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / &
     &0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * cos(theta) - 0.4D1 * (-sig&
     &(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta))) * sig(1) * sin(theta) - 0.4D1 / 0.3D1 * (0&
     &.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * c&
     &os(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2&
     & / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * cos(theta) - 0.4D1 / 0.&
     &3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2&
     & / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) ** 2 / 0.3D1) * sig(1) * cos(theta)) / 0.4D1 + sig(1) * &
     &cos(theta) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1)&
     & * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2) * s&
     &ig(1) * cos(theta) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 &
     &* cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) &
     &* sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.&
     &3D1) * cos(theta) + 0.8D1 * sig(1) ** 2 * cos(theta) * sin(theta) &
     &** 2 - 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) *&
     & cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-&
     &sig(4) * sin(theta) + sig(5) * cos(theta))) * sin(theta) - 0.4D1 /&
     & 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * sin(theta) ** 2 - 0.2D1 / 0.3D1&
     & * sig(1) * cos(theta) ** 2) * sig(1) * cos(theta) - 0.4D1 / 0.3D1&
     & * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 &
     &* (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D&
     &1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** &
     &2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * cos(theta) - 0.4D1 / 0.3D1 &
     &* (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig&
     &(1) * sin(theta) ** 2) * sig(1) * cos(theta) - 0.4D1 / 0.3D1 * (0.&
     &2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 &
     &- sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & ** 2 / 0.3D1) * cos(theta)) / 0.4D1 + (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * &
     &sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** &
     &2) * (sig(2) * cos(theta) + sig(3) * sin(theta)) - 0.8D1 * sig(1) &
     &* cos(theta) * sin(theta) * (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) - 0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * sin(theta) ** 2&
     & - 0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2) * (sig(2) * cos(theta&
     &) + sig(3) * sin(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig(1&
     &) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2) * &
     &(sig(2) * cos(theta) + sig(3) * sin(theta))) / 0.4D1 + (sig(4) * c&
     &os(theta) + sig(5) * sin(theta)) * mu1 / rho0 * (0.8D1 / 0.3D1 * (&
     &0.4D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) &
     &* sin(theta) ** 2) * (sig(4) * cos(theta) + sig(5) * sin(theta)) -&
     & 0.8D1 * sig(1) * cos(theta) * sin(theta) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1) *&
     & sin(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2) * (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) - 0.4D1 / 0.3D1 * (-0.2D1&
     & / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin&
     &(theta) ** 2) * (sig(4) * cos(theta) + sig(5) * sin(theta))) / 0.4&
     &D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(theta&
     &) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) ** g * (&
     &cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6)&
     & + sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(&
     &6)) / (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sig(6)) * exp(s / cv)
      J(1,2) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.8D1 / 0.3D&
     &1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) * sin(theta)) * sig(1) * cos(theta) - 0.4D1 * (cos(theta)&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - (sig(2) * cos(th&
     &eta) + sig(3) * sin(theta)) * sin(theta)) * sig(1) * sin(theta) - &
     &0.4D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * cos(theta)) * sig(1) * cos(theta) - 0.4D1 &
     &/ 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) * sin(theta)) * sig(1) * cos(theta)) / 0.4D1 + cos&
     &(theta) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) **&
     & 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 /&
     & 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.4D1 * (-&
     &sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) +&
     & (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta))) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta&
     &) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) - 0.4D1 / 0.3D1 * (0.2D1 /&
     & 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig&
     &(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2&
     & / 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(theta))) / 0.4D1 +&
     & (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (0.8D1&
     & / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) * sin(theta)) * (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(the&
     &ta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * co&
     &s(theta) + 0.4D1 * (cos(theta) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) - (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(th&
     &eta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 * (-s&
     &ig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta))) * sin(theta) - 0.4D1 / 0.3D1 * (-0.4D1 /&
     & 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta)&
     & - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * c&
     &os(theta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) - 0.4D1 /&
     & 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / &
     &0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 &
     &/ 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(&
     &1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) &
     &* sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * cos(theta) - 0.4D1 / &
     &0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) * sin(theta)) * (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) **&
     & 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3&
     &D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) &
     &* cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * cos(theta)) / 0.4D1 &
     &+ (sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1 / rho0 * (0.8D&
     &1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) * sin(theta)) * (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) + 0.4D1 * (cos(theta) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) - (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin(&
     &theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.4D1 / 0&
     &.3D1 * (-0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) *&
     & sin(theta)) * cos(theta)) * (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) * sin(theta)) * (sig(4) * cos(theta&
     &) + sig(5) * sin(theta))) / 0.4D1)
      J(1,3) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.8D1 / 0.3D&
     &1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta))) * sig(1) * cos(theta) - 0.4D1 * ((-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) * sin(theta) + (sig(2) * cos(th&
     &eta) + sig(3) * sin(theta)) * cos(theta)) * sig(1) * sin(theta) - &
     &0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) * sin(theta)) * sig(1) * cos(theta) - 0.4D1 /&
     & 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta))) * sig(1) * cos(theta)) / 0.4D1 + sin(&
     &theta) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** &
     &2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / &
     &0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.4D1 * (-s&
     &ig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta))) * (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta)&
     & ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(the&
     &ta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / &
     &0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(&
     &1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) &
     &* sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 &
     &/ 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin(theta))) / 0.4D1 + &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (0.8D1 &
     &/ 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta))) * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(thet&
     &a) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sin&
     &(theta) + 0.4D1 * ((-sig(2) * sin(theta) + sig(3) * cos(theta)) * &
     &sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(the&
     &ta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.4D1 * (-si&
     &g(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta)&
     & + sig(5) * cos(theta))) * cos(theta) - 0.4D1 / 0.3D1 * (0.4D1 / 0&
     &.3D1 * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) -&
     & 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sin&
     &(theta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) - 0.4D1 / 0&
     &.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.&
     &3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / &
     &0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1)&
     & ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sin(theta) - 0.4D1 / 0.&
     &3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta))) * (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2&
     & * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1&
     & - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * &
     &cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * sin(theta)) / 0.4D1 + &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1 / rho0 * (0.8D1 &
     &/ 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta))) * (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) + 0.4D1 * ((-sig(2) * sin(theta) + sig(3) * cos(theta)) *&
     & sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(th&
     &eta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.4D1 / 0.3&
     &D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) * sin(theta)) * (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta))) * (sig(4) * cos(theta) +&
     & sig(5) * sin(theta))) / 0.4D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) * (si&
     &g(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) *&
     & sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) * sig(6)) ** g * (sig(1) * cos(theta) ** 2 * sig(6) + si&
     &g(1) * sin(theta) ** 2 * sig(6)) / (sig(1) * cos(theta) * (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(theta&
     &) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) * exp(s &
     &/ cv)
      J(1,4) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.8D1 / 0.3D&
     &1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos&
     &(theta)) * sin(theta)) * sig(1) * cos(theta) - 0.4D1 * (cos(theta)&
     & * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) * sin(theta)) * sig(1) * sin(theta) - &
     &0.4D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) * cos(theta)) * sig(1) * cos(theta) - 0.4D1 &
     &/ 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) * sin(theta)) * sig(1) * cos(theta)) / 0.4D1 + (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (0.8D1 / 0&
     &.3D1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * &
     &cos(theta)) * sin(theta)) * (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) + 0.4D1 * (cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos&
     &(theta)) - (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta&
     &)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 &
     &* (-0.4D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) *&
     & sin(theta) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(&
     &theta)) * cos(theta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) * sin(theta)) * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta))) / 0.4D1 + cos(theta) * mu1 / rho0 * (0.8D1 / &
     &0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0&
     &.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / &
     &0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) &
     &** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(t&
     &heta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) + 0.4D1 * (sig(4) * cos(&
     &theta) + sig(5) * sin(theta)) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2&
     &D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig&
     &(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos&
     &(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 /&
     & 0.3D1 - sig(6) ** 2 / 0.3D1) * (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) &
     &** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0&
     &.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * (sig(4) * cos(thet&
     &a) + sig(5) * sin(theta))) / 0.4D1 + (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) + 0.2D1 / 0&
     &.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta)) &
     &* (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.8D1 / 0.3D1 * (0&
     &.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin&
     &(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2&
     & / 0.3D1 - sig(6) ** 2 / 0.3D1) * cos(theta) + 0.4D1 * (cos(theta)&
     & * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) * sin(theta)) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) - 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(&
     &theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) * si&
     &n(theta) + 0.4D1 * cos(theta) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (-0.&
     &4D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(t&
     &heta) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) * cos(theta)) * (sig(4) * cos(theta) + sig(5) * sin(theta)) - 0.&
     &4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2&
     &D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0&
     &.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 -&
     & sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * cos(theta) - 0.4&
     &D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) * sin(theta)) * (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(&
     &1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 &
     &/ 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1&
     & - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * cos(theta)) / 0&
     &.4D1)
      J(1,5) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.8D1 / 0.3D&
     &1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta))) * sig(1) * cos(theta) - 0.4D1 * ((-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) * cos(theta)) * sig(1) * sin(theta) - &
     &0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) * sin(theta)) * sig(1) * cos(theta) - 0.4D1 /&
     & 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta))) * sig(1) * cos(theta)) / 0.4D1 + (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (0.8D1 / 0.&
     &3D1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta)&
     & + sig(5) * cos(theta))) * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) + 0.4D1 * ((-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin&
     &(theta) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta)&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 *&
     & (0.4D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) * sin(theta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &- 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta))) * (sig(2) * cos(theta) + sig&
     &(3) * sin(theta))) / 0.4D1 + sin(theta) * mu1 / rho0 * (0.8D1 / 0.&
     &3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3&
     &D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.&
     &3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) **&
     & 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(the&
     &ta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) * (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) + 0.4D1 * (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1&
     & / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos(t&
     &heta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     &* 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0&
     &.3D1 - sig(6) ** 2 / 0.3D1) * (sig(4) * cos(theta) + sig(5) * sin(&
     &theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) **&
     & 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3&
     &D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) &
     &* cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * (sig(4) * cos(theta)&
     & + sig(5) * sin(theta))) / 0.4D1 + (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) - 0.2D1 / 0.3&
     &D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) * &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.8D1 / 0.3D1 * (0.2&
     &D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(t&
     &heta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 /&
     & 0.3D1 - sig(6) ** 2 / 0.3D1) * sin(theta) + 0.4D1 * ((-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(4) * cos(thet&
     &a) + sig(5) * sin(theta)) * cos(theta)) * (-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(th&
     &eta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) * cos(&
     &theta) + 0.4D1 * sin(theta) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.4D1&
     & / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & sin(theta)) * (sig(4) * cos(theta) + sig(5) * sin(theta)) - 0.4D1&
     & / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 &
     &/ 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D&
     &1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - si&
     &g(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sin(theta) - 0.4D1 &
     &/ 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta))) * (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) &
     &** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0&
     &.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * sin(theta)) / 0.4D&
     &1)
      J(1,6) = -rho0 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 *&
     & mu1 / rho0 * sig(6) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 * mu1 / rho0 * sig(6) + 0.4D1 / 0.3D1 * (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) ** 2 * mu1 / rho0 * sig(6)&
     &) - (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(theta) &
     &* (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) ** g * (si&
     &g(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) +&
     & sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta))&
     &) / (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(th&
     &eta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * sig(6)) * exp(s / cv)
      J(2,1) = -rho0 * (-sin(theta) * mu1 / rho0 * (0.4D1 / 0.3D1 * (0.2&
     &D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(t&
     &heta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 /&
     & 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * sin(theta) + 0.4D1 * (-si&
     &g(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta)&
     & + sig(5) * cos(theta))) * sig(1) * cos(theta) - 0.8D1 / 0.3D1 * (&
     &0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * &
     &cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** &
     &2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0&
     &.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** &
     &2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta)&
     & + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) ** 2 / 0.3D1) * sig(1) * sin(theta)) / 0.4D1 - sig(1) *&
     & sin(theta) * mu1 / rho0 * (0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1&
     &) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2) * &
     &sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2&
     & * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0&
     &.3D1) * sin(theta) - 0.8D1 * sig(1) ** 2 * cos(theta) ** 2 * sin(t&
     &heta) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) &
     &* cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta))) * cos(theta) - 0.8D1 &
     &/ 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * sin(theta) ** 2 - 0.2D1 / 0.3D&
     &1 * sig(1) * cos(theta) ** 2) * sig(1) * sin(theta) - 0.8D1 / 0.3D&
     &1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3&
     &D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) **&
     & 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sin(theta) + 0.4D1 / 0.3D1&
     & * (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * si&
     &g(1) * sin(theta) ** 2) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (0&
     &.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1&
     & - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(&
     &theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) ** 2 / 0.3D1) * sin(theta)) / 0.4D1 + (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1&
     & * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) &
     &** 2) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.8D1 * sig&
     &(1) * cos(theta) * sin(theta) * (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) + 0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * sin(theta) *&
     &* 2 - 0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2) * (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * s&
     &ig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) / 0.4D1 + (-sig(&
     &4) * sin(theta) + sig(5) * cos(theta)) * mu1 / rho0 * (-0.4D1 / 0.&
     &3D1 * (0.4D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * &
     &sig(1) * sin(theta) ** 2) * (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) - 0.8D1 * sig(1) * cos(theta) * sin(theta) * (sig(4) * cos(&
     &theta) + sig(5) * sin(theta)) + 0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * s&
     &ig(1) * sin(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2&
     &) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.4D1 / 0.3D1 *&
     & (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(&
     &1) * sin(theta) ** 2) * (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &))) / 0.4D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * &
     &sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6))&
     & ** g * (cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & * sig(6) + sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * sig(6)) / (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(th&
     &eta) + sig(3) * sin(theta)) * sig(6)) * exp(s / cv)
      J(2,2) = -rho0 * (-sig(1) * sin(theta) * mu1 / rho0 * (0.4D1 / 0.3&
     &D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) * sin(theta)) * sig(1) * sin(theta) + 0.4D1 * (cos(theta&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * sin(theta)) * sig(1) * cos(theta) -&
     & 0.8D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) * cos(theta)) * sig(1) * sin(theta) + 0.4D1&
     & / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) * sin(theta)) * sig(1) * sin(theta)) / 0.4D1 - si&
     &n(theta) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) &
     &** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 -&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2&
     & / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.4D1 *&
     & (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta))) * (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(the&
     &ta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(th&
     &eta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * c&
     &os(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(thet&
     &a) + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D&
     &1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - &
     &sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) *&
     &* 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) / 0.4&
     &D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 / rho0 * (&
     &-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) * sin(theta)) * (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * &
     &cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D&
     &1) * sin(theta) + 0.4D1 * (cos(theta) * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) - (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sin(theta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.4D1&
     & * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta))) * cos(theta) + 0.8D1 / 0.3D1 * (-0&
     &.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(&
     &theta) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * cos(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - &
     &0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0&
     &.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 +&
     & 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2&
     & - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sin(theta) - 0&
     &.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * s&
     &in(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + &
     &sig(3) * cos(theta)) * sin(theta)) * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - s&
     &ig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) **&
     & 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.&
     &3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * sin(theta)) &
     &/ 0.4D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * mu1 / rho&
     &0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(the&
     &ta) + sig(3) * cos(theta)) * sin(theta)) * (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta)) + 0.4D1 * (cos(theta) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) - (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * sin(theta)) * (sig(4) * cos(theta) + sig(5) * sin(theta)) +&
     & 0.8D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) * cos(theta)) * (-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta)) * (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta))) / 0.4D1)
      J(2,3) = -rho0 * (-sig(1) * sin(theta) * mu1 / rho0 * (0.4D1 / 0.3&
     &D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta))) * sig(1) * sin(theta) + 0.4D1 * ((-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) * sin(theta) + (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * cos(theta)) * sig(1) * cos(theta) -&
     & 0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(thet&
     &a) + sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * sin(theta)) * sig(1) * sin(theta) + 0.4D1 &
     &/ 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta))) * sig(1) * sin(theta)) / 0.4D1 + cos&
     &(theta) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) *&
     &* 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + si&
     &g(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(&
     &4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 &
     &/ 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.4D1 * &
     &(-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta))) * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(thet&
     &a) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta&
     &) + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1&
     & / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - s&
     &ig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) **&
     & 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) / 0.4D&
     &1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 / rho0 * (-&
     &0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * s&
     &in(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta))) * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * c&
     &os(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * s&
     &in(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1&
     &) * cos(theta) + 0.4D1 * ((-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * &
     &cos(theta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.4D1 &
     &* (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta))) * sin(theta) + 0.8D1 / 0.3D1 * (0.4&
     &D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.&
     &8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2&
     &D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0&
     &.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 -&
     & sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * cos(theta) - 0.4&
     &D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta))) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig&
     &(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2&
     & / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D&
     &1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * cos(theta)) / &
     &0.4D1 + (-sig(4) * sin(theta) + sig(5) * cos(theta)) * mu1 / rho0 &
     &* (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta))) * (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) + 0.4D1 * ((-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * cos(theta)) * (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0&
     &.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta)&
     & + sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) * sin(theta)) * (-sig(4) * sin(theta) + sig(5)&
     & * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(th&
     &eta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta))) / 0.4D1) - (lambda1 + 0.2D1 / 0.3D1 &
     &* mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * sig(6)) ** g * (sig(1) * cos(theta) ** 2 *&
     & sig(6) + sig(1) * sin(theta) ** 2 * sig(6)) / (sig(1) * cos(theta&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1)&
     & * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(&
     &6)) * exp(s / cv)
      J(2,4) = -rho0 * (-sig(1) * sin(theta) * mu1 / rho0 * (0.4D1 / 0.3&
     &D1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &* cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) * sin(theta)) * sig(1) * sin(theta) + 0.4D1 * (cos(theta&
     &) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) * sin(theta)) * sig(1) * cos(theta) -&
     & 0.8D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) * cos(theta)) * sig(1) * sin(theta) + 0.4D1&
     & / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(&
     &5) * cos(theta)) * sin(theta)) * sig(1) * sin(theta)) / 0.4D1 + (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 / rho0 * (-0.4D1 &
     &/ 0.3D1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5)&
     & * cos(theta)) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) + 0.4D1 * (cos(theta) * (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) - (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(t&
     &heta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.8D1 / 0.3&
     &D1 * (-0.4D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) * sin(theta) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) * cos(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(th&
     &eta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) * sin(theta)) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta))) / 0.4D1 - sin(theta) * mu1 / rho0 * (-0.&
     &4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2&
     &D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.&
     &2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - s&
     &ig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig&
     &(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * &
     &cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) &
     &* sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig&
     &(2) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     &) * (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.8D1 / 0.3D1 * &
     &(0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 *&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 *&
     & cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) **&
     & 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(4) * sin(theta) + sig(5)&
     & * cos(theta)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 -&
     & sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) &
     &** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / &
     &0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 &
     &- (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta))) / 0.4D1 + (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 /&
     & 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) &
     &+ 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * s&
     &in(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.4D1 &
     &/ 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 /&
     & 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 &
     &/ 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1&
     &) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sin(theta) + 0.4D1 *&
     & (cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta)) * (sig(4) * &
     &cos(theta) + sig(5) * sin(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(th&
     &eta) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta&
     &) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta))) * cos(theta) + 0.8D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta)) * (-sig(&
     &4) * sin(theta) + sig(5) * cos(theta)) - 0.8D1 / 0.3D1 * (0.2D1 / &
     &0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos(thet&
     &a) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2&
     & / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D&
     &1 - sig(6) ** 2 / 0.3D1) * sin(theta) - 0.4D1 * sin(theta) * sig(6&
     &) ** 2 - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) * sin(theta)) * (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6&
     &) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin&
     &(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 &
     &/ 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1&
     & - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * si&
     &n(theta)) / 0.4D1)
      J(2,5) = -rho0 * (-sig(1) * sin(theta) * mu1 / rho0 * (0.4D1 / 0.3&
     &D1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &* sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta))) * sig(1) * sin(theta) + 0.4D1 * ((-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) * cos(theta)) * sig(1) * cos(theta) -&
     & 0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(thet&
     &a) + sig(5) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) * sin(theta)) * sig(1) * sin(theta) + 0.4D1 &
     &/ 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta))) * sig(1) * sin(theta)) / 0.4D1 + (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 / rho0 * (-0.4D1 /&
     & 0.3D1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta))) * (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) + 0.4D1 * ((-sig(4) * sin(theta) + sig(5) * cos(theta)) *&
     & sin(theta) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(th&
     &eta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.8D1 / 0.3D&
     &1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta))) * (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta))) / 0.4D1 + cos(theta) * mu1 / rho0 * (-0.4D&
     &1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1&
     & / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D&
     &1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig&
     &(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * &
     &sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) &
     &* (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.8D1 / 0.3D1 * (0&
     &.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * c&
     &os(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2&
     & / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - s&
     &ig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) **&
     & 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.&
     &3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta))) / 0.4D1 + (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0&
     &.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) - &
     &0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta))) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.4D1 / &
     &0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0&
     &.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / &
     &0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) &
     &** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * cos(theta) + 0.4D1 * (&
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * cos(theta)) * (sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(thet&
     &a) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a))) * sin(theta) + 0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.2D1 / 0.3D1 * (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta)) * (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3&
     &D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos(theta) &
     &** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / &
     &0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 -&
     & sig(6) ** 2 / 0.3D1) * cos(theta) + 0.4D1 * cos(theta) * sig(6) *&
     &* 2 - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta))) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) *&
     &* 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(th&
     &eta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) **&
     & 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0&
     &.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * cos(t&
     &heta)) / 0.4D1)
      J(2,6) = -rho0 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 *&
     & mu1 / rho0 * sig(6) - 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig&
     &(3) * cos(theta)) ** 2 * mu1 / rho0 * sig(6) + 0.4D1 / 0.3D1 * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * mu1 / rho0 * sig(&
     &6)) - (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta) * (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * sin(theta&
     &) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) ** g * (&
     &sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta&
     &))) / (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sig(6)) * exp(s / cv)
      J(3,1) = -sig(6) * mu1 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1)&
     & * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2) * s&
     &ig(6) - 0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * sin(theta) ** 2 &
     &- 0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2) * sig(6) + 0.8D1 / 0.3&
     &D1 * (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * &
     &sig(1) * sin(theta) ** 2) * sig(6)) / 0.4D1 - (lambda1 + 0.2D1 / 0.3D1 &
     &* mu1) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * sig(6)) ** g * (cos(theta) * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) * sig(6) + sin(theta) * (sig(2) * &
     &cos(theta) + sig(3) * sin(theta)) * sig(6)) / (sig(1) * cos(theta)&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) &
     &* sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6&
     &)) * exp(s / cv)
      J(3,2) = -sig(6) * mu1 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta)) * si&
     &g(6) - 0.4D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * cos(theta)) * sig(6) + 0.8D1 / 0.3D&
     &1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) * sin(theta)) * sig(6)) / 0.4D1
      J(3,3) = -sig(6) * mu1 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1&
     & * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) * si&
     &g(6) - 0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(2) * cos(th&
     &eta) + sig(3) * sin(theta)) * sin(theta)) * sig(6) + 0.8D1 / 0.3D1&
     & * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta))) * sig(6)) / 0.4D1 - (lambda1 + 0.2D1 / 0.3D1 * m&
     &u1) * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * sig(6)) ** g * (sig(1) * cos(theta) ** 2 * si&
     &g(6) + sig(1) * sin(theta) ** 2 * sig(6)) / (sig(1) * cos(theta) *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * &
     &sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6))&
     & * exp(s / cv)
      J(3,4) = -sig(6) * mu1 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1&
     & * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta)) * si&
     &g(6) + 0.8D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6&
     &) * cos(theta) - 0.4D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * cos(theta)) * sig(6) - 0.&
     &8D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) * sin(&
     &theta) + 0.8D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) * sin(theta)) * sig(6)) / 0.4D1
      J(3,5) = -sig(6) * mu1 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1&
     & * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) * si&
     &g(6) + 0.8D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6&
     &) * sin(theta) - 0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * sin(theta)) * sig(6) + 0.8&
     &D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) * cos(t&
     &heta) + 0.8D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-&
     &sig(4) * sin(theta) + sig(5) * cos(theta))) * sig(6)) / 0.4D1
      J(3,6) = -mu1 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * c&
     &os(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * s&
     &in(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1&
     &) * sig(6) + 0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     &* 2 * sig(6) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(&
     &theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) &
     &* cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) *&
     & sig(6) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** &
     &2 * sig(6) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1)&
     & ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / &
     &0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 -&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * sig(6)) / 0.4D1 -&
     & sig(6) * mu1 * (0.8D1 * sig(6) ** 2 - 0.4D1 / 0.3D1 * sig(1) ** 2&
     & * cos(theta) ** 2 - 0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) ** 2 + 0.8D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 - 0.4D1 / 0.3D1 * sig(1) ** 2 * sin(theta) **&
     & 2 - 0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &** 2 + 0.8D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) ** 2) / 0.4D1 - (lambda1 + 0.2D1 / 0.3D1 * mu1) * (sig(1) * cos(theta&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1)&
     & * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(&
     &6)) ** g * (sig(1) * cos(theta) * (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) + sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3)&
     & * sin(theta))) / (sig(1) * cos(theta) * (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) * sig(6) + sig(1) * sin(theta) * (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) * sig(6)) * exp(s / cv)
      J(4,1) = -rho0 * (cos(theta) * mu1 / rho0 * (0.4D1 / 0.3D1 * (0.2D&
     &1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(th&
     &eta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) *&
     &* 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / &
     &0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * sin(theta) + 0.4D1 * (-sig&
     &(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + sig(3)&
     & * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta))) * sig(1) * cos(theta) - 0.8D1 / 0.3D1 * (0&
     &.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * c&
     &os(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2&
     & / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * sin(theta) + 0.4D1 / 0.&
     &3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2&
     & / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) ** 2 / 0.3D1) * sig(1) * sin(theta)) / 0.4D1 + sig(1) * &
     &cos(theta) * mu1 / rho0 * (0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1)&
     & * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2) * s&
     &ig(1) * sin(theta) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 &
     &* cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) &
     &* sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.&
     &3D1) * sin(theta) - 0.8D1 * sig(1) ** 2 * cos(theta) ** 2 * sin(th&
     &eta) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) *&
     & cos(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-&
     &sig(4) * sin(theta) + sig(5) * cos(theta))) * cos(theta) - 0.8D1 /&
     & 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * sin(theta) ** 2 - 0.2D1 / 0.3D1&
     & * sig(1) * cos(theta) ** 2) * sig(1) * sin(theta) - 0.8D1 / 0.3D1&
     & * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 &
     &* (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D&
     &1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** &
     &2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sin(theta) + 0.4D1 / 0.3D1 &
     &* (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig&
     &(1) * sin(theta) ** 2) * sig(1) * sin(theta) + 0.4D1 / 0.3D1 * (0.&
     &2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 &
     &- sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & ** 2 / 0.3D1) * sin(theta)) / 0.4D1 + (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 *&
     & sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) **&
     & 2) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.8D1 * sig(1&
     &) * cos(theta) * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) + 0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * sin(theta) ** &
     &2 - 0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2) * (-sig(2) * sin(the&
     &ta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * sig&
     &(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2) &
     &* (-sig(2) * sin(theta) + sig(3) * cos(theta))) / 0.4D1 + (sig(4) &
     &* cos(theta) + sig(5) * sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1&
     & * (0.4D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig&
     &(1) * sin(theta) ** 2) * (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) - 0.8D1 * sig(1) * cos(theta) * sin(theta) * (sig(4) * cos(the&
     &ta) + sig(5) * sin(theta)) + 0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(&
     &1) * sin(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2) *&
     & (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (-&
     &0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) &
     &* sin(theta) ** 2) * (-sig(4) * sin(theta) + sig(5) * cos(theta)))&
     & / 0.4D1)
      J(4,2) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.4D1 / 0.3D&
     &1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) * sin(theta)) * sig(1) * sin(theta) + 0.4D1 * (cos(theta)&
     & * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - (sig(2) * cos(th&
     &eta) + sig(3) * sin(theta)) * sin(theta)) * sig(1) * cos(theta) - &
     &0.8D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) * cos(theta)) * sig(1) * sin(theta) + 0.4D1 &
     &/ 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) * sin(theta)) * sig(1) * sin(theta)) / 0.4D1 + cos&
     &(theta) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) *&
     &* 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig&
     &(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + si&
     &g(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(&
     &4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 &
     &/ 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.4D1 * &
     &(-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta))) * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(thet&
     &a) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta&
     &) + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1&
     & / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - s&
     &ig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) **&
     & 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) / 0.4D&
     &1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (-0&
     &.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) * sin(theta)) * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * co&
     &s(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * s&
     &in(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2)&
     & * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1)&
     & * sin(theta) + 0.4D1 * (cos(theta) * (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) - (sig(2) * cos(theta) + sig(3) * sin(theta)) * s&
     &in(theta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.4D1 *&
     & (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(t&
     &heta) + sig(5) * cos(theta))) * cos(theta) + 0.8D1 / 0.3D1 * (-0.4&
     &D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(th&
     &eta) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & * cos(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.&
     &8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2&
     &D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0&
     &.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 -&
     & sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sin(theta) - 0.4&
     &D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + si&
     &g(3) * cos(theta)) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig&
     &(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2&
     & / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D&
     &1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * sin(theta)) / &
     &0.4D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1 / rho0 *&
     & (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) &
     &* sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta)&
     & + sig(3) * cos(theta)) * sin(theta)) * (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) + 0.4D1 * (cos(theta) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) - (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) * sin(theta)) * (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.&
     &8D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) * cos(theta)) * (-sig(4) * sin(theta) + sig(5)&
     & * cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) * sin(theta)) * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta))) / 0.4D1)
      J(4,3) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.4D1 / 0.3D&
     &1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta))) * sig(1) * sin(theta) + 0.4D1 * ((-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) * sin(theta) + (sig(2) * cos(th&
     &eta) + sig(3) * sin(theta)) * cos(theta)) * sig(1) * cos(theta) - &
     &0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) * sin(theta)) * sig(1) * sin(theta) + 0.4D1 /&
     & 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta))) * sig(1) * sin(theta)) / 0.4D1 + sin(&
     &theta) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) **&
     & 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 /&
     & 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.4D1 * (&
     &-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) &
     &+ (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta))) * (sig(2) * cos(theta) + sig(3) * sin(&
     &theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta&
     &) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-si&
     &g(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 &
     &/ 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - si&
     &g(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3&
     &) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta&
     &)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** &
     &2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) / 0.4D1&
     & + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (-0.&
     &4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta))) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos&
     &(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(&
     &theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) &
     &* cos(theta) + 0.4D1 * ((-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * co&
     &s(theta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.4D1 * &
     &(-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta))&
     & + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta))) * sin(theta) + 0.8D1 / 0.3D1 * (0.4D1&
     & / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta&
     &)) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) *&
     & sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) + 0.8D&
     &1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1&
     & / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2&
     &D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - s&
     &ig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * cos(theta) - 0.4D1&
     & / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta))) * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1&
     &) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 /&
     & 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 &
     &- (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * cos(theta)) / 0.&
     &4D1 + (sig(4) * cos(theta) + sig(5) * sin(theta)) * mu1 / rho0 * (&
     &-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta))) * (-sig(4) * sin(theta) + sig(&
     &5) * cos(theta)) + 0.4D1 * ((-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* cos(theta)) * (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.8D&
     &1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta) + &
     &sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * &
     &cos(theta)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta))) / 0.4D1)
      J(4,4) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.4D1 / 0.3D&
     &1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos&
     &(theta)) * sin(theta)) * sig(1) * sin(theta) + 0.4D1 * (cos(theta)&
     & * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) * sin(theta)) * sig(1) * cos(theta) - &
     &0.8D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) * cos(theta)) * sig(1) * sin(theta) + 0.4D1 &
     &/ 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(th&
     &eta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) * sin(theta)) * sig(1) * sin(theta)) / 0.4D1 + (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (-0.4D1 / &
     &0.3D1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta&
     &)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) + 0.4D1 * (cos(theta) * (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) - (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(the&
     &ta)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.8D1 / 0.3D1&
     & * (-0.4D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) &
     &* sin(theta) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin&
     &(theta)) * cos(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) * sin(theta)) * (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta))) / 0.4D1 + cos(theta) * mu1 / rho0 * (-0.4D&
     &1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1&
     & / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D&
     &1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig&
     &(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(4) * sin(the&
     &ta) + sig(5) * cos(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * &
     &sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) &
     &* (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.8D1 / 0.3D1 * (0&
     &.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * c&
     &os(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2&
     & / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - s&
     &ig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) **&
     & 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.&
     &3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta))) / 0.4D1 + (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.&
     &3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) + 0&
     &.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(&
     &theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.4D1 / 0&
     &.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.&
     &3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0&
     &.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) *&
     &* 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(the&
     &ta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sin(theta) + 0.4D1 * (c&
     &os(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) * sin(theta)) * (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta&
     &) * sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) +&
     & sig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &))) * cos(theta) + 0.8D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta)) * (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) - 0.8D1 / 0.3D1 * (0.2D1 / 0.3&
     &D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos(theta) &
     &** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / &
     &0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 -&
     & sig(6) ** 2 / 0.3D1) * sin(theta) - 0.4D1 * sin(theta) * sig(6) *&
     &* 2 - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta)) * sin(theta)) * (-sig(4) * sin(theta) &
     &+ sig(5) * cos(theta)) + 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) *&
     &* 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(th&
     &eta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) **&
     & 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0&
     &.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - &
     &(-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * sin(t&
     &heta)) / 0.4D1)
      J(4,5) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.4D1 / 0.3D&
     &1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     & sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta))) * sig(1) * sin(theta) + 0.4D1 * ((-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) * cos(theta)) * sig(1) * cos(theta) - &
     &0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) * sin(theta)) * sig(1) * sin(theta) + 0.4D1 /&
     & 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(th&
     &eta) + sig(5) * cos(theta))) * sig(1) * sin(theta)) / 0.4D1 + (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (-0.4D1 / 0&
     &.3D1 * (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta))) * (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) + 0.4D1 * ((-sig(4) * sin(theta) + sig(5) * cos(theta)) * s&
     &in(theta) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(thet&
     &a)) * (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.8D1 / 0.3D1 &
     &* (0.4D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(&
     &4) * sin(theta) + sig(5) * cos(theta))) * (-sig(2) * sin(theta) + &
     &sig(3) * cos(theta))) / 0.4D1 + sin(theta) * mu1 / rho0 * (-0.4D1 &
     &/ 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 /&
     & 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 &
     &/ 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1&
     &) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(4) * sin(theta&
     &) + sig(5) * cos(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * si&
     &n(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) * &
     &(sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.8D1 / 0.3D1 * (0.2&
     &D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig&
     &(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos&
     &(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 /&
     & 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * c&
     &os(theta)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) &
     &* sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig&
     &(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2&
     & / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D&
     &1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (s&
     &ig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) &
     &* sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta))) / 0.4D1 + (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D&
     &1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) - 0.2&
     &D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(the&
     &ta))) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0.4D1 / 0.3&
     &D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D&
     &1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3&
     &D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** &
     &2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * cos(theta) + 0.4D1 * ((-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) + (sig(4) *&
     & cos(theta) + sig(5) * sin(theta)) * cos(theta)) * (sig(4) * cos(t&
     &heta) + sig(5) * sin(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) &
     &* sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig&
     &(2) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     &) * sin(theta) + 0.8D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * sin(theta)) * (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 &
     &* sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(t&
     &heta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos(theta) ** &
     &2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3&
     &D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - si&
     &g(6) ** 2 / 0.3D1) * cos(theta) + 0.4D1 * cos(theta) * sig(6) ** 2&
     & - 0.4D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4&
     &) * sin(theta) + sig(5) * cos(theta))) * (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2&
     & - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta&
     &) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 &
     &/ 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D&
     &1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * cos(thet&
     &a)) / 0.4D1)
      J(4,6) = -rho0 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) * mu1 /&
     & rho0 * sig(6) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * mu1 / rho0 * sig(6) * (-sig(2) * sin(the&
     &ta) + sig(3) * cos(theta)) + 0.4D1 / 0.3D1 * (sig(4) * cos(theta) &
     &+ sig(5) * sin(theta)) * mu1 / rho0 * (-sig(4) * sin(theta) + sig(&
     &5) * cos(theta)) * sig(6))
      J(5,1) = -rho0 * (-sin(theta) * mu1 / rho0 * (0.4D1 * (sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) * sig(6) * sig(1) * cos(theta) - 0&
     &.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) * sig&
     &(1) * sin(theta)) / 0.4D1 - sig(1) * sin(theta) * mu1 / rho0 * (0.&
     &4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * cos(t&
     &heta) - 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig&
     &(6) * sin(theta)) / 0.4D1 + (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * &
     &cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2) * sig(&
     &6) - 0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * sin(theta) ** 2 - 0&
     &.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2) * sig(6) + 0.8D1 / 0.3D1 &
     &* (-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig&
     &(1) * sin(theta) ** 2) * sig(6)) / 0.4D1)
      J(5,2) = -rho0 * (-sin(theta) * mu1 / rho0 * (0.4D1 * (sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) * sig(6) * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) * sig(6) * (-sig(2) * sin(theta) + sig(3) * cos(theta)))&
     & / 0.4D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 / rh&
     &o0 * (0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6)&
     & * cos(theta) - 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) * sig(6) * sin(theta)) / 0.4D1 + (-sig(4) * sin(theta) + sig(5)&
     & * cos(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (&
     &sig(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) + 0.2D1 / &
     &0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta))&
     & * sig(6) - 0.4D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * &
     &cos(theta) + sig(3) * sin(theta)) * cos(theta)) * sig(6) + 0.8D1 /&
     & 0.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(the&
     &ta)) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) * sin(theta)) * sig(6)) / 0.4D1)
      J(5,3) = -rho0 * (cos(theta) * mu1 / rho0 * (0.4D1 * (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) * sig(6) * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos&
     &(theta)) * sig(6) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) &
     &/ 0.4D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 / rho&
     &0 * (0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) &
     &* sin(theta) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) * sig(6) * cos(theta)) / 0.4D1 + (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (s&
     &ig(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) - 0.2D1 / 0&
     &.3D1 * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) &
     &* sig(6) - 0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) * sin(theta)) * sig(6) + 0.8D1 / 0&
     &.3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(thet&
     &a) + sig(3) * cos(theta))) * sig(6)) / 0.4D1)
      J(5,4) = -rho0 * (-sig(1) * sin(theta) * mu1 / rho0 * (0.4D1 * sig&
     &(1) * cos(theta) ** 2 * sig(6) + 0.4D1 * sig(1) * sin(theta) ** 2 &
     &* sig(6)) / 0.4D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) *&
     & mu1 / rho0 * (0.4D1 * cos(theta) * sig(6) * (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) - 0.4D1 * sin(theta) * sig(6) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta))) / 0.4D1 - sin(theta) * mu1 / rh&
     &o0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) *&
     &* 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta))&
     & ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(the&
     &ta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + &
     &sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) &
     &+ 0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6&
     &) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2&
     & + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) **&
     & 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) &
     &** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta&
     &) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0&
     &.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) &
     &+ 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos&
     &(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) *&
     & sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) +&
     & sig(5) * cos(theta)) ** 2 / 0.3D1) * sig(6)) / 0.4D1 + (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 &
     &* (0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * c&
     &os(theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) * sin(theta)) * sig(6) + 0.8D1 * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) * sig(6) * cos(theta) - 0.4D1 / 0.3D1 * (-0.4D1 &
     &/ 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta&
     &) - 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * &
     &cos(theta)) * sig(6) - 0.8D1 * (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) * sig(6) * sin(theta) + 0.8D1 / 0.3D1 * (-0.2D1 / 0.3D1 &
     &* (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) + 0.2D1&
     & / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(thet&
     &a)) * sig(6)) / 0.4D1)
      J(5,5) = -rho0 * ((-sig(2) * sin(theta) + sig(3) * cos(theta)) * m&
     &u1 / rho0 * (0.4D1 * sin(theta) * (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) * sig(6) + 0.4D1 * cos(theta) * (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta)) * sig(6)) / 0.4D1 + cos(theta) * mu1 / rho0&
     & * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** &
     &2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) **&
     & 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     &* 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + &
     &0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6) &
     &- 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 +&
     & 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2&
     & + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) **&
     & 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4&
     &D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) + &
     &0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(t&
     &heta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(thet&
     &a) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) ** 2 / 0.3D1) * sig(6)) / 0.4D1 + (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * &
     &(0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin&
     &(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig&
     &(5) * cos(theta))) * sig(6) + 0.8D1 * (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * sig(6) * sin(theta) - 0.4D1 / 0.3D1 * (0.4D1 / 0&
     &.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) -&
     & 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin&
     &(theta)) * sig(6) + 0.8D1 * (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) * sig(6) * cos(theta) + 0.8D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) - 0.2D1 / &
     &0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)))&
     & * sig(6)) / 0.4D1)
      J(5,6) = -rho0 * (-sig(1) * sin(theta) * mu1 / rho0 * (0.4D1 * (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) * sig(1) * cos(theta) - 0&
     &.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(1) * sin&
     &(theta)) / 0.4D1 + (-sig(2) * sin(theta) + sig(3) * cos(theta)) * &
     &mu1 / rho0 * (0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &* (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.4D1 * (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) * (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta))) / 0.4D1 + (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) * mu1 / rho0 * (0.8D1 * sig(6) ** 2 - 0.4D1 / 0.3D1 * sig(1)&
     & ** 2 * cos(theta) ** 2 - 0.4D1 / 0.3D1 * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) ** 2 + 0.8D1 / 0.3D1 * (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) ** 2 - 0.4D1 / 0.3D1 * sig(1) ** 2 * sin(thet&
     &a) ** 2 - 0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) ** 2 + 0.8D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2) / 0.4D1)
      J(6,1) = -rho0 * (cos(theta) * mu1 / rho0 * (0.4D1 * (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) * sig(6) * sig(1) * cos(theta) - 0.&
     &4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) * sig(&
     &1) * sin(theta)) / 0.4D1 + sig(1) * cos(theta) * mu1 / rho0 * (0.4&
     &D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * cos(th&
     &eta) - 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(&
     &6) * sin(theta)) / 0.4D1 + (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * co&
     &s(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1) * sin(theta) ** 2) * sig(6)&
     & - 0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * sig(1) * sin(theta) ** 2 - 0.2&
     &D1 / 0.3D1 * sig(1) * cos(theta) ** 2) * sig(6) + 0.8D1 / 0.3D1 * &
     &(-0.2D1 / 0.3D1 * sig(1) * cos(theta) ** 2 - 0.2D1 / 0.3D1 * sig(1&
     &) * sin(theta) ** 2) * sig(6)) / 0.4D1)
      J(6,2) = -rho0 * (cos(theta) * mu1 / rho0 * (0.4D1 * (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) * sig(6) * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos&
     &(theta)) * sig(6) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) &
     &/ 0.4D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0&
     & * (0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) *&
     & cos(theta) - 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & * sig(6) * sin(theta)) / 0.4D1 + (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * cos(theta) + 0.2D1 / 0.3&
     &D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sin(theta)) * &
     &sig(6) - 0.4D1 / 0.3D1 * (-0.4D1 / 0.3D1 * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) * sin(theta) - 0.2D1 / 0.3D1 * (sig(2) * cos&
     &(theta) + sig(3) * sin(theta)) * cos(theta)) * sig(6) + 0.8D1 / 0.&
     &3D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)&
     &) * cos(theta) + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) * sin(theta)) * sig(6)) / 0.4D1)
      J(6,3) = -rho0 * (sin(theta) * mu1 / rho0 * (0.4D1 * (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) * sig(6) * (sig(2) * cos(theta) + s&
     &ig(3) * sin(theta)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos&
     &(theta)) * sig(6) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) &
     &/ 0.4D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0&
     & * (0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) *&
     & sin(theta) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & * sig(6) * cos(theta)) / 0.4D1 + (sig(4) * cos(theta) + sig(5) * &
     &sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * (sig&
     &(2) * cos(theta) + sig(3) * sin(theta)) * sin(theta) - 0.2D1 / 0.3&
     &D1 * cos(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) * &
     &sig(6) - 0.4D1 / 0.3D1 * (0.4D1 / 0.3D1 * cos(theta) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) - 0.2D1 / 0.3D1 * (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) * sin(theta)) * sig(6) + 0.8D1 / 0.3&
     &D1 * (-0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & * sin(theta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(2) * sin(theta)&
     & + sig(3) * cos(theta))) * sig(6)) / 0.4D1)
      J(6,4) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.4D1 * sig(&
     &1) * cos(theta) ** 2 * sig(6) + 0.4D1 * sig(1) * sin(theta) ** 2 *&
     & sig(6)) / 0.4D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * m&
     &u1 / rho0 * (0.4D1 * cos(theta) * sig(6) * (sig(2) * cos(theta) + &
     &sig(3) * sin(theta)) - 0.4D1 * sin(theta) * sig(6) * (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta))) / 0.4D1 + cos(theta) * mu1 / rho0&
     & * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** &
     &2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) **&
     & 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) *&
     &* 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + &
     &0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6) &
     &- 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 +&
     & 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2&
     & + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) **&
     & 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4&
     &D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) + &
     &0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(t&
     &heta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * s&
     &in(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(thet&
     &a) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) ** 2 / 0.3D1) * sig(6)) / 0.4D1 + (sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (&
     &0.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(&
     &theta) + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) * sin(theta)) * sig(6) + 0.8D1 * (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * sig(6) * cos(theta) - 0.4D1 / 0.3D1 * (-0.4D1 / 0&
     &.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta) -&
     & 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * cos&
     &(theta)) * sig(6) - 0.8D1 * (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) * sig(6) * sin(theta) + 0.8D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (&
     &sig(4) * cos(theta) + sig(5) * sin(theta)) * cos(theta) + 0.2D1 / &
     &0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sin(theta))&
     & * sig(6)) / 0.4D1)
      J(6,5) = -rho0 * ((sig(2) * cos(theta) + sig(3) * sin(theta)) * mu&
     &1 / rho0 * (0.4D1 * sin(theta) * (sig(2) * cos(theta) + sig(3) * s&
     &in(theta)) * sig(6) + 0.4D1 * cos(theta) * (-sig(2) * sin(theta) +&
     & sig(3) * cos(theta)) * sig(6)) / 0.4D1 + sin(theta) * mu1 / rho0 &
     &* (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2&
     & + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** &
     &2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) **&
     & 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta)&
     & + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0&
     &.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6) -&
     & 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + &
     &0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 &
     &+ 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** &
     &2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) +&
     & sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D&
     &1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) + 0&
     &.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(th&
     &eta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta&
     &) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) ** 2 / 0.3D1) * sig(6)) / 0.4D1 + (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0&
     &.4D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(t&
     &heta) - 0.2D1 / 0.3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta))) * sig(6) + 0.8D1 * (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * sig(6) * sin(theta) - 0.4D1 / 0.3D1 * (0.4D1 / 0.3&
     &D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) - 0&
     &.2D1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * sin(t&
     &heta)) * sig(6) + 0.8D1 * (-sig(4) * sin(theta) + sig(5) * cos(the&
     &ta)) * sig(6) * cos(theta) + 0.8D1 / 0.3D1 * (-0.2D1 / 0.3D1 * (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) * sin(theta) - 0.2D1 / 0.&
     &3D1 * cos(theta) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) *&
     & sig(6)) / 0.4D1)
      J(6,6) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.4D1 * (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * sig(1) * cos(theta) - 0.&
     &4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(1) * sin(&
     &theta)) / 0.4D1 + (sig(2) * cos(theta) + sig(3) * sin(theta)) * mu&
     &1 / rho0 * (0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) * &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.4D1 * (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta))) / 0.4D1 + (sig(4) * cos(theta) + sig(5) * sin(theta&
     &)) * mu1 / rho0 * (0.8D1 * sig(6) ** 2 - 0.4D1 / 0.3D1 * sig(1) **&
     & 2 * cos(theta) ** 2 - 0.4D1 / 0.3D1 * (sig(2) * cos(theta) + sig(&
     &3) * sin(theta)) ** 2 + 0.8D1 / 0.3D1 * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) ** 2 - 0.4D1 / 0.3D1 * sig(1) ** 2 * sin(theta) &
     &** 2 - 0.4D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)&
     &) ** 2 + 0.8D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) ** 2) / 0.4D1)   
        case(2)
 J(1,1) = -rho0 * (mu1 / rho0 * (0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * s&
     &ig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) &
     &** 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0&
     &.3D1) * sig(1) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3) ** 2 + 0.&
     &2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.&
     &3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) - 0.4D1 &
     &/ 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig&
     &(2) ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig&
     &(5) ** 2 / 0.3D1) * sig(1)) / 0.4D1 + sig(1) * mu1 / rho0 * (0.8D1&
     & * sig(1) ** 2 + 0.8D1 / 0.3D1 * sig(2) ** 2 + 0.8D1 / 0.3D1 * sig&
     &(4) ** 2 - 0.4D1 / 0.3D1 * sig(3) ** 2 - 0.4D1 / 0.3D1 * sig(5) **&
     & 2 - 0.4D1 / 0.3D1 * sig(6) ** 2) / 0.4D1 + 0.4D1 / 0.3D1 * sig(2)&
     & ** 2 * mu1 / rho0 * sig(1) + 0.4D1 / 0.3D1 * sig(4) ** 2 * mu1 / &
     &rho0 * sig(1)) + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) * sig&
     &(3) ** 2 * sig(6) ** 2 * (-sig(1) * sig(3) * sig(6) + 0.1D1) - (lambda1&
     & + 0.2D1 / 0.3D1 * mu1) * sig(1) ** 2 * sig(3) ** 3 * sig(6) ** 3 &
     &- rho0 * T0 * (sig(1) * sig(3) * sig(6)) ** g * g / sig(1) * (exp(&
     &s / cv) - 0.1D1)
      J(1,2) = -rho0 * (0.4D1 / 0.3D1 * sig(1) ** 2 * mu1 / rho0 * sig(2&
     &) + mu1 / rho0 * (0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0&
     &.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) &
     &** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(2)&
     & + 0.4D1 * (sig(3) * sig(2) + sig(4) * sig(5)) * sig(3) - 0.4D1 / &
     &0.3D1 * (0.2D1 / 0.3D1 * sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2&
     & - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1&
     & - sig(6) ** 2 / 0.3D1) * sig(2) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 &
     &* sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(3)&
     & ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) * sig(2&
     &)) / 0.4D1 + sig(2) * mu1 / rho0 * (0.8D1 * sig(2) ** 2 + 0.8D1 / &
     &0.3D1 * sig(1) ** 2 + 0.8D1 / 0.3D1 * sig(4) ** 2 + 0.8D1 / 0.3D1 &
     &* sig(3) ** 2 - 0.4D1 / 0.3D1 * sig(5) ** 2 - 0.4D1 / 0.3D1 * sig(&
     &6) ** 2) / 0.4D1 + sig(4) * mu1 / rho0 * (0.16D2 / 0.3D1 * sig(2) &
     &* sig(4) + 0.4D1 * sig(3) * sig(5)) / 0.4D1)
      J(1,3) = -rho0 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * mu1 / rho0 * sig(&
     &3) + sig(2) * mu1 / rho0 * (0.16D2 / 0.3D1 * sig(3) * sig(2) + 0.4&
     &D1 * sig(4) * sig(5)) / 0.4D1 + sig(4) * mu1 / rho0 * (-0.8D1 / 0.&
     &3D1 * sig(3) * sig(4) + 0.4D1 * sig(2) * sig(5)) / 0.4D1) + 0.2D1 &
     &* (lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) ** 2 * sig(3) * sig(6) ** 2 &
     &* (-sig(1) * sig(3) * sig(6) + 0.1D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1)&
     & * sig(1) ** 3 * sig(3) ** 2 * sig(6) ** 3 - rho0 * T0 * (sig(1) *&
     & sig(3) * sig(6)) ** g * g / sig(3) * (exp(s / cv) - 0.1D1)
      J(1,4) = -rho0 * (0.4D1 / 0.3D1 * sig(1) ** 2 * mu1 / rho0 * sig(4&
     &) + sig(2) * mu1 / rho0 * (0.16D2 / 0.3D1 * sig(2) * sig(4) + 0.4D&
     &1 * sig(3) * sig(5)) / 0.4D1 + mu1 / rho0 * (0.8D1 / 0.3D1 * (0.2D&
     &1 / 0.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.&
     &3D1 * sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - si&
     &g(6) ** 2 / 0.3D1) * sig(4) + 0.4D1 * (sig(3) * sig(2) + sig(4) * &
     &sig(5)) * sig(5) + 0.4D1 * sig(4) * sig(6) ** 2 - 0.4D1 / 0.3D1 * &
     &(0.2D1 / 0.3D1 * sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1&
     &) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6&
     &) ** 2 / 0.3D1) * sig(4) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6)&
     & ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(3) ** 2 / &
     &0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) * sig(4)) / 0.4&
     &D1 + sig(4) * mu1 / rho0 * (0.8D1 * sig(4) ** 2 + 0.8D1 / 0.3D1 * &
     &sig(1) ** 2 + 0.8D1 / 0.3D1 * sig(2) ** 2 - 0.4D1 / 0.3D1 * sig(3)&
     & ** 2 + 0.8D1 / 0.3D1 * sig(5) ** 2 + 0.8D1 / 0.3D1 * sig(6) ** 2)&
     & / 0.4D1)
      J(1,5) = -rho0 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * mu1 / rho0 * sig(&
     &5) + sig(2) * mu1 / rho0 * (-0.8D1 / 0.3D1 * sig(2) * sig(5) + 0.4&
     &D1 * sig(3) * sig(4)) / 0.4D1 + sig(4) * mu1 / rho0 * (0.16D2 / 0.&
     &3D1 * sig(4) * sig(5) + 0.4D1 * sig(3) * sig(2)) / 0.4D1)
      J(1,6) = -rho0 * (-0.2D1 / 0.3D1 * sig(1) ** 2 * mu1 / rho0 * sig(&
     &6) - 0.2D1 / 0.3D1 * sig(2) ** 2 * mu1 / rho0 * sig(6) + 0.4D1 / 0&
     &.3D1 * sig(4) ** 2 * mu1 / rho0 * sig(6)) + 0.2D1 * (lambda1 + 0.2D1 / &
     &0.3D1 * mu1) * sig(1) ** 2 * sig(3) ** 2 * sig(6) * (-sig(1) * sig&
     &(3) * sig(6) + 0.1D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) ** 3 *&
     & sig(3) ** 3 * sig(6) ** 2 - rho0 * T0 * (sig(1) * sig(3) * sig(6)&
     &) ** g * g / sig(6) * (exp(s / cv) - 0.1D1)
      J(2,1) = -rho0 * (-0.2D1 / 0.3D1 * sig(3) ** 2 * mu1 / rho0 * sig(&
     &1) - 0.2D1 / 0.3D1 * sig(5) ** 2 * mu1 / rho0 * sig(1)) + 0.2D1 * &
     &(lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) * sig(3) ** 2 * sig(6) ** 2 * &
     &(-sig(1) * sig(3) * sig(6) + 0.1D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) *&
     & sig(1) ** 2 * sig(3) ** 3 * sig(6) ** 3 - rho0 * T0 * (sig(1) * s&
     &ig(3) * sig(6)) ** g * g / sig(1) * (exp(s / cv) - 0.1D1)
      J(2,2) = -rho0 * (sig(3) * mu1 / rho0 * (0.16D2 / 0.3D1 * sig(3) *&
     & sig(2) + 0.4D1 * sig(4) * sig(5)) / 0.4D1 + sig(5) * mu1 / rho0 *&
     & (-0.8D1 / 0.3D1 * sig(2) * sig(5) + 0.4D1 * sig(3) * sig(4)) / 0.&
     &4D1)
      J(2,3) = -rho0 * (mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * &
     &sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4)&
     & ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / &
     &0.3D1) * sig(3) + 0.4D1 * (sig(3) * sig(2) + sig(4) * sig(5)) * si&
     &g(2) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3) ** 2 + 0.2D1 / 0.3D&
     &1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(&
     &4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(3) - 0.4D1 / 0.3D1 * &
     &(0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 /&
     & 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 /&
     & 0.3D1) * sig(3)) / 0.4D1 + sig(3) * mu1 / rho0 * (0.8D1 * sig(3) &
     &** 2 - 0.4D1 / 0.3D1 * sig(1) ** 2 + 0.8D1 / 0.3D1 * sig(2) ** 2 -&
     & 0.4D1 / 0.3D1 * sig(4) ** 2 + 0.8D1 / 0.3D1 * sig(5) ** 2 - 0.4D1&
     & / 0.3D1 * sig(6) ** 2) / 0.4D1 + sig(5) * mu1 / rho0 * (0.16D2 / &
     &0.3D1 * sig(3) * sig(5) + 0.4D1 * sig(2) * sig(4)) / 0.4D1) + 0.2D&
     &1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) ** 2 * sig(3) * sig(6) ** &
     &2 * (-sig(1) * sig(3) * sig(6) + 0.1D1) - (lambda1 + 0.2D1 / 0.3D1 * mu&
     &1) * sig(1) ** 3 * sig(3) ** 2 * sig(6) ** 3 - rho0 * T0 * (sig(1)&
     & * sig(3) * sig(6)) ** g * g / sig(3) * (exp(s / cv) - 0.1D1)
      J(2,4) = -rho0 * (sig(3) * mu1 / rho0 * (-0.8D1 / 0.3D1 * sig(3) *&
     & sig(4) + 0.4D1 * sig(2) * sig(5)) / 0.4D1 + sig(5) * mu1 / rho0 *&
     & (0.16D2 / 0.3D1 * sig(4) * sig(5) + 0.4D1 * sig(3) * sig(2)) / 0.&
     &4D1)
      J(2,5) = -rho0 * (sig(3) * mu1 / rho0 * (0.16D2 / 0.3D1 * sig(3) *&
     & sig(5) + 0.4D1 * sig(2) * sig(4)) / 0.4D1 + mu1 / rho0 * (-0.4D1 &
     &/ 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) **&
     & 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) **&
     & 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(5) + 0.4D1 * (sig(3) * sig&
     &(2) + sig(4) * sig(5)) * sig(4) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 *&
     & sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 -&
     & sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) &
     &* sig(5) + 0.4D1 * sig(5) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 /&
     & 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 -&
     & sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) &
     &* sig(5)) / 0.4D1 + sig(5) * mu1 / rho0 * (0.8D1 * sig(5) ** 2 - 0&
     &.4D1 / 0.3D1 * sig(1) ** 2 - 0.4D1 / 0.3D1 * sig(2) ** 2 + 0.8D1 /&
     & 0.3D1 * sig(4) ** 2 + 0.8D1 / 0.3D1 * sig(3) ** 2 + 0.8D1 / 0.3D1&
     & * sig(6) ** 2) / 0.4D1)
      J(2,6) = -rho0 * (-0.2D1 / 0.3D1 * sig(3) ** 2 * mu1 / rho0 * sig(&
     &6) + 0.4D1 / 0.3D1 * sig(5) ** 2 * mu1 / rho0 * sig(6)) + 0.2D1 * &
     &(lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) ** 2 * sig(3) ** 2 * sig(6) * &
     &(-sig(1) * sig(3) * sig(6) + 0.1D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) *&
     & sig(1) ** 3 * sig(3) ** 3 * sig(6) ** 2 - rho0 * T0 * (sig(1) * s&
     &ig(3) * sig(6)) ** g * g / sig(6) * (exp(s / cv) - 0.1D1)
      J(3,1) = 0.2D1 / 0.3D1 * sig(6) ** 2 * mu1 * sig(1) + 0.2D1 * (lambda1 &
     &+ 0.2D1 / 0.3D1 * mu1) * sig(1) * sig(3) ** 2 * sig(6) ** 2 * (-si&
     &g(1) * sig(3) * sig(6) + 0.1D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) * sig&
     &(1) ** 2 * sig(3) ** 3 * sig(6) ** 3 - rho0 * T0 * (sig(1) * sig(3&
     &) * sig(6)) ** g * g / sig(1) * (exp(s / cv) - 0.1D1)
      J(3,2) = 0.2D1 / 0.3D1 * sig(6) ** 2 * mu1 * sig(2)
      J(3,3) = 0.2D1 / 0.3D1 * sig(6) ** 2 * mu1 * sig(3) + 0.2D1 * (lambda1 &
     &+ 0.2D1 / 0.3D1 * mu1) * sig(1) ** 2 * sig(3) * sig(6) ** 2 * (-si&
     &g(1) * sig(3) * sig(6) + 0.1D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) * sig&
     &(1) ** 3 * sig(3) ** 2 * sig(6) ** 3 - rho0 * T0 * (sig(1) * sig(3&
     &) * sig(6)) ** g * g / sig(3) * (exp(s / cv) - 0.1D1)
      J(3,4) = -0.4D1 / 0.3D1 * sig(6) ** 2 * mu1 * sig(4)
      J(3,5) = -0.4D1 / 0.3D1 * sig(6) ** 2 * mu1 * sig(5)
      J(3,6) = -mu1 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0&
     &.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) &
     &** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6)&
     & + 0.4D1 * sig(4) ** 2 * sig(6) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 *&
     & sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 -&
     & sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) &
     &* sig(6) + 0.4D1 * sig(5) ** 2 * sig(6) + 0.8D1 / 0.3D1 * (0.2D1 /&
     & 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 -&
     & sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) &
     &* sig(6)) / 0.4D1 - sig(6) * mu1 * (0.8D1 * sig(6) ** 2 - 0.4D1 / &
     &0.3D1 * sig(1) ** 2 - 0.4D1 / 0.3D1 * sig(2) ** 2 + 0.8D1 / 0.3D1 &
     &* sig(4) ** 2 - 0.4D1 / 0.3D1 * sig(3) ** 2 + 0.8D1 / 0.3D1 * sig(&
     &5) ** 2) / 0.4D1 + 0.2D1 * (lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) ** &
     &2 * sig(3) ** 2 * sig(6) * (-sig(1) * sig(3) * sig(6) + 0.1D1) - (&
     &lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) ** 3 * sig(3) ** 3 * sig(6) ** &
     &2 - rho0 * T0 * (sig(1) * sig(3) * sig(6)) ** g * g / sig(6) * (ex&
     &p(s / cv) - 0.1D1)
      J(4,1) = -rho0 * (0.2D1 * sig(1) * mu1 / rho0 * (sig(3) * sig(2) +&
     & sig(4) * sig(5)) - 0.2D1 / 0.3D1 * sig(2) * mu1 / rho0 * sig(3) *&
     & sig(1) - 0.2D1 / 0.3D1 * sig(4) * mu1 / rho0 * sig(5) * sig(1))
      J(4,2) = -rho0 * (sig(1) ** 2 * mu1 / rho0 * sig(3) + mu1 / rho0 *&
     & (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * &
     &sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - &
     &sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(3) + 0.4D1 * (sig&
     &(3) * sig(2) + sig(4) * sig(5)) * sig(2) + 0.8D1 / 0.3D1 * (0.2D1 &
     &/ 0.3D1 * sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 &
     &/ 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 &
     &/ 0.3D1) * sig(3) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 -&
     & sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 -&
     & sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) * sig(3)) / 0.4D1 + si&
     &g(2) * mu1 / rho0 * (0.16D2 / 0.3D1 * sig(3) * sig(2) + 0.4D1 * si&
     &g(4) * sig(5)) / 0.4D1 + sig(4) * mu1 / rho0 * (-0.8D1 / 0.3D1 * s&
     &ig(2) * sig(5) + 0.4D1 * sig(3) * sig(4)) / 0.4D1)
      J(4,3) = -rho0 * (sig(1) ** 2 * mu1 / rho0 * sig(2) + sig(2) * mu1&
     & / rho0 * (0.8D1 * sig(3) ** 2 - 0.4D1 / 0.3D1 * sig(1) ** 2 + 0.8&
     &D1 / 0.3D1 * sig(2) ** 2 - 0.4D1 / 0.3D1 * sig(4) ** 2 + 0.8D1 / 0&
     &.3D1 * sig(5) ** 2 - 0.4D1 / 0.3D1 * sig(6) ** 2) / 0.4D1 + sig(4)&
     & * mu1 / rho0 * (0.16D2 / 0.3D1 * sig(3) * sig(5) + 0.4D1 * sig(2)&
     & * sig(4)) / 0.4D1)
      J(4,4) = -rho0 * (sig(1) ** 2 * mu1 / rho0 * sig(5) + sig(2) * mu1&
     & / rho0 * (-0.8D1 / 0.3D1 * sig(3) * sig(4) + 0.4D1 * sig(2) * sig&
     &(5)) / 0.4D1 + mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig&
     &(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) **&
     & 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3&
     &D1) * sig(5) + 0.4D1 * (sig(3) * sig(2) + sig(4) * sig(5)) * sig(4&
     &) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3) ** 2 + 0.2D1 / 0.3D1 *&
     & sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(4) &
     &** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(5) + 0.4D1 * sig(5) * si&
     &g(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) &
     &** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4) &
     &** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) * sig(5)) / 0.4D1 + sig(4) * m&
     &u1 / rho0 * (0.16D2 / 0.3D1 * sig(4) * sig(5) + 0.4D1 * sig(3) * s&
     &ig(2)) / 0.4D1)
      J(4,5) = -rho0 * (sig(1) ** 2 * mu1 / rho0 * sig(4) + sig(2) * mu1&
     & / rho0 * (0.16D2 / 0.3D1 * sig(3) * sig(5) + 0.4D1 * sig(2) * sig&
     &(4)) / 0.4D1 + sig(4) * mu1 / rho0 * (0.8D1 * sig(5) ** 2 - 0.4D1 &
     &/ 0.3D1 * sig(1) ** 2 - 0.4D1 / 0.3D1 * sig(2) ** 2 + 0.8D1 / 0.3D&
     &1 * sig(4) ** 2 + 0.8D1 / 0.3D1 * sig(3) ** 2 + 0.8D1 / 0.3D1 * si&
     &g(6) ** 2) / 0.4D1)
      J(4,6) = -rho0 * (-0.2D1 / 0.3D1 * sig(2) * mu1 / rho0 * sig(3) * &
     &sig(6) + 0.4D1 / 0.3D1 * sig(4) * mu1 / rho0 * sig(5) * sig(6))
      J(5,1) = 0.2D1 / 0.3D1 * sig(5) * mu1 * sig(1) * sig(6)
      J(5,2) = -rho0 * (sig(3) * mu1 / rho0 * sig(4) * sig(6) - 0.2D1 / &
     &0.3D1 * sig(5) * mu1 / rho0 * sig(6) * sig(2))
      J(5,3) = -rho0 * (mu1 / rho0 * (0.4D1 * sig(4) * sig(6) * sig(2) +&
     & 0.4D1 * sig(5) * sig(6) * sig(3)) / 0.4D1 + sig(3) * mu1 / rho0 *&
     & sig(5) * sig(6) / 0.3D1)
      J(5,4) = -rho0 * (sig(2) * mu1 / rho0 * sig(3) * sig(6) + 0.4D1 / &
     &0.3D1 * sig(4) * mu1 / rho0 * sig(5) * sig(6))
      J(5,5) = -rho0 * (sig(3) ** 2 * mu1 / rho0 * sig(6) + mu1 / rho0 *&
     & (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * &
     &sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - &
     &sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D1 * sig(&
     &4) ** 2 * sig(6) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3) ** 2 + &
     &0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / &
     &0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D&
     &1 * sig(5) ** 2 * sig(6) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6)&
     & ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(3) ** 2 / &
     &0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) * sig(6)) / 0.4&
     &D1 + 0.4D1 / 0.3D1 * sig(5) ** 2 * mu1 / rho0 * sig(6))
      J(5,6) = -rho0 * (sig(3) * mu1 / rho0 * (0.4D1 * sig(2) * sig(4) +&
     & 0.4D1 * sig(3) * sig(5)) / 0.4D1 + sig(5) * mu1 / rho0 * (0.8D1 *&
     & sig(6) ** 2 - 0.4D1 / 0.3D1 * sig(1) ** 2 - 0.4D1 / 0.3D1 * sig(2&
     &) ** 2 + 0.8D1 / 0.3D1 * sig(4) ** 2 - 0.4D1 / 0.3D1 * sig(3) ** 2&
     & + 0.8D1 / 0.3D1 * sig(5) ** 2) / 0.4D1)
      J(6,1) = -0.4D1 / 0.3D1 * sig(4) * mu1 * sig(1) * sig(6)
      J(6,2) = -rho0 * (mu1 / rho0 * (0.4D1 * sig(4) * sig(6) * sig(2) +&
     & 0.4D1 * sig(5) * sig(6) * sig(3)) / 0.4D1 + sig(2) * mu1 / rho0 *&
     & sig(4) * sig(6) / 0.3D1)
      J(6,3) = -rho0 * (sig(5) * mu1 / rho0 * sig(6) * sig(2) - 0.2D1 / &
     &0.3D1 * sig(3) * mu1 / rho0 * sig(4) * sig(6))
      J(6,4) = -rho0 * (sig(1) ** 2 * mu1 / rho0 * sig(6) + sig(2) ** 2 &
     &* mu1 / rho0 * sig(6) + mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.&
     &3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 * &
     &sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) *&
     &* 2 / 0.3D1) * sig(6) + 0.4D1 * sig(4) ** 2 * sig(6) - 0.4D1 / 0.3&
     &D1 * (0.2D1 / 0.3D1 * sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - &
     &sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - &
     &sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D1 * sig(5) ** 2 * sig(6) + 0.8&
     &D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - &
     &sig(2) ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - &
     &sig(5) ** 2 / 0.3D1) * sig(6)) / 0.4D1 + 0.4D1 / 0.3D1 * sig(4) **&
     & 2 * mu1 / rho0 * sig(6))
      J(6,5) = -rho0 * (sig(2) * mu1 / rho0 * sig(3) * sig(6) + 0.4D1 / &
     &0.3D1 * sig(4) * mu1 / rho0 * sig(5) * sig(6))
      J(6,6) = -rho0 * (sig(1) ** 2 * mu1 / rho0 * sig(4) + sig(2) * mu1&
     & / rho0 * (0.4D1 * sig(2) * sig(4) + 0.4D1 * sig(3) * sig(5)) / 0.&
     &4D1 + sig(4) * mu1 / rho0 * (0.8D1 * sig(6) ** 2 - 0.4D1 / 0.3D1 *&
     & sig(1) ** 2 - 0.4D1 / 0.3D1 * sig(2) ** 2 + 0.8D1 / 0.3D1 * sig(4&
     &) ** 2 - 0.4D1 / 0.3D1 * sig(3) ** 2 + 0.8D1 / 0.3D1 * sig(5) ** 2&
     &) / 0.4D1)           
        case default
            print *, 'EOS_mode in A2sigmaJacobianSGEOS', EOS_mode, 'not implemented!'
        end select
    end subroutine A2sigmaJacobianSGEOS
    ! *******************************************************************************************************
    subroutine A2sigmaComputeFSGEOS(fp,sig,theta,lambda1,mu1,rho0,g,s,cv,p0,EOS_mode, T0)
        implicit none
        real,intent(out) :: fp(6)
        real, intent(in) :: sig(6),theta,lambda1,mu1,rho0,g,s,cv,p0, T0
        integer, intent(in) :: EOS_mode
        ! Computed from MAPLE
        select case(EOS_mode)
            case(1)
     fp(1) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.8D1 / 0.3D1&
     & * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 &
     &* (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1&
     & * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 &
     &* sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * cos(theta) - 0.4D1&
     & * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta))) * sig(1) * sin(theta) - 0.4D1 / 0.&
     &3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3&
     &D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0&
     &.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) &
     &** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * s&
     &in(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * cos(theta) - 0.&
     &4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(the&
     &ta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2)&
     & * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) ** 2 / 0.3D1) * sig(1) * cos(theta)) / 0.4D1 + (&
     &sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (0.8D1 /&
     & 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / &
     &0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 /&
     & 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1)&
     & ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(&
     &theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) * &
     &sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) * (-&
     &sig(2) * sin(theta) + sig(3) * cos(theta)) - 0.4D1 / 0.3D1 * (0.2D&
     &1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(&
     &2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos(&
     &theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / &
     &0.3D1 - sig(6) ** 2 / 0.3D1) * (sig(2) * cos(theta) + sig(3) * sin&
     &(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) *&
     &* 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.&
     &3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * (sig(2) * cos(theta&
     &) + sig(3) * sin(theta))) / 0.4D1 + (sig(4) * cos(theta) + sig(5) &
     &* sin(theta)) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig&
     &(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3&
     &D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) &
     &** 2 / 0.3D1) * (sig(4) * cos(theta) + sig(5) * sin(theta)) + 0.4D&
     &1 * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta))) * (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) + 0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)&
     &) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * s&
     &in(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5) &
     &* cos(theta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(&
     &2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * co&
     &s(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1&
     &) * (sig(4) * cos(theta) + sig(5) * sin(theta)) - 0.4D1 / 0.3D1 * &
     &(0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3&
     &D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) ** 2 / 0.3D1) * (sig(4) * cos(theta) + sig(5) * sin(theta))) /&
     & 0.4D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) / g * (sig(1) * cos(theta) * &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * s&
     &in(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) &
     &** g * exp(s / cv) + (lambda1 + 0.2D1 / 0.3D1 * mu1) / g - p0
      fp(2) = -rho0 * (-sig(1) * sin(theta) * mu1 / rho0 * (0.4D1 / 0.3D&
     &1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1&
     & * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D&
     &1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2&
     & * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(&
     &theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * sin(theta) + 0.4D&
     &1 * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta)&
     & + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(the&
     &ta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * si&
     &n(theta) + sig(5) * cos(theta))) * sig(1) * cos(theta) - 0.8D1 / 0&
     &.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.&
     &3D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / &
     &0.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1)&
     & ** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * &
     &sin(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(the&
     &ta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * sin(theta) + 0&
     &.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(th&
     &eta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2&
     &) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * si&
     &n(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta&
     &) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) ** 2 / 0.3D1) * sig(1) * sin(theta)) / 0.4D1 + &
     &(-sig(2) * sin(theta) + sig(3) * cos(theta)) * mu1 / rho0 * (-0.4D&
     &1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1&
     & / 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D&
     &1 / 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig&
     &(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3&
     &) * cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(2) * sin(the&
     &ta) + sig(3) * cos(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * &
     &sin(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2&
     &) * sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) &
     &* (sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.8D1 / 0.3D1 * (0&
     &.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-s&
     &ig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (&
     &-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * c&
     &os(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta&
     &)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2&
     & / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(&
     &1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 &
     &/ 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1&
     & - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (si&
     &g(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) *&
     & sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * (-sig(2) * sin(&
     &theta) + sig(3) * cos(theta))) / 0.4D1 + (-sig(4) * sin(theta) + s&
     &ig(5) * cos(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D&
     &1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(&
     &theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos&
     &(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** &
     &2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.&
     &3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - &
     &sig(6) ** 2 / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)&
     &) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) *&
     & cos(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta))) * (sig(4) * cos(theta) + &
     &sig(5) * sin(theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** &
     &2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / &
     &0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.4D1 * (-&
     &sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) ** 2 - 0.4D1 /&
     & 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) *&
     &* 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * co&
     &s(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(thet&
     &a) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + si&
     &g(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) *&
     & cos(theta)) ** 2 / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta))) / 0.4D1) - (lambda1 + 0.2D1 / 0.3D1 * mu1) / g * (sig(1) * cos&
     &(theta) * (-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + &
     &sig(1) * sin(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) &
     &* sig(6)) ** g * exp(s / cv) + (lambda1 + 0.2D1 / 0.3D1 * mu1) / g - p0
      fp(3) = -sig(6) * mu1 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) &
     &** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(theta) + si&
     &g(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(theta) + s&
     &ig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 -&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig(6) ** 2&
     & / 0.3D1) * sig(6) + 0.4D1 * (sig(4) * cos(theta) + sig(5) * sin(t&
     &heta)) ** 2 * sig(6) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** &
     &2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + si&
     &g(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4)&
     & * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / &
     &0.3D1) * sig(6) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(the&
     &ta)) ** 2 * sig(6) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 &
     &- sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta)&
     & ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 /&
     & 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1&
     & - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * sig(6)) /&
     & 0.4D1 - (lambda1 + 0.2D1 / 0.3D1 * mu1) / g * (sig(1) * cos(theta) * (&
     &-sig(2) * sin(theta) + sig(3) * cos(theta)) * sig(6) + sig(1) * si&
     &n(theta) * (sig(2) * cos(theta) + sig(3) * sin(theta)) * sig(6)) *&
     &* g * exp(s / cv) + (lambda1 + 0.2D1 / 0.3D1 * mu1) / g - p0
      fp(4) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.4D1 / 0.3D1&
     & * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 &
     &* (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1&
     & * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 &
     &* sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(t&
     &heta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta))&
     & ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * sin(theta) + 0.4D1&
     & * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * cos(thet&
     &a)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4) * sin&
     &(theta) + sig(5) * cos(theta))) * sig(1) * cos(theta) - 0.8D1 / 0.&
     &3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3&
     &D1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0&
     &.3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) &
     &** 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * s&
     &in(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(thet&
     &a)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) * sin(theta) + 0.&
     &4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(the&
     &ta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2)&
     & * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin&
     &(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta)&
     & + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig&
     &(5) * cos(theta)) ** 2 / 0.3D1) * sig(1) * sin(theta)) / 0.4D1 + (&
     &sig(2) * cos(theta) + sig(3) * sin(theta)) * mu1 / rho0 * (-0.4D1 &
     &/ 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 /&
     & 0.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 &
     &/ 0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1&
     &) ** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(&
     &theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(2) * sin(theta&
     &) + sig(3) * cos(theta)) + 0.4D1 * (-sig(1) ** 2 * cos(theta) * si&
     &n(theta) + (sig(2) * cos(theta) + sig(3) * sin(theta)) * (-sig(2) &
     &* sin(theta) + sig(3) * cos(theta)) + (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) * (-sig(4) * sin(theta) + sig(5) * cos(theta))) * &
     &(sig(2) * cos(theta) + sig(3) * sin(theta)) + 0.8D1 / 0.3D1 * (0.2&
     &D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig&
     &(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-s&
     &ig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2 * cos&
     &(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta))&
     & ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 /&
     & 0.3D1 - sig(6) ** 2 / 0.3D1) * (-sig(2) * sin(theta) + sig(3) * c&
     &os(theta)) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1)&
     & ** 2 * cos(theta) ** 2 / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / &
     &0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 -&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * s&
     &in(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1) * (-sig(2) * sin(th&
     &eta) + sig(3) * cos(theta))) / 0.4D1 + (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 *&
     & sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D1 * (sig(2) * cos(the&
     &ta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3D1 * (sig(4) * cos(th&
     &eta) + sig(5) * sin(theta)) ** 2 - sig(1) ** 2 * sin(theta) ** 2 /&
     & 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 / 0.3D1&
     & - (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 / 0.3D1 - sig&
     &(6) ** 2 / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) +&
     & 0.4D1 * (-sig(1) ** 2 * cos(theta) * sin(theta) + (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) * (-sig(2) * sin(theta) + sig(3) * co&
     &s(theta)) + (sig(4) * cos(theta) + sig(5) * sin(theta)) * (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta))) * (sig(4) * cos(theta) + sig&
     &(5) * sin(theta)) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 *&
     & sin(theta) ** 2 + 0.2D1 / 0.3D1 * (-sig(2) * sin(theta) + sig(3) &
     &* cos(theta)) ** 2 + 0.2D1 / 0.3D1 * (-sig(4) * sin(theta) + sig(5&
     &) * cos(theta)) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3D1 - (si&
     &g(2) * cos(theta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (sig(4) * &
     &cos(theta) + sig(5) * sin(theta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3&
     &D1) * (-sig(4) * sin(theta) + sig(5) * cos(theta)) + 0.4D1 * (-sig&
     &(4) * sin(theta) + sig(5) * cos(theta)) * sig(6) ** 2 - 0.4D1 / 0.&
     &3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2&
     & / 0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(t&
     &heta) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) &
     &+ sig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5&
     &) * sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * co&
     &s(theta)) ** 2 / 0.3D1) * (-sig(4) * sin(theta) + sig(5) * cos(the&
     &ta))) / 0.4D1)
      fp(5) = -rho0 * (-sig(1) * sin(theta) * mu1 / rho0 * (0.4D1 * (sig&
     &(4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * sig(1) * cos(th&
     &eta) - 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(&
     &6) * sig(1) * sin(theta)) / 0.4D1 + (-sig(2) * sin(theta) + sig(3)&
     & * cos(theta)) * mu1 / rho0 * (0.4D1 * (sig(4) * cos(theta) + sig(&
     &5) * sin(theta)) * sig(6) * (sig(2) * cos(theta) + sig(3) * sin(th&
     &eta)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig&
     &(6) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) / 0.4D1 + (-si&
     &g(4) * sin(theta) + sig(5) * cos(theta)) * mu1 / rho0 * (-0.4D1 / &
     &0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0&
     &.3D1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / &
     &0.3D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) &
     &** 2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * &
     &cos(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(th&
     &eta)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D1 * (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6) - 0.4D1 / 0.3&
     &D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D&
     &1 * (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.&
     &3D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) *&
     &* 2 * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * si&
     &n(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta&
     &)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D1 * (-sig(4)&
     & * sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) + 0.8D1 / 0.3D1&
     & * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / &
     &0.3D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(thet&
     &a) + sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + s&
     &ig(3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) *&
     & sin(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(t&
     &heta)) ** 2 / 0.3D1) * sig(6)) / 0.4D1)
      fp(6) = -rho0 * (sig(1) * cos(theta) * mu1 / rho0 * (0.4D1 * (sig(&
     &4) * cos(theta) + sig(5) * sin(theta)) * sig(6) * sig(1) * cos(the&
     &ta) - 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6&
     &) * sig(1) * sin(theta)) / 0.4D1 + (sig(2) * cos(theta) + sig(3) *&
     & sin(theta)) * mu1 / rho0 * (0.4D1 * (sig(4) * cos(theta) + sig(5)&
     & * sin(theta)) * sig(6) * (sig(2) * cos(theta) + sig(3) * sin(thet&
     &a)) + 0.4D1 * (-sig(4) * sin(theta) + sig(5) * cos(theta)) * sig(6&
     &) * (-sig(2) * sin(theta) + sig(3) * cos(theta))) / 0.4D1 + (sig(4&
     &) * cos(theta) + sig(5) * sin(theta)) * mu1 / rho0 * (-0.4D1 / 0.3&
     &D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 * cos(theta) ** 2 + 0.2D1 / 0.3D&
     &1 * (sig(2) * cos(theta) + sig(3) * sin(theta)) ** 2 + 0.2D1 / 0.3&
     &D1 * (sig(4) * cos(theta) + sig(5) * sin(theta)) ** 2 - sig(1) ** &
     &2 * sin(theta) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(3) * cos&
     &(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(theta&
     &)) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D1 * (sig(4) &
     &* cos(theta) + sig(5) * sin(theta)) ** 2 * sig(6) - 0.4D1 / 0.3D1 &
     &* (0.2D1 / 0.3D1 * sig(1) ** 2 * sin(theta) ** 2 + 0.2D1 / 0.3D1 *&
     & (-sig(2) * sin(theta) + sig(3) * cos(theta)) ** 2 + 0.2D1 / 0.3D1&
     & * (-sig(4) * sin(theta) + sig(5) * cos(theta)) ** 2 - sig(1) ** 2&
     & * cos(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) + sig(3) * sin(t&
     &heta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * sin(theta)) &
     &** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D1 * (-sig(4) * &
     &sin(theta) + sig(5) * cos(theta)) ** 2 * sig(6) + 0.8D1 / 0.3D1 * &
     &(0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 * cos(theta) ** 2 / 0.3&
     &D1 - sig(1) ** 2 * sin(theta) ** 2 / 0.3D1 - (sig(2) * cos(theta) &
     &+ sig(3) * sin(theta)) ** 2 / 0.3D1 - (-sig(2) * sin(theta) + sig(&
     &3) * cos(theta)) ** 2 / 0.3D1 - (sig(4) * cos(theta) + sig(5) * si&
     &n(theta)) ** 2 / 0.3D1 - (-sig(4) * sin(theta) + sig(5) * cos(thet&
     &a)) ** 2 / 0.3D1) * sig(6)) / 0.4D1)
            case(2)
    fp(1) = -rho0 * (sig(1) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.2D1 / 0&
     &.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 *&
     & sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) &
     &** 2 / 0.3D1) * sig(1) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3) *&
     &* 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) *&
     &* 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(1) &
     &- 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3&
     &D1 - sig(2) ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3&
     &D1 - sig(5) ** 2 / 0.3D1) * sig(1)) / 0.4D1 + sig(2) * mu1 / rho0 &
     &* (0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * &
     &sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - &
     &sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(2) + 0.4D1 * (sig&
     &(3) * sig(2) + sig(4) * sig(5)) * sig(3) - 0.4D1 / 0.3D1 * (0.2D1 &
     &/ 0.3D1 * sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 &
     &/ 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 &
     &/ 0.3D1) * sig(2) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 -&
     & sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 -&
     & sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) * sig(2)) / 0.4D1 + si&
     &g(4) * mu1 / rho0 * (0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 &
     &+ 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(&
     &3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig&
     &(4) + 0.4D1 * (sig(3) * sig(2) + sig(4) * sig(5)) * sig(5) + 0.4D1&
     & * sig(4) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3) &
     &** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) &
     &** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(4)&
     & - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.&
     &3D1 - sig(2) ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.&
     &3D1 - sig(5) ** 2 / 0.3D1) * sig(4)) / 0.4D1) + (lambda1 + 0.2D1 / 0.3D&
     &1 * mu1) * sig(1) ** 2 * sig(3) ** 2 * sig(6) ** 2 * (-sig(1) * si&
     &g(3) * sig(6) + 0.1D1) - rho0 * T0 * (sig(1) * sig(3) * sig(6)) **&
     & g * (exp(s / cv) - 0.1D1)
      fp(2) = -rho0 * (sig(3) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / &
     &0.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 &
     &* sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6)&
     & ** 2 / 0.3D1) * sig(3) + 0.4D1 * (sig(3) * sig(2) + sig(4) * sig(&
     &5)) * sig(2) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3) ** 2 + 0.2D&
     &1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D&
     &1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(3) - 0.4D1 / &
     &0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2&
     &) ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5&
     &) ** 2 / 0.3D1) * sig(3)) / 0.4D1 + sig(5) * mu1 / rho0 * (-0.4D1 &
     &/ 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) **&
     & 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) **&
     & 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(5) + 0.4D1 * (sig(3) * sig&
     &(2) + sig(4) * sig(5)) * sig(4) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 *&
     & sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 -&
     & sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) &
     &* sig(5) + 0.4D1 * sig(5) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 /&
     & 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 -&
     & sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) &
     &* sig(5)) / 0.4D1) + (lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) ** 2 * si&
     &g(3) ** 2 * sig(6) ** 2 * (-sig(1) * sig(3) * sig(6) + 0.1D1) - rh&
     &o0 * T0 * (sig(1) * sig(3) * sig(6)) ** g * (exp(s / cv) - 0.1D1)
      fp(3) = -sig(6) * mu1 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) &
     &** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 -&
     & sig(3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) &
     &* sig(6) + 0.4D1 * sig(4) ** 2 * sig(6) - 0.4D1 / 0.3D1 * (0.2D1 /&
     & 0.3D1 * sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 /&
     & 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 /&
     & 0.3D1) * sig(6) + 0.4D1 * sig(5) ** 2 * sig(6) + 0.8D1 / 0.3D1 * &
     &(0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 /&
     & 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 /&
     & 0.3D1) * sig(6)) / 0.4D1 + (lambda1 + 0.2D1 / 0.3D1 * mu1) * sig(1) **&
     & 2 * sig(3) ** 2 * sig(6) ** 2 * (-sig(1) * sig(3) * sig(6) + 0.1D&
     &1) - rho0 * T0 * (sig(1) * sig(3) * sig(6)) ** g * (exp(s / cv) - &
     &0.1D1)
      fp(4) = -rho0 * (sig(1) ** 2 * mu1 / rho0 * (sig(3) * sig(2) + sig&
     &(4) * sig(5)) + sig(2) * mu1 / rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0&
     &.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 *&
     & sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) &
     &** 2 / 0.3D1) * sig(3) + 0.4D1 * (sig(3) * sig(2) + sig(4) * sig(5&
     &)) * sig(2) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3) ** 2 + 0.2D1&
     & / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1&
     & - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(3) - 0.4D1 / 0&
     &.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2)&
     & ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5)&
     & ** 2 / 0.3D1) * sig(3)) / 0.4D1 + sig(4) * mu1 / rho0 * (-0.4D1 /&
     & 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** &
     &2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** &
     &2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(5) + 0.4D1 * (sig(3) * sig(&
     &2) + sig(4) * sig(5)) * sig(4) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * &
     &sig(3) ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - &
     &sig(2) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) *&
     & sig(5) + 0.4D1 * sig(5) * sig(6) ** 2 - 0.4D1 / 0.3D1 * (0.2D1 / &
     &0.3D1 * sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - &
     &sig(3) ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) *&
     & sig(5)) / 0.4D1)
      fp(5) = -rho0 * (sig(3) * mu1 / rho0 * (0.4D1 * sig(4) * sig(6) * &
     &sig(2) + 0.4D1 * sig(5) * sig(6) * sig(3)) / 0.4D1 + sig(5) * mu1 &
     &/ rho0 * (-0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0.2D1 / &
     &0.3D1 * sig(2) ** 2 + 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) ** 2 / &
     &0.3D1 - sig(5) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D&
     &1 * sig(4) ** 2 * sig(6) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3)&
     & ** 2 + 0.2D1 / 0.3D1 * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2)&
     & ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6&
     &) + 0.4D1 * sig(5) ** 2 * sig(6) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 &
     &* sig(6) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(3)&
     & ** 2 / 0.3D1 - sig(4) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) * sig(6&
     &)) / 0.4D1)
      fp(6) = -rho0 * (sig(1) ** 2 * mu1 / rho0 * sig(4) * sig(6) + sig(&
     &2) * mu1 / rho0 * (0.4D1 * sig(4) * sig(6) * sig(2) + 0.4D1 * sig(&
     &5) * sig(6) * sig(3)) / 0.4D1 + sig(4) * mu1 / rho0 * (-0.4D1 / 0.&
     &3D1 * (0.2D1 / 0.3D1 * sig(1) ** 2 + 0.2D1 / 0.3D1 * sig(2) ** 2 +&
     & 0.2D1 / 0.3D1 * sig(4) ** 2 - sig(3) ** 2 / 0.3D1 - sig(5) ** 2 /&
     & 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D1 * sig(4) ** 2 * sig&
     &(6) - 0.4D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(3) ** 2 + 0.2D1 / 0.3D1&
     & * sig(5) ** 2 - sig(1) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(4&
     &) ** 2 / 0.3D1 - sig(6) ** 2 / 0.3D1) * sig(6) + 0.4D1 * sig(5) **&
     & 2 * sig(6) + 0.8D1 / 0.3D1 * (0.2D1 / 0.3D1 * sig(6) ** 2 - sig(1&
     &) ** 2 / 0.3D1 - sig(2) ** 2 / 0.3D1 - sig(3) ** 2 / 0.3D1 - sig(4&
     &) ** 2 / 0.3D1 - sig(5) ** 2 / 0.3D1) * sig(6)) / 0.4D1)            
                case default
                print *, 'EOS mode in A2sigmaComputeFSGEOS not implemented',EOS_mode
                end select
    end subroutine A2sigmaComputeFSGEOS
    ! *******************************************************************************************************
    PURE SUBROUTINE compute_l_m_dedc(lam, mu, dEdc, c, devG, detA, ODE)
        use expintegrator_type, only : T_ode_parameters,nvarode,ncc,nccc
        ! Elasticity constant homogeneization and computation of the derivative of 
        ! energy with respect to the damage parameter.
        REAL(8),                 INTENT(OUT) :: lam, mu, dEdc
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: devG
        TYPE(T_ode_parameters),  INTENT(IN)  :: ODE
        REAL(8),                 INTENT(IN)  :: c, detA
        REAL(8)                              :: mus, Ks, K, dM, dK
        REAL(8), PARAMETER                   :: zero = 1.0d-100
        mus = c*ODE%mu1 + (1.0d0 - c)*ODE%mu2
        Ks  = c*ODE%K1  + (1.0d0 - c)*ODE%K2
        mu  = ODE%mu1*ODE%mu2/(mus + zero)
        K   = ODE%K1*ODE%K2/(Ks + zero)
        lam = K - 2.0d0/3.0d0*mu
        dM  = (ODE%mu1 - ODE%mu2)*mu/mus
        dK  = (ODE%K1 - ODE%K2)*K/Ks
        dEdc = -(0.5d0*dK*(1.0d0 - detA)**2 + 0.25d0*dM*sum(devG**2))/ODE%rho0 
    END SUBROUTINE compute_l_m_dedc
    
    PURE SUBROUTINE compute_l_m_mix(lam, mu,lam1,mu1,lam12,mu12,c)
        IMPLICIT NONE
        REAL(8), INTENT(IN)  :: lam1,mu1,lam12,mu12,c
        REAL(8), INTENT(OUT) :: lam, mu
        REAL(8)             :: K1,K2,Ks,mus,K,lam2,mu2
        REAL(8), PARAMETER                   :: zero = 1.0d-100
        
        mu2     = mu1*mu12   
        lam2 = lam1+lam12*(mu1-mu2) 
        
        K1     = lam1 + 2.0d0*mu1/3.0d0
        K2     = lam2 + 2.0d0*mu2/3.0d0
        mus = c*mu1 + (1.0d0 - c)*mu2
        Ks  = c*K1  + (1.0d0 - c)*K2
        mu  = mu1*mu2/(mus + zero)
        K   = K1*K2/(Ks + zero)
        lam = K - 2.0d0/3.0d0*mu
    END SUBROUTINE compute_l_m_mix
    ! *******************************************************************************************************
    PURE SUBROUTINE compute_hydrodynamic_pressure(hydro_pressure, detA, lam, mu,S,gamma,cv, p0,rho,EOS,T0)
        REAL(8), INTENT(OUT) :: hydro_pressure
        REAL(8), INTENT(IN)  :: detA, lam, mu, S,gamma,cv, p0,T0,rho
        integer, intent(in)  :: EOS
        REAL(8)              :: bulk_modulus
        bulk_modulus = lam + 2.0d0/3.0d0*mu
        !hydro_pressure = -bulk_modulus*detA**2*(1.0d0 - detA)   ! Old equation of state
        select case(EOS)
            case(0)
                hydro_pressure = -bulk_modulus*detA**2*(1.0d0 - detA)
            case(1)
                hydro_pressure = bulk_modulus/gamma*(detA)**gamma*EXP(S/cv) - bulk_modulus/gamma + p0
            case(2)
                hydro_pressure = -bulk_modulus*detA**2*(1.0d0 - detA)+  rho*detA**gamma*T0*(EXP(S/cv)-1.)
            case default
                !print *,'EOS ',EOS, ' not implemented in compute_hydrodynamic_pressure' ! For pure functions this is not allowed =)
        end select
    END SUBROUTINE compute_hydrodynamic_pressure
    ! *******************************************************************************************************
    PURE SUBROUTINE compute_temperature(temperature, detA, lam, mu,rho,gamma,cv, p0,S,EOS,T0)
        REAL(8), INTENT(OUT) :: temperature
        REAL(8), INTENT(IN)  :: detA, lam, mu, rho,gamma,cv, p0, S,T0
        INTEGER, INTENT(IN)  :: EOS
        REAL(8)              :: bulk_modulus
        bulk_modulus = lam + 2.0d0/3.0d0*mu
        select case(EOS)
            case(0)
                temperature = 0.
            case(1)
                temperature = 1.0/cv/(gamma-1.0)/gamma*bulk_modulus/rho*detA**gamma*EXP(S/cv)
            case(2)
                temperature = T0*detA**gamma*EXP(S/cv)
        end select
        !temperature = 0.0
    END SUBROUTINE compute_temperature
    ! *******************************************************************************************************
    PURE SUBROUTINE compute_s_of_T(S,temperature, detA, lam, mu,rho,gamma,cv, p0,EOS,T0)
        REAL(8), INTENT(OUT) :: S
        REAL(8), INTENT(IN)  :: detA, lam, mu, rho,gamma,cv, p0,T0,temperature
        INTEGER, INTENT(IN)  :: EOS
        REAL(8)              :: bulk_modulus
        bulk_modulus = lam + 2.0d0/3.0d0*mu
        select case(EOS)
            case(0)
                S = 0.
            case(1)
                S=log(temperature*cv*(gamma-1.0)*gamma*bulk_modulus*rho/detA**gamma)*cv
            case(2)
                S=log(temperature/(T0*detA**gamma))*cv
        end select
    END SUBROUTINE compute_s_of_T
    ! *******************************************************************************************************
    PURE SUBROUTINE compute_coulomb(Yeq2, Ydev,Yp,  lam, mu, G, devG, detA, ODE)
        use expintegrator_type, only : T_ode_parameters,nvarode,ncc,nccc
        USE expintegrator_linalg, only: compute_matrix_exponential, solve, determinant_lu, build_jacobian_matrix_mask, &
                            build_reduced_jacobian_matrix, compute_first_order_linear_ode_solution, select_solution_from_mask, halven_rank42_rows, &
                            flatten_rows, build_by_rows, cross_product, diadic_product, trace
        USE, INTRINSIC :: ieee_arithmetic
        ! Computation of total stress
        REAL(8),                 INTENT(OUT) :: Yeq2, Yp, Ydev
        REAL(8),                 INTENT(IN)  :: lam, mu
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: G, devG
        REAL(8),                 INTENT(IN)  :: detA
        TYPE(T_ode_parameters),  INTENT(IN)  :: ODE
        REAL(8), DIMENSION(3,3)              :: sigma_visc, sigma, dev_sigma
        REAL(8)                              :: hydro_pressure, total_pressure, J2_dev,Ypos,Yneg,xi
        REAL(8), DIMENSION(3,3), PARAMETER   :: ID = reshape(&
            [1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0], [3,3])
        sigma_visc = -detA*mu*matmul(G, devG) ! note: non trace free, but without hydrodynamic pressure, changed!
        CALL compute_hydrodynamic_pressure(hydro_pressure, detA, lam, mu,ODE%S,ODE%gamma,ODE%cv, ODE%p0,ODE%rho0, ODE%EOS, ODE%T0)
        sigma = sigma_visc + hydro_pressure*ID
        total_pressure = trace(sigma)/3.0d0 ! note: total pressure = hydro_pressure + tr(sigma_visc)/3
        dev_sigma = sigma - total_pressure*ID
        J2_dev = 0.5d0*trace(matmul(dev_sigma, dev_sigma))
        Ydev = sqrt(3.0d0*J2_dev)
        Yp = total_pressure
	Yeq2=Ydev
	end subroutine compute_coulomb

    PURE SUBROUTINE compute_total_stress(Y, Ydev,Yp,  lam, mu, G, devG, detA, ODE)
        use expintegrator_type, only : T_ode_parameters,nvarode,ncc,nccc
        USE expintegrator_linalg, only: compute_matrix_exponential, solve, determinant_lu, build_jacobian_matrix_mask, &
                            build_reduced_jacobian_matrix, compute_first_order_linear_ode_solution, select_solution_from_mask, halven_rank42_rows, &
                            flatten_rows, build_by_rows, cross_product, diadic_product, trace
        USE, INTRINSIC :: ieee_arithmetic
        ! Computation of total stress
        REAL(8),                 INTENT(OUT) :: Y, Yp, Ydev
        REAL(8),                 INTENT(IN)  :: lam, mu
        REAL(8), DIMENSION(3,3), INTENT(IN)  :: G, devG
        REAL(8),                 INTENT(IN)  :: detA
        TYPE(T_ode_parameters),  INTENT(IN)  :: ODE
        REAL(8), DIMENSION(3,3)              :: sigma_visc, sigma, dev_sigma
        REAL(8)                              :: hydro_pressure, total_pressure, J2_dev,Ypos,Yneg,xi
        REAL(8), DIMENSION(3,3), PARAMETER   :: ID = reshape(&
            [1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0], [3,3])
        sigma_visc = -detA*mu*matmul(G, devG) ! note: non trace free, but without hydrodynamic pressure, changed!
        CALL compute_hydrodynamic_pressure(hydro_pressure, detA, lam, mu,ODE%S,ODE%gamma,ODE%cv, ODE%p0,ODE%rho0, ODE%EOS, ODE%T0)
        sigma = sigma_visc + hydro_pressure*ID
        total_pressure = trace(sigma)/3.0d0 ! note: total pressure = hydro_pressure + tr(sigma_visc)/3
        dev_sigma = sigma - total_pressure*ID
        J2_dev = 0.5d0*trace(matmul(dev_sigma, dev_sigma))
        Ydev = sqrt(3.0d0*J2_dev)
        Yp = total_pressure
        IF (ODE%Yeq_mode .eq. 1) THEN
            Y = ODE%Yeq_A*Ydev + ODE%Yeq_B*abs(Yp) ! default 0.25, 0.8
            Y = Y*ODE%alpha**2
        ELSE IF (ODE%Yeq_mode .eq. 2) THEN
            ! default ODE%Yeq_A = 0.5
            ! default ODE%Yeq_s0 = 3.2e+9
            ! ODE%Yeq_smoothing = 1.0d-4
            Y       = ODE%Yeq_B*abs(Yp) + ODE%Yeq_A*(Ydev - ODE%Yeq_s0)*0.5d0*(1.0d0 + erf((Ydev - ODE%Yeq_s0)/ODE%Yeq_smoothing))
            Y = Y*ODE%alpha**2
        ELSE IF (ODE%Yeq_mode .eq. 3) THEN
            !Y = ODE%Yeq_C !!! posto per materiale debole a trazione
            Y       = ODE%Yeq_A*Ydev+ODE%Yeq_B*ABS(Yp)*0.5d0*(1.0d0 + erf((-Yp + ODE%Yeq_s0)/ODE%Yeq_smoothing))
            Y = Y*ODE%alpha**2
        ELSE IF (ODE%Yeq_mode .eq. 4) THEN
            ! ----------------------------------------------------------------------------------------------------------------
            ! Generalized mode that contains all the possible failure configurations. 
 
            xi      = 0.5d0*(1.0d0 + erf((Ydev - ODE%Yeq_s0)/min(1.e+5,ODE%Y0*1.e-3))) ! Smooth coefficient for the shear contribution (take 1e-3 times the illness stress but never over 1e+5
            Ypos    = ODE%Yeq_B*abs(Yp)+xi*ODE%Yeq_A*(Ydev - ODE%Yeq_s0)
            Yneg    = ODE%Yeq_C*abs(Yp)+ODE%Yeq_A*(Ydev - ODE%Yeq_s0)
            xi      = 0.5d0*(1.0d0 + erf(Yp/min(1.e+5,ODE%Y0*1.e-3)))                ! Smooth the positive/negative pressure corner
			!if(Yp>0) then
			!	xi=1.
			!else
			!	xi=0.
			!end if
			!print *, 'Here', ODE%Yeq_C
            Y       = xi*Ypos+(1.-xi)*Yneg
            Y = Y*ODE%alpha
        ELSE IF (ODE%Yeq_mode .eq. 5) THEN
            xi= 0.5d0*(1.0d0 + erf((Ydev - ODE%Yeq_s0)/min(1.e+5,ODE%Y0*1.e-3)))
            Y = ODE%Yeq_A*Ydev + ODE%Yeq_B*abs(Yp) ! default 0.25, 0.8
            Y = Y*ODE%alpha**2
        ELSE
            Y = -1.0d300
        END IF
			!Y=max(0., Y)
    END SUBROUTINE compute_total_stress
    ! *******************************************************************************************************
    subroutine StaticLimiterEQ99(dmpresult,lxb,ldx)
        use MainVariables, only : ICType
        implicit none
        logical :: dmpresult
        real    :: lxb(3),ldx(3)
        
        SELECT CASE(ICType)
            CASE('SSCRACK')
                if(abs(lxb(1))<10000.0+ldx(1) .and. abs(lxb(2))<ldx(2)) then
                    dmpresult = .FALSE.
                end if
            CASE('StiffInclusion')
                ! Manually activate the limiter on the strong material jumps
                if(     ( (abs(lxb(1)-0.5)<(ldx(1)) .or. abs(lxb(1)+0.5)<(ldx(1)))  .and. (lxb(2)>-0.1-(ldx(2)) .and. lxb(2)<0.1+(ldx(2)))) ) then
                     dmpresult = .FALSE.
                end if
                if(     ( (abs(lxb(2)-0.1)<(ldx(2)) .or. abs(lxb(2)+0.1)<(ldx(2)))  .and. (lxb(1)>-0.5-(ldx(1)) .and. lxb(1)<0.5+(ldx(1)))) ) then
                     dmpresult = .FALSE.
                end if
            CASE('TPV3','NMGPRDB')
                if(abs(lxb(2))<ldx(2)) then
                    dmpresult = .FALSE.
                end if
                !if(abs(lxb(2)) < 700) then
                !    dmpresult = .FALSE.   
                !end if
            CASE('NLRUPTURE')
                if(abs(lxb(2))<ldx(2)) then
                    dmpresult = .FALSE.
                end if
        END SELECT
        
    end subroutine StaticLimiterEQ99
    !*********************************************************************************
    function SmoothInterface(r,ICsig,epsilon,smooth_order_in)
        implicit none
        real    :: r
        real    :: SmoothInterface,smooth_order
        real    :: eta,ICsig,xi 
        real, optional :: epsilon
        real, optional :: smooth_order_in
        
        if(.not. present(epsilon)) then
            epsilon=1.e-9    
        end if
        if(.not. present(smooth_order_in)) then
            smooth_order=4   
        end if
        smooth_order=smooth_order_in

        eta=0.0
        smooth_order=4.0
        ! =============== WAY 1 ================================
        if(r>(1+eta)*ICsig) then
            xi=1    
        elseif(r<-(1-eta)*ICsig) then
            xi=0
        else
            xi = (r+(1-eta)*ICsig)/(2.0*ICsig) 
        end if
        ! =============== WAY 2 ================================
        SmoothInterface  = (1.0)*(1-xi)**smooth_order + (0.0+epsilon)*xi**smooth_order 
        if(abs(smooth_order) .le. 1.e-3) then
            xi = 0.5+0.5*ERF(r/(2*ICsig))  
            SmoothInterface  = (1.0-epsilon)*(1-xi) + (0.0+epsilon)*xi
        end if
    end function SmoothInterface
	
SUBROUTINE LinSolve(N,A,b,x)
  IMPLICIT NONE
  INTEGER       :: N
  REAL          :: A(N,N), b(N), x(N)  
  !
  INTEGER       :: i,j,ml(1) 
  REAL          :: temp(N+1),piv 
  REAL          :: C(N+1,N)  
  !
  C(1:N,:)     = TRANSPOSE(A(:,:))
  C(N+1,:)     = b 
  !    
  ! Forward elimination with column pivoting 
  ! 
  DO i = 1, N
     ! If pivot element is zero, then swap rows 
     ml = MAXLOC(ABS(C(i,i:N))) 
     j = i - 1 + ml(1) 
     temp   = C(:,j) 
     C(:,j) = C(:,i)
     C(:,i) = temp      
     IF(C(i,i).EQ.0.) THEN
        !PRINT *, 'ERROR. Matrix is singular!'
        !DO j = 1, N
        !   PRINT *, A(j,:) 
        !ENDDO
        !STOP
x = 0. 
RETURN 
     ENDIF
     piv    = 1./C(i,i)
     C(:,i) = C(:,i)*piv 
     DO j = i+1, N 
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  ! Back substitution
  !
  DO i = N,1,-1   
     DO j = i-1,1,-1
        C(:,j) = C(:,j) - C(i,j)*C(:,i)
     ENDDO
  ENDDO
  !
  x = C(N+1,:)  
  !
END SUBROUTINE LinSolve
	
	
    end module SpecificVarEqn99
#endif
