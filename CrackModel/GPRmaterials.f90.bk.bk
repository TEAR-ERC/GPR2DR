#ifndef HEADER_GPRMATERIALS
#define HEADER_GPRMATERIALS
!#include "MainVariables.f90"
#if defined(EQNTYPEC99) || defined(EQNTYPEB99) || defined(EQNTYPE99) || defined(EQNTYPED99)
! This module contains the parameters of several materials as well as several tools in order to check the material response
MODULE GPRmaterials
    IMPLICIT NONE
    PRIVATE
    
    PUBLIC :: AssignMaterialProperties
    PUBLIC :: AssignMaterialPropertiesMix
    PUBLIC :: getLameCoefficients
    PUBLIC :: getLameCoefficientsMix

CONTAINS
    FUNCTION  getLameCoefficients(MATERIAL)  
        IMPLICIT NONE
        CHARACTER(LEN=*), INTENT(IN)   :: MATERIAL         ! Identify the material
        REAL             :: getLameCoefficients(3)
        ! Local variables
        REAL                           :: ChWSpeed(2)      ! The two characteristic wave speeds
        REAL                           :: LamCoeff(3)      ! (lambda, mu, rho)
        REAL                           :: DamCoeff(2)      ! (mu{2->1} lambda{2->1})
        REAL                           :: RuptureCoeff(8)  ! (K,Y0,Y1,alpha1,beta1,alpha2,beta2)
        REAL                           :: StressCoeff(4)   ! YEQ A,B,C and s0


        call GetMaterialParameters(ChWSpeed,LamCoeff,DamCoeff,RuptureCoeff,StressCoeff, MATERIAL) 

       getLameCoefficients  = LamCoeff(:)


    END FUNCTION getLameCoefficients
    
    FUNCTION  getLameCoefficientsMix(MATERIALS,mixXI,nMat)  
        IMPLICIT NONE
        INTEGER, INTENT(IN)            :: nMAT
        REAL            , INTENT(IN)   :: mixXI(nMAT)          ! State vector to be assigned
        CHARACTER(LEN=*), INTENT(IN)   :: MATERIALS(nMAT)  ! Identify the material
        REAL             :: getLameCoefficientsMix(3)
        REAL                           :: ChWSpeed(2,nMAT)      ! The two characteristic wave speeds
        REAL                           :: LamCoeff(3,nMAT)      ! (lambda, mu, rho)
        REAL                           :: DamCoeff(2,nMAT)      ! (mu{2->1} lambda{2->1})
        REAL                           :: RuptureCoeff(8,nMAT)  ! (K,Y0,Y1,alpha1,beta1,alpha2,beta2)
        REAL                           :: StressCoeff(4,nMAT)   ! YEQ A,B,C and s0
        ! Local variables:
        REAL                           :: locXI(nMAT)
        REAL                           :: XITOT,XImean(nMAT)
        INTEGER                        :: iMAT
        LOGICAL                        :: FLAG
        REAL                           :: LChWSpeed(2)     
        REAL                           :: LLamCoeff(3)     
        REAL                           :: LDamCoeff(2)     
        REAL                           :: LRuptureCoeff(8) 
        REAL                           :: LStressCoeff(4)  

        do iMAT=1,nMAT
            locXI(iMAT)=max(0., min(1.,mixXI(iMAT)))
            !if(locXI(iMAT)>0.5) then
            !    locXI(iMAT)=1.0   
            !else
            !    locXI(iMAT)=0.0  
            !endif
        END DO
		if(locXI(1)>0.5) then
			locXI(1)=1.0
			locXI(2)=1.0-locXI(1)
		else
		    locXI(2)=1.0
			locXI(1)=1.0-locXI(2)
		end if
        XITOT=sum(locXI)
        XImean=locXI/XITOT
        
        do iMAT=1,nMAT
            call GetMaterialParameters(ChWSpeed(:,iMAT),LamCoeff(:,iMAT),DamCoeff(:,iMAT),RuptureCoeff(:,iMAT),StressCoeff(:,iMAT), MATERIALS(iMAT)) 
        end do
        
        LChWSpeed     = 0.
        LLamCoeff     = 0.
        
        do iMAT=1,nMAT
          LChWSpeed     = LChWSpeed     + XImean(iMat)*ChWSpeed(:,iMat)     
          LLamCoeff     = LLamCoeff     + XImean(iMat)*LamCoeff(:,iMat)       
        end do

        getLameCoefficientsMix  = LLamCoeff(:)

    END FUNCTION getLameCoefficientsMix
    
    PURE SUBROUTINE AssignMaterialPropertiesMix(Q,MATERIALS,mixXI,nMAT,RUPTUREONLY)   
        USE MainVariables, only : nVar
        IMPLICIT NONE   
        INTEGER, INTENT(IN)            :: nMAT
#ifdef EQNTYPED99
        REAL            , INTENT(INOUT):: Q(14)          ! For D99 it contains just the damage parameters
#else
        REAL            , INTENT(INOUT):: Q(nVar)          ! State vector to be assigned       
#endif
        REAL            , INTENT(IN)   :: mixXI(nMAT)          ! State vector to be assigned
        CHARACTER(LEN=*), INTENT(IN)   :: MATERIALS(nMAT)  ! Identify the material
        LOGICAL, INTENT(IN), OPTIONAL  :: RUPTUREONLY      ! Assign only the rupture properties
        REAL                           :: ChWSpeed(2,nMAT)      ! The two characteristic wave speeds
        REAL                           :: LamCoeff(3,nMAT)      ! (lambda, mu, rho)
        REAL                           :: DamCoeff(2,nMAT)      ! (mu{2->1} lambda{2->1})
        REAL                           :: RuptureCoeff(8,nMAT)  ! (K,Y0,Y1,alpha1,beta1,alpha2,beta2)
        REAL                           :: StressCoeff(4,nMAT)   ! YEQ A,B,C and s0
        ! Local variables:
        REAL                           :: locXI(nMAT)
        REAL                           :: XITOT,XImean(nMAT)
        INTEGER                        :: iMAT
        LOGICAL                        :: FLAG
        REAL                           :: LChWSpeed(2)     
        REAL                           :: LLamCoeff(3)     
        REAL                           :: LDamCoeff(2)     
        REAL                           :: LRuptureCoeff(8) 
        REAL                           :: LStressCoeff(4)  
        
        IF(.NOT. PRESENT(RUPTUREONLY)) THEN
            FLAG = .FALSE.
        ELSE
            FLAG = RUPTUREONLY
        END IF
        do iMAT=1,nMAT
            locXI(iMAT)=max(0., min(1.,mixXI(iMAT)))
            !if(locXI(iMAT)>0.5) then
            !    locXI(iMAT)=1.0   
            !else
            !    locXI(iMAT)=0.0  
            !endif
        END DO
		if(locXI(1)>0.5) then
			locXI(1)=1.0
			locXI(2)=1.0-locXI(1)
			else
		    locXI(2)=1.0
			locXI(1)=1.0-locXI(2)
		end if
        XITOT=sum(locXI)
        XImean=locXI/XITOT
        
        do iMAT=1,nMAT
            call GetMaterialParameters(ChWSpeed(:,iMAT),LamCoeff(:,iMAT),DamCoeff(:,iMAT),RuptureCoeff(:,iMAT),StressCoeff(:,iMAT), MATERIALS(iMAT)) 
        end do
        
        LChWSpeed     = 0.
        LLamCoeff     = 0.
        LDamCoeff     = 0.
        LRuptureCoeff = 0.
        LStressCoeff  = 0.
        
        do iMAT=1,nMAT
          LChWSpeed     = LChWSpeed     + XImean(iMat)*ChWSpeed(:,iMat)     
          LLamCoeff     = LLamCoeff     + XImean(iMat)*LamCoeff(:,iMat)     
          LDamCoeff     = LDamCoeff     + XImean(iMat)*DamCoeff(:,iMat)     
          LRuptureCoeff = LRuptureCoeff + XImean(iMat)*RuptureCoeff(:,iMat) 
          LStressCoeff  = LStressCoeff  + XImean(iMat)*StressCoeff(:,iMat)     
        end do
#if defined(EQNTYPEC99)
        IF(.NOT. FLAG) THEN
            Q(1)     = LLamCoeff(1)
            Q(20)    = LLamCoeff(2)
            Q(19)    = LLamCoeff(3)
        END IF
        !        
        Q(22)    = LDamCoeff(1)
        Q(23)    = LDamCoeff(2)
        !
        Q(24:31) = LRuptureCoeff(1:8)
        Q(32:35) = LStressCoeff(1:4)
#elif defined(EQNTYPED99)
        Q(1:2)    = LDamCoeff(1:2)
        Q(3:10) = LRuptureCoeff(1:8)
        Q(11:14) = LStressCoeff(1:4)
#else
        !print *, 'Warning, function "AssignMaterialProperties" not implemented for the current PDE system!'
        stop
#endif
    END SUBROUTINE AssignMaterialPropertiesMix
  
	
    SUBROUTINE AssignMaterialProperties(Q,MATERIAL,RUPTUREONLY)  
        USE MainVariables, only : nVar
        IMPLICIT NONE
        REAL            , INTENT(INOUT):: Q(nVar)          ! State vector to be assigned
        CHARACTER(LEN=*), INTENT(IN)   :: MATERIAL         ! Identify the material
        LOGICAL, INTENT(IN), OPTIONAL  :: RUPTUREONLY      ! Assign only the rupture properties
        ! Local variables
        REAL                           :: ChWSpeed(2)      ! The two characteristic wave speeds
        REAL                           :: LamCoeff(3)      ! (lambda, mu, rho)
        REAL                           :: DamCoeff(2)      ! (mu{2->1} lambda{2->1})
        REAL                           :: RuptureCoeff(8)  ! (K,Y0,Y1,alpha1,beta1,alpha2,beta2)
        REAL                           :: StressCoeff(4)   ! YEQ A,B,C and s0
        LOGICAL                        :: FLAG

        IF(.NOT. PRESENT(RUPTUREONLY)) THEN
            FLAG = .FALSE.
        ELSE
            FLAG = RUPTUREONLY
        END IF
        
        call GetMaterialParameters(ChWSpeed,LamCoeff,DamCoeff,RuptureCoeff,StressCoeff, MATERIAL) 
#if defined(EQNTYPEC99)
        IF(.NOT. FLAG) THEN
            Q(1)     = LamCoeff(1)
            Q(20)    = LamCoeff(2)
            Q(19)    = LamCoeff(3)
        END IF
        !        
        Q(22)    = DamCoeff(1)
        Q(23)    = DamCoeff(2)
        !
        Q(24:31) = RuptureCoeff(1:8)
        Q(32:35) = StressCoeff(1:4)
#elif defined(EQNTYPED99)
        ! Just assign the Lamè properties
        Q(1)     = LamCoeff(1)
        Q(20)    = LamCoeff(2)
        Q(19)    = LamCoeff(3)
#else
        print *, 'Warning, function "AssignMaterialProperties" not implemented for the current PDE system!'
        stop
#endif
        
    END SUBROUTINE AssignMaterialProperties

    PURE SUBROUTINE GetMaterialParameters(ChWSpeed,LamCoeff,DamCoeff,RuptureCoeff,StressCoeff, MATERIAL)  
        IMPLICIT NONE
        REAL, INTENT(OUT)              :: ChWSpeed(2)      ! The two characteristic wave speeds
        REAL, INTENT(OUT)              :: LamCoeff(3)      ! (lambda, mu, rho)
        REAL, INTENT(OUT)              :: DamCoeff(2)      ! ( mu{2->1}, lambda{2->1})
        REAL, INTENT(OUT)              :: RuptureCoeff(8)  ! (K,Y0,Y1,alpha1,beta1,alpha2,beta2)
        REAL, INTENT(OUT)              :: StressCoeff(4)   ! YEQ A,B,C and s0
        CHARACTER(LEN=*), INTENT(IN)   :: MATERIAL

        !!!! ----------------------------------------------------------------------------------------------------------------!!!!
        !!!! Generalized mode that contains all the possible failure configurations (Yeq_mode=4).                            !!!!
        !!!! How to set the StressCoeff vector:                                                                              !!!!
        !!!! In the quasi-static test by defining the quasi-static critical stress by Y0                                     !!!!
        !!!! we set the parameters according to the type of the material, in particular:                                     !!!!
        !!!! 1) Shear based rupture criteria                                                                                 !!!!
        !!!!       > s0=0                                                                                                    !!!!
        !!!!       > Yb=Yc=0                                                                                                 !!!!
        !!!!       > Critical point (Ya): Ysh=Y0/Ya, Yp=infinity                                                             !!!!
        !!!! 2) Pressure based rupture criteria                                                                              !!!!
        !!!!       > s0=0                                                                                                    !!!!
        !!!!       > Ya=0                                                                                                    !!!!
        !!!!       > Critical point (Yb, Yc): Ysh=infinity Yp=Y0/Yb, Yp=-Y0/Yc                                               !!!!
        !!!! 3) Linear Pressure/Shear criteria                                                                               !!!!
        !!!!       > s0=0                                                                                                    !!!!
        !!!!       > Yc=Yb                                                                                                   !!!!
        !!!!       > Critical point (Ya, Yb): Ycr=Ya*Ysh+Yb*Yp                                                               !!!!
        !!!! 4) General shear/pressure and traction criteria                                                                 !!!!
        !!!!       > Critical point (Ya, Yb, Yc, s0): Line connecting (Y0/Yb,0), (Y0/Yb,s0), (0,Y0/Ya+s0), ((Y0-Ya*s0)/Yc,0) !!!!
        !!!! ----------------------------------------------------------------------------------------------------------------!!!!
        
        SELECT CASE(MATERIAL)   ! Select the material type
		CASE('ROCKTEST')
			! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 2670.0                                              ! RHO        !
            ChWSpeed    = (/ 3463.999983, 6000.000000 /)                      ! CS, CL     !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            !DamCoeff    = (/    1./1.16666   ,   0.666666   /)          ! MU12,LAM12               !
			DamCoeff    = (/    1./1.16666  ,   0.666666   /)          ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 1.0           ! K                                            !
            RuptureCoeff(2) = 1.8e8        ! Y0                                           !
            RuptureCoeff(3) = 0.01e12        ! Y1                                           !
            RuptureCoeff(4) = 32.5          ! aexp                                         !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 0.0           ! beta1                                        !
            RuptureCoeff(7) = 36.25         ! alpha2                                       !
            RuptureCoeff(8) = 1.e-5         ! beta2                                        !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            ! Shear based criteria                                                         !
            StressCoeff(1)  = 1.5           ! Yeq_A                                        !
            StressCoeff(2)  = -0.2           ! Yeq_B                                        !
            StressCoeff(3)  = 4.0           ! Yeq_C                                        !
            StressCoeff(4)  = 0.0           ! s0                                           !
			!StressCoeff(1)  = 0.3           ! Yeq_A                                        !
            !StressCoeff(2)  = 0.0           ! Yeq_B                                        !
            !StressCoeff(3)  = 4.7           ! Yeq_C                                        !
            !StressCoeff(4)  = 0.0           ! s0    
            ! -----------------------------------------------------------------------------!
        CASE('ROCK1')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 2670.0                                              ! RHO        !
            ChWSpeed    = (/ 3463.999983, 6000.000000 /)                      ! CS, CL     !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1./1.16666   ,   0.666666   /)          ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 1.0           ! K                                            !
            RuptureCoeff(2) = 1.8e8        ! Y0                                           !
            RuptureCoeff(3) = 0.01e9        ! Y1                                           !
            RuptureCoeff(4) = 32.5          ! aexp                                         !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 0.0           ! beta1                                        !
            RuptureCoeff(7) = 36.25         ! alpha2                                       !
            RuptureCoeff(8) = 1.e-6         ! beta2                                        !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            ! Shear based criteria                                                         !
            StressCoeff(1)  = 1.0           ! Yeq_A                                        !
            StressCoeff(2)  = 0.0           ! Yeq_B                                        !
            StressCoeff(3)  = 0.0           ! Yeq_C                                        !
            StressCoeff(4)  = 0.0           ! s0                                           !
            ! -----------------------------------------------------------------------------!
        CASE('ROCK3')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 2670.0                                              ! RHO        !
            ChWSpeed    = (/ 3463.999983, 6000.000000 /)                      ! CS, CL     !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1./1.16666   ,   0.666666   /)          ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 0.2           ! K                                            !
            RuptureCoeff(2) = 2.4e8        ! Y0                                           !
            RuptureCoeff(3) = 0.01e9        ! Y1                                           !
            RuptureCoeff(4) = 32.5          ! aexp                                         !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 0.0           ! beta1                                        !
            RuptureCoeff(7) = 36.25         ! alpha2                                       !
            RuptureCoeff(8) = 1.e-6         ! beta2                                        !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            ! Shear based criteria                                                         !
            StressCoeff(1)  = 1.0           ! Yeq_A                                        !
            StressCoeff(2)  = 0.0           ! Yeq_B                                        !
            StressCoeff(3)  = 0.0           ! Yeq_C                                        !
            StressCoeff(4)  = 0.0           ! s0                                           !
            ! -----------------------------------------------------------------------------!
        CASE('ROCK2')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 1440.0                                              ! RHO        !
            ChWSpeed    = (/ 2074.7, 3424.29 /)                      ! CS, CL              !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1.E-3    ,   0.0  /)          ! MU12,LAM12                 !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 1.0           ! K                                            !
            RuptureCoeff(2) = 5.5e5        ! Y0                                            !
            RuptureCoeff(3) = 0.01e5         ! Y1                                          !
            RuptureCoeff(4) = 52.5         ! aexp                                          !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 0.0           ! beta1                                        !
            RuptureCoeff(7) = 36.25         ! alpha2                                       !
            RuptureCoeff(8) = 1.e-6         ! beta2                                        !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            StressCoeff(1)  = 1.0           ! Yeq_A                                        !
            StressCoeff(2)  = 8.0           ! Yeq_B                                        !
            StressCoeff(3)  = 0.0           ! Yeq_C                                        !
            StressCoeff(4)  = 0.0           ! s0                                           !
            ! -----------------------------------------------------------------------------! 
        CASE('ROCK4')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 2670.0                                              ! RHO        !
            ChWSpeed    = (/ 3463.999983, 6000.000000 /)                      ! CS, CL     !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1./1.16666   ,   0.666666   /)          ! MU12,LAM12               !
            DamCoeff    = (/    0.001   ,   0.666666   /)          ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 1.0           ! K                                            !
            RuptureCoeff(2) = 9e6        ! Y0                                           !
            RuptureCoeff(3) = 0.01e9        ! Y1                                           !
            RuptureCoeff(4) = 52.5          ! aexp                                         !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 0.0           ! beta1                                        !
            RuptureCoeff(7) = 36.25         ! alpha2                                       !
            RuptureCoeff(8) = 1.e-6         ! beta2                                        !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            ! Shear based criteria                                                         !
            StressCoeff(1)  = 1.0           ! Yeq_A                                        !
            StressCoeff(2)  = 0.2           ! Yeq_B                                        !
            StressCoeff(3)  = 0.0           ! Yeq_C                                        !
            StressCoeff(4)  = 0.0           ! s0                                           !
            ! -----------------------------------------------------------------------------!
        CASE('ROCK5')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 2670.0                                              ! RHO        !
            ChWSpeed    = (/ 3463.999983, 6000.000000 /)                      ! CS, CL     !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1./1.16666   ,   0.666666   /)          ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 10.0           ! K                                            !
            RuptureCoeff(2) = 1.75e8        ! Y0                                           !
            RuptureCoeff(3) = 0.01e9        ! Y1                                           !
            RuptureCoeff(4) = 42.5          ! aexp                                         !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 0.0           ! beta1                                        !
            RuptureCoeff(7) = 36.25         ! alpha2                                       !
            RuptureCoeff(8) = 1.e-6         ! beta2                                        !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            ! Shear based criteria                                                         !
            StressCoeff(1)  = 1.0           ! Yeq_A                                        !
            StressCoeff(2)  = 0.0           ! Yeq_B                                        !
            StressCoeff(3)  = 0.0           ! Yeq_C                                        !
            StressCoeff(4)  = 1.e+30           ! s0                                           !
            ! -----------------------------------------------------------------------------!
        CASE('ROCKP')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 2670.0                                              ! RHO        !
            ChWSpeed    = (/ 3463.999983, 6000.000000 /)                      ! CS, CL     !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1./1.16666   ,   0.666666   /)          ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 10.0           ! K                                            !
            RuptureCoeff(2) = 1.8e8        ! Y0                                           !
            RuptureCoeff(3) = 0.01e9        ! Y1                                           !
            RuptureCoeff(4) = 42.5          ! aexp                                         !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 0.0           ! beta1                                        !
            RuptureCoeff(7) = 36.25         ! alpha2                                       !
            RuptureCoeff(8) = 1.e-6         ! beta2                                        !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            ! Shear based criteria                                                         !
            StressCoeff(1)  = 1.0           ! Yeq_A                                        !
            StressCoeff(2)  = 0.0           ! Yeq_B                                        !
            StressCoeff(3)  = 0.0           ! Yeq_C                                        !
            StressCoeff(4)  = 0.0           ! s0                                           !
            !  EQN%Yeq_mode= 4
            ! -----------------------------------------------------------------------------!
        CASE('PYREXGLASS')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 2230.0                                              ! RHO        !
            ChWSpeed    = (/ 3690.0, 6050.0 /)                      ! CS, CL               !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = 0.5*(/    1.E-2    ,   2./3.  /)           ! MU12,LAM12               !1.E-2
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 1.0           ! K                                            !
            RuptureCoeff(2) = 1.2e9        ! Y0                                            !
            RuptureCoeff(3) = 0.01e9         ! Y1                                          !
            RuptureCoeff(4) = 32.5         ! aexp                                          !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 22.31e-9      ! beta1                                        !
            RuptureCoeff(7) = 34.8         ! alpha2                                        !
            RuptureCoeff(8) = 223.07e-9         ! beta2                                    !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            StressCoeff(1)  =  0.9   !0.38           ! Yeq_A                                        !
            StressCoeff(2)  =  0.05    !0.40           ! Yeq_B                                        !
            StressCoeff(3)  =  0.0    !0.0           ! Yeq_C                                        !
            StressCoeff(4)  =  1.8e+9 !2.8e+9           ! s0                                          !
            ! -----------------------------------------------------------------------------!     
		CASE('PYREXGLASS2D')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 2230.0                                              ! RHO        !
            ChWSpeed    = (/ 3690.0, 6050.0 /)                      ! CS, CL               !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1.E-2    ,   2./3.  /)          ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 1.0           ! K                                            !
            RuptureCoeff(2) = 1.6e9        ! Y0                                            !
            RuptureCoeff(3) = 0.01e9         ! Y1                                          !
            RuptureCoeff(4) = 32.5         ! aexp                                          !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 22.31e-9      ! beta1                                        !
            RuptureCoeff(7) = 34.8         ! alpha2                                        !
            RuptureCoeff(8) = 223.07e-9         ! beta2                                    !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            StressCoeff(1)  = 0.38           ! Yeq_A                                        !
            StressCoeff(2)  = 0.40           ! Yeq_B                                        !
            StressCoeff(3)  = 0.0           ! Yeq_C                                        !
            StressCoeff(4)  = 2.8e+9           ! s0                                           ! 
        CASE('PYREXGLASS1D')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 2230.0                                              ! RHO        !
            ChWSpeed    = (/ 3690.0, 6050.0 /)                      ! CS, CL               !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1.E-2    ,   2./3.  /)          ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 1.0           ! K                                            !
            RuptureCoeff(2) = 1.6e9        ! Y0                                            !
            RuptureCoeff(3) = 0.01e9         ! Y1                                          !
            RuptureCoeff(4) = 32.5         ! aexp                                          !
            ! -------------------------------------                                        !
            RuptureCoeff(5) = 36.25         ! alpha1                                       !
            RuptureCoeff(6) = 22.31e-9      ! beta1                                        !
            RuptureCoeff(7) = 34.8         ! alpha2                                        !
            RuptureCoeff(8) = 223.07e-9         ! beta2                                    !
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            StressCoeff(1)  =  0.45   !0.38           ! Yeq_A                                        !
            StressCoeff(2)  =  0.8    !0.40           ! Yeq_B                                        !
            StressCoeff(3)  =  0.0    !0.0           ! Yeq_C                                        !
            StressCoeff(4)  =  2.0e+9 !2.8e+9           ! s0                                           !
            ! -----------------------------------------------------------------------------!  
        CASE('COPPER1')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff(1) = 8930.0                                              ! RHO        !
            ChWSpeed    = (/ 2325.0, 4760.0 /)                      ! CS, CL               !
            LamCoeff(2) = LamCoeff(1)*ChWSpeed(1)**2                          ! MU         !
            LamCoeff(3) = LamCoeff(1)*(ChWSpeed(2)**2 - 2.0d0*ChWSpeed(1)**2) ! LAM        !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1./1.16666   ,   0.666666   /)                                   ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 0.0           ! K                                            !
            RuptureCoeff(2) = 1E+22         ! Y0                                           !
            RuptureCoeff(3) = 1E+22         ! Y1                                           !
            RuptureCoeff(4) = 1.0           ! aexp                                         !
            ! -----------------------------------------------------------------------------!                     
            RuptureCoeff(5) = 40.            ! alpha1                       
            RuptureCoeff(6) = 0.            ! beta1                           
            RuptureCoeff(7) = 40.            ! alpha2                       
            RuptureCoeff(8) = 0.            ! beta2                        
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            StressCoeff(1)     = 1.0           ! Yeq_A                                        !
            ! -----------------------------------------------------------------------------!              
        CASE('UNBREAKABLE')
            ! -----------------------------------------------------------------------------!
            ! Set Lame parameters and characteristic wave speed ---------------------------!
            LamCoeff    = 0.                                                  ! RHO        !
            ChWSpeed    = 0.                                        ! CS, CL               !
            LamCoeff    = 0.                                                               !
            ! -----------------------------------------------------------------------------!
            ! Set the degradation parameters for the mixture ------------------------------!
            DamCoeff    = (/    1./1.16666   ,   0.666666   /)                                   ! MU12,LAM12               !
            ! -----------------------------------------------------------------------------!
            ! Set the rupture coefficients ------------------------------------------------!
            RuptureCoeff(1) = 0.0           ! K                                            !
            RuptureCoeff(2) = 1E+22         ! Y0                                           !
            RuptureCoeff(3) = 1E+22         ! Y1                                           !
            RuptureCoeff(4) = 1.0           ! aexp                                         !
            ! -----------------------------------------------------------------------------!                     
            RuptureCoeff(5) = 50.            ! alpha1                       
            RuptureCoeff(6) = 0.            ! beta1                           
            RuptureCoeff(7) = 50.            ! alpha2                       
            RuptureCoeff(8) = 0.            ! beta2                        
            ! -----------------------------------------------------------------------------!
            ! Set the equivalent stress coefficients --------------------------------------!
            StressCoeff(1)     = 1.0           ! Yeq_A                                        !
            ! -----------------------------------------------------------------------------!  
        CASE DEFAULT
         !   PRINT *, 'Warning!!! the selected material "',MATERIAL,'" does not exist in module GPRmaterials!'
        END SELECT
    END SUBROUTINE GetMaterialParameters
END MODULE  GPRmaterials
#endif
#endif
