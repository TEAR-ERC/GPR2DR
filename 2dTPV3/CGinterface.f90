	! ================ Specific complex geometries routines =======================
!RECURSIVE subroutine LoadCGFile(MyOffset,MyDomain,STRLEN,PARSETUP,CX,CY,ISBINARY)
RECURSIVE subroutine LoadCGFile(MyOffset,MyDomain,CX,CY,ISBINARY)
		USE :: MainVariables, only : nDim,USECG
		USE :: ComplexGeometry
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
		!INTEGER				 , INTENT(IN) :: STRLEN
		!CHARACTER(LEN=STRLEN), INTENT(IN) :: PARSETUP
		REAL				 , INTENT(IN) :: MyOffset(3),MyDomain(3),CX,CY
		INTEGER							  :: ISBINARY
		
		CHARACTER(LEN=100) :: CGEOMFile
		
        integer     :: jcg,nx_cg_new,ny_cg_new,i,j, ix(2)
        real        :: leng(2,2),center(2)
        integer     :: n_new_in(2)
        real, allocatable   :: x_cg_new(:),y_cg_new(:),z_cg_new(:,:)
        real            :: h, phi(4), xi,gamma
        real            :: minx,maxx,miny,maxy
		logical			:: invert_coordinates, binary_input
		real			::  scalefactor(3)
		! Input parameters
		leng(1:2,1)=MyDomain(1:2)+MyOffset(1:2)
		leng(1:2,2)=-MyOffset(1:2)
		print *, "**********************************************************"
		print *, MyDomain
		print *, MyOffset
		print *, "**********************************************************"
	if(nDim .eq. 3) then
		USECG=.true.
		!print *,Cx,Cy, ISBINARY
		!pause
		!leng=15000.0
		!n_new_in=(/200, 200/)			! Number of elements for the min sub tri function
		!
		!CGEOMFile="CG.dat"			! DTM file
		!center=(/0.0, 0.0/)			! UTM coordinates of the center (with respect to the DTM data file)
		!binary_input=.false.
		
		CGEOMFile="DTM/trient_003_44_48_9_13.bin"			! DTM file
		center=(/4405.905971174,2551.552691730/)			! UTM coordinates of the center (with respect to the DTM data file)
		binary_input=.true.
		
		!CGEOMFile= TRIM (PARSETUP)
		CENTER(1)=CX
		CENTER(2)=CY
		IF(ISBINARY.EQ.1) THEN
			BINARY_INPUT=.TRUE.
		ELSE
			BINARY_INPUT=.FALSE.
		END IF
		n_new_in=(/200, 200/)
		
		
		!leng=15000.0
		!center=(/600.0, 5110.0/)			! UTM coordinates of the center (with respect to the DTM data file)
		!n_new_in=(/200, 200/)			! Number of elements for the min sub tri function
		!CGEOMFile="alps_01.txt"			! DTM file
		
		if(binary_input) then
			open(8, file=trim(CGEOMFile) ,form='unformatted')
				read(8) nx_cg
				read(8) ny_cg
				allocate(x_cg(nx_cg),y_cg(ny_cg),z_cg(nx_cg,ny_cg))
				read(8) scalefactor(1:3)
				read(8) x_cg(1:nx_cg)
				read(8) y_cg(1:ny_cg)
				read(8) z_cg      
			close(8)			
		else
			open(8, file=trim(CGEOMFile), action='read')
				read(8,*) nx_cg
				read(8,*) ny_cg
				allocate(x_cg(nx_cg),y_cg(ny_cg),z_cg(nx_cg,ny_cg))
				read(8,*) scalefactor(1:3)
				read(8,*) x_cg(1:nx_cg)
				read(8,*) y_cg(1:ny_cg)
				do jcg=1,ny_cg
					read(8,*) z_cg(1:nx_cg,jcg)       
				end do
			close(8)
		end if
		print *, 'Min-Max of z (DTM)=',minval(z_cg), maxval(z_cg)
		!print *, 'nx,ny=',nx_cg, ny_cg,scalefactor,z_cg(5000,5000)
		center(1)=center(1)*scalefactor(1);
		center(2)=center(2)*scalefactor(2);
		x_cg=x_cg*scalefactor(1)
		y_cg=y_cg*scalefactor(2)
		z_cg=z_cg*scalefactor(3)
        if(y_cg(2)<y_cg(1)) then
			invert_coordinates = .true.
			y_cg(1:nx_cg)=y_cg(nx_cg:1:-1)
			z_cg(nx_cg,1:ny_cg)=z_cg(nx_cg,ny_cg:1:-1)
		else
			invert_coordinates = .false.
		end if        
        ! Connect the input parameters 
        maxx=center(1)+leng(1,1)
        minx=center(1)-leng(1,2)
        maxy=center(2)+leng(2,1)
        miny=center(2)-leng(2,2)
        
        nx_cg_new=n_new_in(1);
        ny_cg_new=n_new_in(2);
        allocate(x_cg_new(nx_cg_new),y_cg_new(ny_cg_new),z_cg_new(nx_cg_new,ny_cg_new))
        h=(maxx-minx)/(nx_cg_new-1)
        x_cg_new(1)=minx
        do i=2,nx_cg_new
            x_cg_new(i)=x_cg_new(i-1)+h   
        end do
        h=(maxy-miny)/(ny_cg_new-1)
        y_cg_new(1)=miny
        do i=2,ny_cg_new
            y_cg_new(i)=y_cg_new(i-1)+h   
        end do

        do i=1,nx_cg_new
            do j=1,ny_cg_new
                ix=lookatindex_cg(x_cg_new(i),y_cg_new(j))
                phi(1)=z_cg(ix(1),ix(2))
                phi(2)=z_cg(ix(1)+1,ix(2))
                phi(3)=z_cg(ix(1),ix(2)+1)
                phi(4)=z_cg(ix(1)+1,ix(2)+1)
                xi=(x_cg_new(i)-x_cg(ix(1)))/(x_cg(ix(1)+1)-x_cg(ix(1)))
                gamma=(y_cg_new(j)-y_cg(ix(2)))/(y_cg(ix(2)+1)-y_cg(ix(2)))
                z_cg_new(i,j)=(1-xi)*(1-gamma)*phi(1)+xi*(1-gamma)*phi(2)+gamma*(1-xi)*phi(3)+xi*gamma*phi(4)
            end do
        end do
        x_cg_new=x_cg_new-center(1) ! The center for the DTM is (0,0) in the fortran code
        y_cg_new=y_cg_new-center(2) ! The center for the DTM is (0,0) in the fortran code
        
        deallocate(x_cg,y_cg,z_cg)
        nx_cg=nx_cg_new;
        ny_cg=ny_cg_new;        
        allocate(x_cg(nx_cg),y_cg(ny_cg),z_cg(nx_cg,ny_cg))
        x_cg=x_cg_new
        y_cg=y_cg_new!-243.9     ! Move by 1 element since CG is shifted with respect to the real DTM
		if(invert_coordinates) then
			do i=1,nx_cg
				z_cg(i,1:ny_cg)=z_cg_new(i,ny_cg:1:-1)
			end do
		else
			do i=1,nx_cg
				z_cg(i,1:ny_cg)=z_cg_new(i,1:ny_cg)
			end do
		end if
		
        xL_cg=x_cg(1)
        xR_cg=x_cg(nx_cg)
        yL_cg=y_cg(1)
        yR_cg=y_cg(ny_cg)
        dx_cg=(xR_cg-xL_cg)/nx_cg
        dy_cg=(yR_cg-yL_cg)/ny_cg
        zMin_cg=minval(z_cg)
        zMax_cg=maxval(z_cg)
		print *, nx_cg, ny_cg
		deallocate(x_cg_new,y_cg_new,z_cg_new)
		print *, 'Min-Max of x=',minval(x_cg), maxval(x_cg)
		print *, 'Min-Max of y=',minval(y_cg), maxval(y_cg)
		print *, 'Min-Max of z=',minval(z_cg), maxval(z_cg)
		!print *, maxy,miny
		!stop
	end if
end subroutine LoadCGFile

