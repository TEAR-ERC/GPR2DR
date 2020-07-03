module ComplexGeometry
    ! Data for the complex 3D profile
    real, allocatable                   :: x_cg(:), y_cg(:),z_cg(:,:)
	real    							:: xL_cg, xR_cg,yL_cg,yR_cg,dx_cg,dy_cg,zMax_cg,zMin_cg
    integer                             :: nx_cg,ny_cg
	
	
	public :: DistanceFromSurfaceCG
    contains
    !*********************************************************************************
    ! Three dimensional complex geometry
    !*********************************************************************************
     real function DistanceFromSurfaceCG(x_in,y_in,z_in,di_size) ! Look at the minimum on the loaded configuration
        implicit none
        real    :: x_in,y_in,z_in
        ! Local variables
        real    :: minvals(3), rr, np(3), np_prj(3), point(3), xx(3,3), nv(3),u(3),v(3),cv(3),tv(3), sign, np_norm,np_nv_norm,di_size
        integer :: minpos(3,2),i,j,k,l,shift(8,2,2)
        integer :: ix_mid,iy_mid,minpos_deb(2)
        point=(/x_in,y_in,z_in /)
        

		
      if(z_in<zMin_cg-10.0*di_size) then
           DistanceFromSurfaceCG=-1.e+14 
           return
       end if
       if(z_in>zMax_cg+10.0*di_size) then
           DistanceFromSurfaceCG=1.e+14 
           return
       end if
		!print *, "here"
        ! Fast way
        ix_mid=floor((x_in-xL_cg)/dx_cg)
        iy_mid=floor((y_in-yL_cg)/dy_cg)
        minvals=1.e+14
        minpos=0
        do i=max(ix_mid-40,1),min(nx_cg,ix_mid+40)
            do j=max(iy_mid-40,1),min(ny_cg,iy_mid+40)
               rr=sqrt((point(1)-x_cg(i))**2+(point(2)-y_cg(j))**2+(point(3)-z_cg(i,j))**2)
               if(rr<minvals(1)) then
                   minvals(1)=rr
                   minpos(1,:)=(/i,j/);
                end if
            end do
        end do 
		
		
				! Debug
		!print *,ix_mid,iy_mid,nx_cg,ny_cg
		!i=ix_mid
		!do j=max(iy_mid-40,1),min(ny_cg,iy_mid+40)
		!   print *,z_cg(i,j)
		!end do
        !
		!pause
		
		if(point(3)-z_cg(minpos(1,1),minpos(1,2))>0) then
			DistanceFromSurfaceCG=minvals(1)
		else
			DistanceFromSurfaceCG=-minvals(1)
		end if
		if(minvals(1)>2.0*di_size) then
			!rr=RadiusFromCG(point(1),point(2),point(3))
			return
		end if
        !if(maxval(abs(minpos_deb-minpos(1,:)))>0) then
        !    print *, 'Warning Distance search!'
        !    print *, minpos(1,:),minpos_deb
        !    print *, 'ix=',ix_mid,'i_est=',(x_in-xL_cg)/dx_cg
        !    print *, x_in,xL_cg,dx_cg
        !    print *, 'ix=',iy_mid,'i_est=',(y_in-yL_cg)/dy_cg
        !    print *, y_in,yL_cg,dy_cg
        !    pause
        !end if
		
		!rr=RadiusFromCG(point(1),point(2),point(3))
		!if(rr>0) then
		!	sign=+1.0    
		!else
		!	sign=-1.0     
		!end if
		
        !DistanceFromSurfaceCG=1.e+14
        !shift(1,1,:)=(/ 0,1  /)
        !shift(1,2,:)=(/ -1,0  /)
        !shift(2,1,:)=(/ 1,0  /)
        !shift(2,2,:)=(/ 0,1  /) 
        !shift(3,1,:)=(/ 0,-1  /)
        !shift(3,2,:)=(/ 1,0  /)
        !shift(4,1,:)=(/ -1,0  /)
        !shift(4,2,:)=(/ 0,-1  /)
        shift(1,1,:)=(/ 0,1  /)
        shift(1,2,:)=(/ -1,1  /)
        shift(2,1,:)=(/ 1,0  /)
        shift(2,2,:)=(/ 0,1  /) 
        shift(3,1,:)=(/ 1,-1  /)
        shift(3,2,:)=(/ 1,0  /)
        shift(4,1,:)=(/ 0,-1  /)
        shift(4,2,:)=(/ 1,-1  /)
        shift(5,1,:)=(/ -1,0  /)
        shift(5,2,:)=(/ 0,-1  /)
        shift(6,1,:)=(/ -1,1  /)
        shift(6,2,:)=(/ -1,0  /)
		
		np_nv_norm=100.0
        do l=1,6
            minpos(2,:)=minpos(1,:)+shift(l,1,:)
            minpos(3,:)=minpos(1,:)+shift(l,2,:)
            if(minval(minpos) .eq. 0 .or. maxval(minpos(:,1))>nx_cg  .or. maxval(minpos(:,2))>ny_cg  ) then
                cycle
            end if    

            ! COnstruct the triangle composed by the three minimum
            do k=1,3
                xx(k,:)=[x_cg(minpos(k,1)),y_cg(minpos(k,2)),z_cg(minpos(k,1),minpos(k,2))];
            end do
            ! COmpute the normal and the tangential vectors
            u=xx(3,:)-xx(1,:)
            v=xx(2,:)-xx(1,:)
            nv(1)=u(2)*v(3)-u(3)*v(2)
            nv(2)=u(3)*v(1)-u(1)*v(3)
            nv(3)=u(1)*v(2)-u(2)*v(1)
            nv=-nv/sqrt(sum(nv**2))  ! Normalize nv
        
            np_norm=dot_product((point-xx(1,:)),nv)
            np_prj=point-np_norm*nv
            
            if(np_norm>0) then
                sign=+1.0    
            else
                sign=-1.0     
            end if
            np=np_prj
            call GetTirangleMinPoint(xx(:,1),xx(:,2),xx(:,3),np)
            if(abs(np_norm)<5.e-1*di_size .and. maxval(np-np_prj)>1.e-6) then ! Normal vector perturbation problem
                cycle    
            end if
            
            rr=sqrt(sum(  (np-point)**2  ))
            if(rr.le. abs(DistanceFromSurfaceCG)) then
                if(abs(rr-abs(DistanceFromSurfaceCG))<1.e-6 .and. abs(np_norm)<abs(np_nv_norm)) then
                    cycle    
                else
                    DistanceFromSurfaceCG=sign*rr
                    np_nv_norm=np_norm
                end if
            end if
        end do
    end function DistanceFromSurfaceCG
    
    
    subroutine GetTirangleMinPoint(XX,YY,ZZ,PP)
        implicit none
        real    :: XX(3),YY(3),ZZ(3),PP(3)
        real    :: xi, gamma, xi_e, gamma_e, x_out(3), alpha
        
    gamma=((PP(1)-XX(1))*(YY(2)-YY(1))-(XX(2)-XX(1))*(PP(2)-YY(1)))/((XX(3)-XX(1))*(YY(2)-YY(1))+(YY(1)-YY(3))*(XX(2)-XX(1)))
    xi=((PP(1)-XX(1))*(YY(3)-YY(1))-(XX(3)-XX(1))*(PP(2)-YY(1)))/((XX(2)-XX(1))*(YY(3)-YY(1))+(YY(1)-YY(2))*(XX(3)-XX(1)))
    
    if(xi>0 .and. gamma > 0 .and. xi+gamma<1) then
        xi_e=xi
        gamma_e=gamma
    elseif(xi<=0) then
      if(gamma<=0) then
          xi_e=0
          gamma_e=0
      elseif(gamma <=1) then
          xi_e=0
          gamma_e=gamma      
      else
          xi_e=0
          gamma_e=1          
      end if
    elseif(gamma<=0) then
      if(xi <=1) then
          xi_e=xi
          gamma_e=0         
      else
          xi_e=1
          gamma_e=0;        
      end if       
    else
        alpha=0.5*(gamma-xi+1)
        if(alpha<=0) then
          xi_e=1
          gamma_e=0              
        elseif(alpha>=1) then
          xi_e=0
          gamma_e=1             
        else
          xi_e=0.5*(1-alpha)
          gamma_e=0.5*(1+alpha)              
        end if
    end if
    
    x_out(1)=XX(1)+(XX(2)-XX(1))*xi_e+(XX(3)-XX(1))*gamma_e
    x_out(2)=YY(1)+(YY(2)-YY(1))*xi_e+(YY(3)-YY(1))*gamma_e
    x_out(3)=ZZ(1)+(ZZ(2)-ZZ(1))*xi_e+(ZZ(3)-ZZ(1))*gamma_e
    
    PP=x_out
    end subroutine GetTirangleMinPoint
    
	
	
    real function RadiusFromCG(x_in,y_in,z_in)
		USE, INTRINSIC :: ISO_C_BINDING
        implicit none
        real    :: x_in,y_in,z_in,z_out
        real    :: i,j,phi(4),xi,gamma
        integer :: ix(2)
        
        ix=lookatindex_cg_fast(x_in,y_in)
        phi(1)=z_cg(ix(1),ix(2))
        phi(2)=z_cg(ix(1)+1,ix(2))
        phi(3)=z_cg(ix(1),ix(2)+1)
        phi(4)=z_cg(ix(1)+1,ix(2)+1)
        xi=(x_in-x_cg(ix(1)))/(x_cg(ix(1)+1)-x_cg(ix(1)))
        gamma=(y_in-y_cg(ix(2)))/(y_cg(ix(2)+1)-y_cg(ix(2)))
        z_out=(1-xi)*(1-gamma)*phi(1)+xi*(1-gamma)*phi(2)+gamma*(1-xi)*phi(3)+xi*gamma*phi(4)
        !z_out=z_cg(ix(1),ix(2))+z_cg(ix(1)+1,ix(2))+z_cg(ix(1),ix(2)+1)+z_cg(ix(1)+1,ix(2)+1)
        ! Reverse for negative representation (z down)
        !z_out=-z_out
        !RadiusFromCG=z_out-z_in
        RadiusFromCG=-z_out+z_in
    end function RadiusFromCG

    
    function lookatindex_cg(x_in,y_in)
        implicit none
        real    :: x_in,y_in
        integer    :: lookatindex_cg(2)
        integer :: i,j
        
        lookatindex_cg=-1
        do i=1,nx_cg-1
            if(x_cg(i).le.x_in .and. x_in .le. x_cg(i+1)) then
               lookatindex_cg(1)=i
               exit
            end if
        end do
        do j=1,ny_cg-1
            if(y_cg(j).le.y_in .and. y_in .le. y_cg(j+1)) then
               lookatindex_cg(2)=j
               exit
            end if
        end do   
        if(minval(lookatindex_cg(:))<0) then
            print *, 'LookIndex error for x_in=',x_in, ' and y_in=',y_in, '! Please choose a larger CG domain'
            stop
        end if
    end function lookatindex_cg
    
    function lookatindex_cg_fast(x_in,y_in)
        implicit none
        real    :: x_in,y_in
        integer    :: lookatindex_cg_fast(2)
        integer :: i,j
        
		
		
		lookatindex_cg_fast(1)=floor((x_in-xL_cg)/dx_cg)+1
        lookatindex_cg_fast(2)=floor((y_in-yL_cg)/dy_cg)+1
		!print *, 'DyL_cg=',y_in-yL_cg,x_in-xL_cg,lookatindex_cg_fast
        if(minval(lookatindex_cg_fast(1:2))<0) then
            print *, 'LookIndex error for x_in=',x_in, ' and y_in=',y_in, '! Please choose a larger CG domain'
			print *, 'x_L=',xL_cg,'y_L=',yL_cg,'dx=',dx_cg,'dy=',dy_cg
            stop
        end if
    end function lookatindex_cg_fast
end module ComplexGeometry
