subroutine PHI_DIRECTO(nx,ny,dx,dy,dt,u,v,phi,R)

	implicit none
		
	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy,dt
	real(kind=8),dimension(ny+2,nx+2),intent(in) :: &
		& u,v
	real(kind=8),intent(out) :: R
	real(kind=8),dimension(ny+2,nx+2),intent(inout) :: & 
		& phi
	
	
	integer :: i,j,k,contador
	real(kind=8) :: DIV
	real(kind=8),dimension(ny+2,nx+2) :: phi_


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	phi_ = phi

	do i=2,ny+1
		do j=2,nx+1
			DIV = (u(i,j+1)-u(i,j))/dx + (v(i,j)-v(i-1,j))/dy
			phi_(i,j) = - DIV * (0.75_8/dt) / ( 1._8/dx**2._8 + 1._8/dy**2._8 )
		end do
	end do
	
	phi = phi_

end subroutine
