subroutine CORRECCION_PRESION(nx,ny,dx,dy,P,phi,u,v)

	implicit none
	
	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy
	real(kind=8),dimension((nx+2),(ny+2)),intent(in) :: phi,u,v
	real(kind=8),dimension((nx+2),(ny+2)),intent(inout) :: P
	
	integer :: i,j

	do i=1,nx+1
		do j=1,ny+1
			P(i,j) = (P(i,j)+phi(i,j))/(dx*dy) &
				& + (u(i,j)-u(i,j+1))*dx &
				& + (v(i-1,j)-v(i,j))*dy
		end do
	end do
	
end subroutine
