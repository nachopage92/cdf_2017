subroutine CORRECCION_PRESION(nx,ny,dx,dy,P,phi,u,v)

	implicit none
	
	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy
	real(kind=8),dimension((ny+2),(nx+2)),intent(in) :: phi,u,v
	real(kind=8),dimension((ny+2),(nx+2)),intent(inout) :: P
	
	integer :: i,j

	do i=2,ny+2
		do j=1,nx+1
			P(i,j) = P(i,j) + phi(i,j) &
				& + (u(i,j)-u(i,j+1))*dy &
				& + (v(i-1,j)-v(i,j))*dx
		end do
	end do
	
	P(1,:) = P(2,:)
	P(ny+2,:) = P(ny+1,:)
	P(:,nx+2) = P(:,nx+1)
	
end subroutine
