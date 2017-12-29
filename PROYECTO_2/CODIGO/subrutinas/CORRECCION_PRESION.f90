subroutine CORRECCION_PRESION(nx,ny,dx,dy,P,phi,u,v)

	implicit none
	
	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy
	real(kind=8),dimension((ny+2),(nx+2)),intent(in) :: phi,u,v
	real(kind=8),dimension((ny+2),(nx+2)),intent(inout) :: P
	
	integer :: i,j

!:::::::::::::::::::::::::::::::::::::::::::::::::

	do i=2,ny+1
		do j=2,nx+1
			P(i,j) = P(i,j) + phi(i,j) &
				& - (u(i,j+1)-u(i,j))*dy &
				& - (v(i,j)-v(i-1,j))*dx
		end do
	end do
	
	P(1,:) = P(2,:)
	P(:,1) = P(:,2)
	P(ny+2,:) = P(ny+1,:)
	P(:,nx+2) = P(:,nx+1)
				
end subroutine
