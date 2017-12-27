subroutine CORRECCION_VELOCIDAD(nx,ny,dx,dy,dt,u_pred,v_pred,phi,u_0,u_1,v_0,v_1)

	implicit none
	
	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy,dt
	real(kind=8),dimension(ny+2,nx+2),intent(in) :: phi
	real(kind=8),dimension(ny+2,nx+2),intent(in) :: u_pred,v_pred
	real(kind=8),dimension(ny+2,nx+2),intent(inout) :: u_0,u_1,v_0,v_1
	
	integer :: i,j

	u_0 = u_1
	v_0 = v_1
		
	do i=2,ny+1
		do j=3,nx+1
			u_1(i,j) = u_pred(i,j) - dt*(phi(i,j)-phi(i,j-1))/dx
		end do
	end do
	
	
	do i=2,ny
		do j=3,nx+1
			v_1(i,j) = v_pred(i,j) - dt*(phi(i+1,j)-phi(i,j))/dy
		end do
	end do
	
end subroutine
