subroutine CORRECCION_PRESION(nx,ny,dx,dy,Re,P,phi,u,v)

	implicit none
	
	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy,Re
	real(kind=8),dimension((ny+2),(nx+2)),intent(in) :: phi,u,v
	real(kind=8),dimension((ny+2),(nx+2)),intent(inout) :: P
	
	integer :: i,j
	real(kind=8) :: DIV

!:::::::::::::::::::::::::::::::::::::::::::::::::

	do i=2,ny+1
		do j=2,nx+1
			DIV = (u(i,j+1)-u(i,j))/dx + (v(i,j)-v(i-1,j))/dy
			P(i,j) = P(i,j) + DIV/Re - phi(i,j) 
		end do
	end do
				
end subroutine
