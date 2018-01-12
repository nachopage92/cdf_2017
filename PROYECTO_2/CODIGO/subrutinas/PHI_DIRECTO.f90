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
	real(kind=8) :: DIV,R_vec(ny*(nx+1))


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	do i=2,ny+1
		do j=2,nx+1
			DIV = (u(i,j+1)-u(i,j))/dx + (v(i,j)-v(i-1,j))/dy
			phi(i,j) = - DIV * (3._8/(4._8*dt)) / ( 1._8/dx**2._8 + 1._8/dy**2._8 )
		end do
	end do
	
	
	
		!CALCULO DE RESTO R (CONVERGENCIA)
	
	contador = 0
	
	do i=2,ny+1
		do j=1,nx+1
	
			contador = contador + 1

			DIV = (u(i,j+1)-u(i,j))/dx + (v(i,j)-v(i-1,j))/dy

			R_vec(contador) = &
			& (phi(i,j-1)-2._8*phi(i,j)+phi(i,j+1))/dx**2._8 &
			& + (phi(i-1,j)-2._8*phi(i,j)+phi(i+1,j))/dy**2._8 &
			& - 3._8/(2._8*dt)*DIV

		end do
	end do
			
	R = maxval(abs(R_vec))

end subroutine
