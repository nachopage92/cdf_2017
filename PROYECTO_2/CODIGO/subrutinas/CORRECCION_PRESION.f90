subroutine CORRECCION_PRESION(nx,ny,dx,dy,dt,Re,P,u,v,R,phi)

	implicit none
	
	!entrada
	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy,dt,Re
	real(kind=8),dimension((ny+2),(nx+2)),intent(in) :: u,v
	
	!entrada/salida
	real(kind=8),dimension((ny+2),(nx+2)),intent(inout) :: P
	
	!salida
	real(kind=8),intent(out) :: R
	real(kind=8),dimension((ny+2),(nx+2)),intent(out) :: phi
	
	!local
	integer :: i,j,k,contador
	real(kind=8),dimension(ny*(nx+1)) :: R_vec

	phi = 0._8

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	!resolucion ecuacion de Poisson para PHI

	do i=2,ny+1
		do j=1,nx+1
			phi(i,j) = - (3._8/(4._8*dt)) * DIV(i,j,u,v,dx,dy) / ( 1._8/dx**2._8 + 1._8/dy**2._8 )
		end do
	end do
	

!:::::::::::::::::::::::::::::::::::::::::::::::::

	!correccion de la presion
	
	do i=2,ny+1
		do j=1,nx+1
			P(i,j) = P(i,j) - DIV(i,j,u,v,dx,dy)/Re + phi(i,j)
		end do
	end do


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::	

	!calculo de convergencia
	
	contador = 0
	do i=2,ny+1
		do j=1,nx+1
			contador = contador + 1
			R_vec(contador) =  LAP(i,j,phi,dx,dy) 
		end do
	end do
			
	R = maxval(abs(R_vec))
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	contains
	
	real(kind=8) function DIV(i,j,u,v,dx,dy)
		integer :: i,j
		real(kind=8) :: dx,dy
		real(kind=8),dimension(:,:) :: u,v
		DIV = (u(i,j+1)-u(i,j))/dx + (v(i,j)-v(i-1,j))/dy
	end function
	
	real(kind=8) function LAP(i,j,phi,dx,dy)
		integer :: i,j
		real(kind=8) :: dx,dy
		real(kind=8),dimension(:,:) :: phi
		LAP= (phi(i,j-1)-2._8*phi(i,j)+phi(i,j+1))/dx**2._8 &
			& + (phi(i-1,j)-2._8*phi(i,j)+phi(i+1,j))/dy**2._8
	end function
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end subroutine
