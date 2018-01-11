subroutine PREDICCION_VELOCIDAD(nx,ny,dx,dy,dt,Re,P,u_1,u_0,v_1,v_0,u_pred,v_pred)

	implicit none

	INTERFACE
			
		subroutine  CALCULO_RHS(u,v,u0,v0,P,i,j,dx,dy,dt,str,RHS)
			implicit none
			!entrada
			integer,intent(in) :: i,j
			real(kind=8)::dx,dy,dt
			real(kind=8),dimension(:,:),intent(in) ::&
				& u,v,u0,v0,P
			character(len=1),intent(in) :: str
			!salida
			real(kind=8),intent(out)::RHS
		end subroutine
				
		subroutine thomas(n,d,a,c,x,b)
			implicit none
			integer,intent(in)::n
			double precision,dimension(n),intent(in) :: d,a,c,b
			double precision,dimension(n),intent(out)	:: x
		end subroutine

	
	END INTERFACE

	integer,intent(in)::nx,ny
	real(kind=8),intent(in)::dx,dy,dt,Re
	real(kind=8),dimension(ny+2,nx+2),intent(in):: &
		& P , u_1 , u_0 , v_1 , v_0
	real(kind=8),dimension(ny+2,nx+2),intent(out):: &
		& u_pred , v_pred
		
		!variable locales
	integer	:: i,j,k,contador
	real(kind=8) :: &
		& RHS,alfa,beta
		
	real(kind=8),dimension(nx-2):: &
		& RHSx_u,diagux_d,diagux_c,diagux_u,dVx_u
		
	real(kind=8),dimension(ny):: &
		& diaguy_d,diaguy_c,diaguy_u,dVy_u
	
!	real(kind=8),dimension(nx-1):: &
!		& RHS_v,diagv_d,diagv_c,diagv_u,dV_v0,dV_v
!	real(kind=8),dimension(ny-1):: &
!		& RHS_v,diagv_d,diagv_c,diagv_u,dV_v0,dV_v
!	real(kind=8),dimension(ny):: &
!		& RHS_v,diagv_d,diagv_c,diagv_u,dV_v0,dV_v
		
	real(kind=8),dimension(ny+2,nx+2) :: dV
	
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!					CALCULO RHS_U

!	-----------------------------------------------

!			CALCULO EN LA DIRECCION X

 dV = 0._8

!	en cada fila...
	do i=2,ny+1
	
		contador = 0
		
!		...recorre las columnas (direccion x)
		do j=3,nx+1
			contador = contador + 1
			call CALCULO_RHS(u_1,v_1,u_0,v_0,P,i,j,dx,dy,dt,'u',RHS)
			RHSx_u(contador) = RHS
		end do
		
		!considera las 
		RHSx_u(2) = RHSx_u(2) - 4._8*dt*dx/(3._8*Re*dy)*u_1(i,2)
		RHSx_u(nx) = RHSx_u(nx) - 4._8*dt*dx/(3._8*Re*dy)*u_1(i,nx)

!			CALCULO DE dV_u EN CADA FILA	
!		diagonales de la matriz tridiagonal	
		alfa = -2._8*dt*dy/(3._8*dx*Re)
		beta =  1._8*dx*dy + 4._8*dt*dy/(3._8*dx*Re)
		diagux_u = alfa
		diagux_c = beta
		diagux_d = alfa
!		diagux_c(nx-1)=beta+alfa
	
!			primer paso predictor u-> resuelve dV_u0
		call thomas(nx-2,diagux_d,diagux_c,diagux_u,dVx_u,RHSx_u)
		dV(i,3:nx) = dVx_u(:)
		write(3,*) dVx_u(:)
		
!	fin do direccion y (fila)
	end do

	do i=1,ny+2
		write(1,*) dV(i,:)
	end do

	
!!			CALCULO EN LA DIRECCION Y

!!	en cada columna...
!	do j=3,nx+1

!!			CALCULO DE dV_u EN CADA FILA	
!!		diagonales de la matriz tridiagonal	
!		alfa = -2._8*dt*dx/(3._8*dy*Re)
!		beta =  1._8*dx*dy + 4._8*dt*dx/(3._8*dy*Re)
!		diaguy_u = alfa
!		diaguy_c = beta
!		diaguy_d = alfa
!!		diaguy_c(1)=beta-alfa
!!		diaguy_c(ny)=beta+alfa
!	
!!			primer paso predictor u-> resuelve dV_u0
!		call thomas(ny,diaguy_d,diaguy_c,diaguy_u,dVy_u,dV(2:ny,j))
!		u_pred(2:ny+1,j) = dVy_u(:)
!				
!!	fin do direccion y (fila)
!	end do

!!	PENDIENTE!!!
!	v_pred=0._8

!return

!!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end subroutine
