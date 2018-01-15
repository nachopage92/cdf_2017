!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine PREDICCION_VELOCIDAD(P,u_1,u_0,v_1,v_0,u_pred,v_pred)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	use variables

	INTERFACE
	
		subroutine  CALCULO_RHS(u,v,u0,v0,P,i,j,str,RHS)
			use variables
			implicit none
			!entrada
			integer,intent(in) :: i,j
			real(kind=8),dimension(ny+2,nx+2),intent(in) ::&
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
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	real(kind=8),dimension(ny+2,nx+2),intent(in):: &
		& P , u_1 , u_0 , v_1 , v_0
	real(kind=8),dimension(ny+2,nx+2),intent(out):: &
		& u_pred , v_pred
		
		!variable locales
	integer	:: i,j,k,contador
	real(kind=8) :: RHS,alfa,beta
	
	real(kind=8),dimension(nx-1):: &
		& diagux_d,diagux_c,diagux_u,dVx_u
	real(kind=8),dimension(ny):: &
		& diaguy_d,diaguy_c,diaguy_u,dVy_u,RHSy_u
		
	real(kind=8),dimension(nx-2):: &
		& diagvx_d,diagvx_c,diagvx_u,dVx_v
	real(kind=8),dimension(ny-1):: &
		& diagvy_d,diagvy_c,diagvy_u,dVy_v,RHSy_v
		
	real(kind=8),dimension(ny+2,nx+2) :: dV
	
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!					CALCULO RHS_U

!	-----------------------------------------------

!			CALCULO EN LA DIRECCION Y

!	-----------------------------------------------

	dV = 0._8

!		recorre cada columna (x1,x2,...)
	do j=3,nx+1
		
		contador = 0
		
!			recorre cada fila (y1,y2,...)
		do i=2,ny+1
			contador = contador + 1
			call CALCULO_RHS(u_1,v_1,u_0,v_0,P,i,j,'u',RHS)
			RHSy_u(contador) = RHS
		end do
		
!			se consideran los esfuerzos cortante
!			debido a la condicion de no deslizamiento
		RHSy_u(1) = RHSy_u(1) - 4._8*pi*dx*(Ly-dx*0.5_8)/(dy*Re) * u_1(2,j) 

!				CALCULO DE dV_u EN CADA COLUMNA	
!			diagonales de la matriz tridiagonal	
		alfa = -2._8*dt*dx/(3._8*dy*Re)
		beta =  dx*dy + 4._8*dt*dx/(3._8*dy*Re)
		diaguy_u = alfa
		diaguy_c = beta
		diaguy_d = alfa
		
!			primer paso predictor u
		call thomas(ny,diaguy_d,diaguy_c,diaguy_u,dVy_u,RHSy_u)
		dV(2:ny+1,j) = dVy_u(:)
			
	end do


!	-----------------------------------------------

!			CALCULO EN LA DIRECCION X
 
!	-----------------------------------------------
 
!		se recorre cada fila
	do i=2,ny+1

!				CALCULO DE dV_u EN CADA FILA	
!			diagonales de la matriz tridiagonal	
		alfa = -2._8*dt/(3._8*dx**2._8*Re)
		beta =  1._8 + 4._8*dt/(3._8*dx**2._8*Re)
		diagux_u = alfa
		diagux_c = beta
		diagux_d = alfa
			
!			segundo paso predictor u
		call thomas(nx,diagux_d,diagux_c,diagux_u,dVx_u,dV(i,3:nx+1))
		dV(i,3:nx+1) = dVx_u(:)
		
	end do

!		calculo de la componente horizontal de la velocidad
	u_pred = u_1 + dV
	
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!					CALCULO RHS_V

!	-----------------------------------------------

!			CALCULO EN LA DIRECCION Y

	dV = 0._8

!		se recorre cada columna
	do j=3,nx
		
		contador = 0
		
!			se calcula cada fila
		do i=2,ny
			contador = contador + 1
			call CALCULO_RHS(u_1,v_1,u_0,v_0,P,i,j,'v',RHS)
			RHSy_v(contador) = RHS
		end do

!				CALCULO DE dV_v EN CADA COLUMNA	
!			diagonales de la matriz tridiagonal	
		alfa = -2._8*dt*dx/(3._8*dy*Re)
		beta =  dx*dy + 4._8*dt*dx/(3._8*dy*Re)
		diagvy_u = alfa
		diagvy_c = beta
		diagvy_d = alfa
		
!			primer paso predictor v
		call thomas(ny-1,diagvy_d,diagvy_c,diagvy_u,dVy_v,RHSy_v)
		dV(2:ny,j) = dVy_v(:)
			
	end do


!	-----------------------------------------------

!			CALCULO EN LA DIRECCION X
 
!		se recorre cada fila
	do i=2,ny

!				CALCULO DE dV_u EN CADA FILA	
!			diagonales de la matriz tridiagonal	
		alfa = -2._8*dt/(3._8*dx**2._8*Re)
		beta =  1._8 + 4._8*dt/(3._8*dx**2._8*Re)
		diagvx_u = alfa
		diagvx_c = beta
		diagvx_d = alfa
			
!			segundo paso predictor v
		call thomas(nx-2,diagvx_d,diagvx_c,diagvx_u,dVx_v,dV(i,3:nx))
		dV(i,3:nx) = dVx_v(:)
		
	end do
	
!		calculo de la componente vertical de la velocidad
	v_pred = v_1 + dV


!!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end subroutine
