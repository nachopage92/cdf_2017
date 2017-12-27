subroutine PREDICCION_VELOCIDAD(nx,ny,dx,dy,dt,Re,P,u_1,u_0,v_1,v_0,u_pred,v_pred)

	implicit none

	INTERFACE
	
		subroutine WPENS(phi,i,j,phi_W,phi_P,phi_E,phi_N,phi_S)
			integer,intent(in)::i,j
			real(kind=8),intent(in),dimension(:,:) :: phi
			real(kind=8),intent(out)::phi_W,phi_P,phi_E,phi_N,phi_S
		end subroutine
		
		subroutine  CALCULO_RHS(&
				& u_W,u_E,u_P,v_N,v_S,v_P, &
				& u0_W,u0_E,u0_P,v0_N,v0_S,v0_P, &
				& phi_W,phi_E,phi_S,phi_N,phi_P,&
				& phi0_W,phi0_E,phi0_S,phi0_N,phi0_P,&
				& P_1,P_2,dx,dy,dist,dt,Re,RHS)
			implicit none
			real(kind=8),intent(in) :: &
				& u_W,u_E,u_P,v_S,v_N,v_P, &				
				& u0_W,u0_E,u0_P,v0_S,v0_N,v0_P
			real(kind=8),intent(in) :: &
				& phi_W,phi_E,phi_S,phi_N,phi_P, &
				& phi0_W,phi0_E,phi0_S,phi0_N,phi0_P
			real(kind=8),intent(in) :: P_1,P_2		
			real(kind=8),intent(in) :: dx,dy,dist,dt
			real(kind=8),intent(in)::Re
			real(kind=8),intent(out) :: RHS
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
	integer	:: i,j,contador
	real(kind=8) :: P_W,P_E,P_N,P_S,&
		& u_W, u_P, u_E, u_N, u_S , &
		& u0_W,u0_P,u0_E,u0_N,u0_S, &
		& v_W, v_P, v_E, v_N, v_S , &
		& v0_W,v0_P,v0_E,v0_N,v0_S, &
		& RHS
	real(kind=8),dimension((nx-1)*ny):: &
		& RHS_u,diagu_d,diagu_c,diagu_u,dV_u0,dV_u
	
	real(kind=8),dimension((nx-1)*(ny-1)):: &
		& RHS_v,diagv_d,diagv_c,diagv_u,dV_v0,dV_v
		
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!					CALCULO RHS

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!			CALCULO DE RHS EN DIRECCION X PARA U

!	recorre en direccion y (fila)
	contador = 0
	do i=2,ny+1
!		recorre en direccion x (columna)
		do j=3,nx+1
		
			contador = contador + 1
				
!			VARIABLES WE-NS a partir del nodo central
			call WPENS(u_1,i,j,	u_W, u_P, u_E, u_N, u_S )
			call WPENS(u_0,i,j,	u0_W,u0_P,u0_E,u0_N,u0_S)
			call WPENS(v_1,i,j,	v_W, v_P, v_E, v_N, v_S )
			call WPENS(v_0,i,j,	v0_W,v0_P,v0_E,v0_N,v0_S)
		
!			GRADIENTE DE PRESION
			P_W = P(i,j)
			P_E = P(i,j+1)
			P_S = P(i,j)
			P_N = P(i+1,j)
		
!			CALCULO DE RHS EN DIRECCION X PARA U
			call CALCULO_RHS(&
				& u_W,u_E,u_P,v_N,v_S,v_P, &
				& u0_W,u0_E,u0_P,v0_N,v0_S,v0_P, &
				& u_W,u_E,u_S,u_N,u_P,&
				& u0_W,u0_E,u0_S,u0_N,u0_P,&
				& P_W,P_E,dx,dy,dx,dt,Re,RHS)
			RHS_u(contador) = RHS
			
!		fin do direccion x (columna)
		end do
!	fin do direccion y (fila)
	end do

!	-----------------------------------------------

!			CALCULO DE RHS EN DIRECCION X PARA V

!	recorre en direccion y (fila)
	contador = 0
	do i=2,ny
!		recorre en direccion x (columna)
		do j=3,nx+1
		
			contador = contador + 1
				
!			VARIABLES WE-NS a partir del nodo central
			call WPENS(u_1,i,j,	u_W, u_P, u_E, u_N, u_S )
			call WPENS(u_0,i,j,	u0_W,u0_P,u0_E,u0_N,u0_S)
			call WPENS(v_1,i,j,	v_W, v_P, v_E, v_N, v_S )
			call WPENS(v_0,i,j,	v0_W,v0_P,v0_E,v0_N,v0_S)
		
!			GRADIENTE DE PRESION
			P_W = P(i,j)
			P_E = P(i,j+1)
			P_S = P(i,j)
			P_N = P(i+1,j)

!			CALCULO DE RHS EN DIRECCION X PARA V
			call CALCULO_RHS(&
				& u_W,u_E,u_P,v_N,v_S,v_P, &
				& u0_W,u0_E,u0_P,v0_N,v0_S,v0_P, &
				& v_W,v_E,v_S,v_N,v_P,&
				& v0_W,v0_E,v0_S,v0_N,v0_P,&
				& P_S,P_N,dx,dy,dx,dt,Re,RHS)
			RHS_v(contador) = RHS	
			
!		fin do direccion x (columna)
		end do
!	fin do direccion y (fila)
	end do

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!			CALCULO DE dV_u

!		direccion x
	diagu_d = -2._8*dt*dy/(3._8*dx*Re)
	diagu_c = 1._8*dx*dy + 4._8*dt*dy/(3._8*dx*Re)
	diagu_u = -2._8*dt*dy/(3._8*dx*Re)
!		primer paso predictor u-> resuelve dV_u0
	call thomas((nx-1)*ny,diagu_d,diagu_c,diagu_u,dV_u0,RHS_u)
	
	
!		direccion y
	diagu_d = -2._8*dt*dx/(3._8*dy*Re)
	diagu_c = 1._8*dx*dy + 4._8*dt*dx/(3._8*dy*Re)
	diagu_u = -2._8*dt*dx/(3._8*dy*Re)
!		segundo paso predictor -> resuelve dV_u
	call thomas((nx-1)*ny,diagu_d,diagu_c,diagu_u,dV_u,dV_u0)

!	-----------------------------------------------

!			CALCULO DE dV_v

!		direccion x
	diagv_d = -2._8*dt*dy/(3._8*dx*Re)
	diagv_c = 1._8*dx*dy + 4._8*dt*dy/(3._8*dx*Re)
	diagv_u = -2._8*dt*dy/(3._8*dx*Re)
!		segundo paso predictor -> resuelve dV_v0
	call thomas((nx-1)*(ny-1),diagv_d,diagv_c,diagv_u,dV_v,dV_v0)
	
!		direccion y
	diagv_d = -2._8*dt*dx/(3._8*dy*Re)
	diagv_c = 1._8*dx*dy + 4._8*dt*dx/(3._8*dy*Re)
	diagv_u = -2._8*dt*dx/(3._8*dy*Re)
!		segundo paso predictor -> resuelve dV_v
	call thomas((nx-1)*(ny-1),diagv_d,diagv_c,diagv_u,dV_v,dV_v0)


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!			PREDICTOR DEL CAMPO DE VELOCIDADES
!				traspaso variables de vector a matriz

	contador = 0
!	recorre en direccion y (fila)
	do i=2,ny+1
!		recorre en direccion x (columna)
		do j=3,nx+1
			contador = contador + 1
			u_pred(i,j) = u_1(i,j) + dV_u(contador)
		end do
	end do
	
!	-----------------------------------------------
	
!	recorre en direccion y (fila)
	contador = 0
	do i=2,ny
!		recorre en direccion x (columna)
		do j=3,nx+1
			contador = contador + 1
			v_pred(i,j) = v_1(i,j) + dV_v(contador)
		end do
	end do
	
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end subroutine
