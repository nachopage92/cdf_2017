!	CALCULO DEL RHS DEL PRIMER PASO
!		primer paso para u		
!		primer paso para v

subroutine RIGHT_HAND_SIDE(nx,ny,dx,dy,dt,Re,P,u_1,u_0,v_1,v_0,RHS_y,RHS_x)

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
	
	END INTERFACE

	integer,intent(in)::nx,ny
	real(kind=8),intent(in)::dx,dy,dt,Re
	real(kind=8),dimension(nx+2,ny+2),intent(in):: &
		& P , u_1 , u_0 , v_1 , v_0
	real(kind=8),dimension(nx*ny),intent(out):: RHS_x,RHS_y
	
	!variable locales
	integer	:: i,j,contador
	real(kind=8) :: P_W,P_E,P_N,P_S,&
		& u_W, u_P, u_E, u_N, u_S , &
		& u0_W,u0_P,u0_E,u0_N,u0_S, &
		& v_W, v_P, v_E, v_N, v_S , &
		& v0_W,v0_P,v0_E,v0_N,v0_S, &
		& RHS

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	contador = 0
	
!	recorre en direccion y (fila)
	do i=2,nx+1
!		recorre en direccion x (columna)
		do j=2,ny+1
		
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
			RHS_x(contador) = RHS
			
!			CALCULO DE RHS EN DIRECCION X PARA V
			call CALCULO_RHS(&
				& u_W,u_E,u_P,v_N,v_S,v_P, &
				& u0_W,u0_E,u0_P,v0_N,v0_S,v0_P, &
				& v_W,v_E,v_S,v_N,v_P,&
				& v0_W,v0_E,v0_S,v0_N,v0_P,&
				& P_S,P_N,dx,dy,dx,dt,Re,RHS)
			RHS_y(contador) = RHS			
			
!		fin do direccion x (columna)
		end do
!	fin do direccion y (fila)
	end do

!FIN PRIMER PASO

end subroutine
