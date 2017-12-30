subroutine PREDICCION_VELOCIDAD(nx,ny,dx,dy,dt,Re,P,u_1,u_0,v_1,v_0,u_pred,v_pred)

	implicit none

	INTERFACE
	
		subroutine WPENS_U(u,v,P,i,j,F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S)
			implicit none
			integer,intent(in)::i,j
			real(kind=8),intent(in),dimension(:,:) :: u,v,P
			real(kind=8),intent(out):: &
				&F_w,F_e,F_n,F_s,u_W,u_E,u_N,u_S,u_P
		end subroutine
		
		subroutine WPENS_V(u,v,P,i,j,F_w,F_e,F_n,F_s,v_W,v_P,v_E,v_N,v_S)
			implicit none
			integer,intent(in)::i,j
			real(kind=8),intent(in),dimension(:,:) :: u,v,P
			real(kind=8),intent(out):: &
				&F_w,F_e,F_n,F_s,v_W,v_E,v_N,v_S,v_P
		end subroutine
		
		subroutine  CALCULO_RHS(&
			&phi_W,phi_P,phi_E,phi_N,phi_S,&
			&phi0_W,phi0_P,phi0_E,phi0_N,phi0_S,&
			&F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S,&
			&F0_w,F0_e,F0_n,F0_s,&
			&P_U,P_D,dx,dy,dist1,dist2,dt,Re,RHS)
			implicit none
			real(kind=8),intent(in):: &
			&phi_W,phi_P,phi_E,phi_N,phi_S,&
			&phi0_W,phi0_P,phi0_E,phi0_N,phi0_S,&
			&F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S,&
			&F0_w,F0_e,F0_n,F0_s,&
			&P_U,P_D,dx,dy,dist1,dist2,dt,Re
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
	integer	:: i,j,k,contador
	real(kind=8) :: &
		& u_W, u_P, u_E, u_N, u_S , &
		& u0_W,u0_P,u0_E,u0_N,u0_S, &
		& v_W, v_P, v_E, v_N, v_S , &
		& v0_W,v0_P,v0_E,v0_N,v0_S, &
		& F_w,F_e,F_n,F_s, &
		& P_w,P_e,P_n,P_s, &
		& F0_w,F0_e,F0_n,F0_s, &
		& P0_w,P0_e,P0_n,P0_s, &
		& RHS,alfa,beta
		
	real(kind=8),dimension(nx-1):: &
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
	
	real(kind=8) :: verif1((nx-1)*ny,(nx-1)*ny) , provi1((nx-1)*ny) ,&
					& verif2((nx-1)*(ny-1),(nx-1)*(ny-1)) , provi2((nx-1)*(ny-1))
		
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!					CALCULO RHS_U

!	-----------------------------------------------

!			CALCULO EN LA DIRECCION X

!	en cada fila...
	do i=2,ny+1
	
		contador = 0
	
!		...recorre las columnas (direccion x)
		do j=3,nx+1
		
			contador = contador + 1
			
			call WPENS_U(u_1,v_1,P,i,j,F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S)
			call WPENS_U(u_0,v_0,P,i,j,F0_w,F0_e,F0_n,F0_s,u0_W,u0_P,u0_E,u0_N,u0_S)
			P_e = P(i,j)
			P_w = P(i,j-1)
			
			call CALCULO_RHS(&
				&u_W,u_P,u_E,u_N,u_S,&
				&u0_W,u0_P,u0_E,u0_N,u0_S,&
				&F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S,&
				&F0_w,F0_e,F0_n,F0_s,&
				&P_e,P_w,dx,dy,dy,dx,dt,Re,RHS)
				
				RHSx_u(contador) = RHS
				
								
!		fin do direccion x (columna)
		end do

!			CALCULO DE dV_u EN CADA FILA	
!		diagonales de la matriz tridiagonal	
		alfa = -2._8*dt*dy/(3._8*dx*Re)
		beta =  1._8*dx*dy + 2._8*dt*dy/(3._8*dx*Re)
		diagux_u = alfa
		diagux_c = beta
		diagux_d = alfa
		diagux_c(nx-1)=beta+alfa
	
!			primer paso predictor u-> resuelve dV_u0
		call thomas(nx-1,diagux_d,diagux_c,diagux_u,dVx_u,RHSx_u)
		dV(i,3:nx+1) = dVx_u(:)
		
!	fin do direccion y (fila)
	end do

	do i=1,ny+2
		write(1,*) dV(i,:)
	end do

!	en cada columna...
	do j=3,nx+1

!			CALCULO DE dV_u EN CADA FILA	
!		diagonales de la matriz tridiagonal	
		alfa = -2._8*dt*dx/(3._8*dy*Re)
		beta =  1._8*dx*dy + 2._8*dt*dx/(3._8*dy*Re)
		diaguy_u = alfa
		diaguy_c = beta
		diaguy_d = alfa
		diaguy_c(1)=beta-alfa
		diaguy_c(ny)=beta+alfa
	
!			primer paso predictor u-> resuelve dV_u0
		call thomas(ny,diaguy_d,diaguy_c,diaguy_u,dVy_u,dV(2:ny,j))
		u_pred(2:ny+1,j) = dVy_u(:)
				
!	fin do direccion y (fila)
	end do

	do i=1,ny+2
		write(2,*) dV(i,:)
	end do

!	PENDIENTE!!!
	v_pred=0._8

return

!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end subroutine
