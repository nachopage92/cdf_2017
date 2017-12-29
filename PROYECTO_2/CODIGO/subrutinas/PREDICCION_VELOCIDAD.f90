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
			&P_U,P_D,dx,dy,dist,dt,Re,RHS)
			implicit none
			real(kind=8),intent(in):: &
			&phi_W,phi_P,phi_E,phi_N,phi_S,&
			&phi0_W,phi0_P,phi0_E,phi0_N,phi0_S,&
			&F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S,&
			&F0_w,F0_e,F0_n,F0_s,&
			&P_U,P_D,dx,dy,dist,dt,Re
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
		
	real(kind=8),dimension((nx-1)*ny):: &
		& RHS_u,diagu_d,diagu_c,diagu_u,dV_u0,dV_u
	
	real(kind=8),dimension((nx-1)*(ny-1)):: &
		& RHS_v,diagv_d,diagv_c,diagv_u,dV_v0,dV_v
		
	real(kind=8),dimension(ny+2,nx+2) :: dV
		
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!					CALCULO RHS_U

!	-----------------------------------------------

!	recorre en direccion y (fila)
	contador = 0
	do i=2,ny+1
!		recorre en direccion x (columna)
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
				&P_e,P_w,dx,dy,dx,dt,Re,RHS)
				
				RHS_u(contador) = -RHS	
				
!		fin do direccion x (columna)
		end do
!	fin do direccion y (fila)
	end do



	dV = 0._8
	contador = 0
	do i=2,ny+1
		do j=3,nx+1
			contador = contador + 1
			dV(i,j) = RHS_u(contador) 
		end do
	end do
	
	do i=1,ny+2
		write(100,*) dV(i,:)
	end do





!	-----------------------------------------------


!			CALCULO DE dV_u

!		direccion x
	alfa = -2._8*dt*dy/(3._8*dx*Re)
	beta =  1._8*dx*dy + 2._8*dt*dy/(3._8*dx*Re)
	do i=1,(nx-1)*ny
		if ( mod(i,(nx-1)) .eq. 0 ) then
			diagu_u(i) = 0._8
			diagu_c(i) = alfa + beta
			diagu_d(i) = 0._8
		else
			diagu_u(i) = alfa
			diagu_c(i) = beta
			diagu_d(i) = alfa
		end if
	end do
!		primer paso predictor u-> resuelve dV_u0
	call thomas((nx-1)*ny,diagu_d,diagu_c,diagu_u,dV_u0,RHS_u)


	dV = 0._8
	contador = 0
	do i=2,ny+1
		do j=3,nx+1
			contador = contador + 1
			dV(i,j) = dV_u0(contador) 
		end do
	end do
	
	do i=1,ny+2
		write(200,*) dV(i,:)
	end do

!	-----------------------------------------------


!		traspaso dV_u0 
!			(vector en direccion x)
!			 	a dV (matriz)
	dV = 0._8
	contador = 0
	do i=2,ny+1
		do j=3,nx+1
			contador = contador + 1
			dV(i,j) = dV_u0(contador) 
		end do
	end do
!		traspaso dV (matriz) 
!			a dV_u0
!				(vector en direccion y)
	contador = 0
	do j=3,nx+1
!		recorre en direccion x (columna)
		do i=2,ny+1
			contador = contador + 1
			dV_u0(contador) = dV(i,j)
		end do
	end do
	
	
!	-----------------------------------------------
	
!		direccion y
	alfa = -2._8*dt*dx/(3._8*dy*Re)
	beta = 1._8*dx*dy + 4._8*dt*dx/(3._8*dy*Re)
	do i=1,(nx-1)*ny
		if ( i .eq. 1 .or. mod(i,(ny+1)) .eq. 0 ) then
			diagu_u(i) = alfa
			diagu_c(i) = beta - alfa
			diagu_d(i) = alfa
		else if ( mod(i,ny) .eq. 0 ) then
			diagu_u(i) = 0._8
			diagu_c(i) = beta + alfa
			diagu_d(i) = 0._8
		else
			diagu_u(i) = alfa
			diagu_c(i) = beta
			diagu_d(i) = alfa
		end if
	end do
	
!		segundo paso predictor -> resuelve dV_u
	call thomas((nx-1)*ny,diagu_d,diagu_c,diagu_u,dV_u,dV_u0)



	dV = 0._8
	contador = 0
	do j=3,nx+1
		do i=2,ny+1 
			contador = contador + 1
			dV(i,j) = dV_u(contador) 
		end do
	end do
	
	do i=1,ny+2
		write(300,*) dV(i,:)
	end do


!	-----------------------------------------------

!			PREDICTOR DEL CAMPO DE VELOCIDADES
!				traspaso variables de vector a matriz

	contador = 0
	do j=3,nx+1
		do i=2,ny+1
			contador = contador + 1
			u_pred(i,j) = u_1(i,j) + dV_u(contador)
		end do
	end do
	
	do i=1,ny+2
		write(400,*) u_pred(i,:)
	end do
	
!	return
	
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!					CALCULO RHS_V

!	-----------------------------------------------

!	recorre en direccion y (fila)
	contador = 0
	do i=2,ny
!		recorre en direccion x (columna)
		do j=3,nx+1
		
			contador = contador + 1
			
			call WPENS_V(u_1,v_1,P,i,j,F_w,F_e,F_n,F_s,v_W,v_P,v_E,v_N,v_S)
			call WPENS_V(u_0,v_0,P,i,j,F0_w,F0_e,F0_n,F0_s,v0_W,v0_P,v0_E,v0_N,v0_S)
			!gradiente de presion
			P_n = P(i+1,j)
			P_s = P(i,j)
			
			call CALCULO_RHS(&
				&v_W,v_P,v_E,v_N,v_S,&
				&v0_W,v0_P,v0_E,v0_N,v0_S,&
				&F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S,&
				&F0_w,F0_e,F0_n,F0_s,&
				&P_n,P_s,dx,dy,dy,dt,Re,RHS)
				
				RHS_v(contador) = -RHS	
				
!		fin do direccion x (columna)
		end do
!	fin do direccion y (fila)
	end do
	
!	-----------------------------------------------

!			CALCULO DE dV_v

!		direccion x
	alfa = -2._8*dt*dy/(3._8*dx*Re)
	beta = 1._8*dx*dy + 4._8*dt*dy/(3._8*dx*Re)
	do i=1,(nx-1)*(ny-1)
		if ( mod(i,(nx-1)) .eq. 0 ) then
			diagv_u(i) = 0._8
			diagv_c(i) = alfa + beta
			diagv_d(i) = 0._8
		else
			diagv_u(i) = alfa
			diagv_c(i) = beta
			diagv_d(i) = alfa
		end if	
	end do
		
!		primer paso predictor u-> resuelve dV_u0
	call thomas((nx-1)*(ny-1),diagv_d,diagv_c,diagv_u,dV_v0,RHS_v)
	
!	-----------------------------------------------


!		traspaso dV_u0 
!			(vector en direccion x)
!			 	a dV (matriz)
	dV = 0._8
	contador = 0
	do i=2,ny
		do j=3,nx+1
			contador = contador + 1
			dV(i,j) = dV_v0(contador) 
		end do
	end do
!		traspaso dV (matriz)
!			a dV_u0
!			 	(vector en direccion y)
	contador = 0
	do j=3,nx+1
		do i=2,ny
			contador = contador + 1
			dV_v0(contador) = dV(i,j)
		end do
	end do


!	-----------------------------------------------


!		direccion y
	alfa = -2._8*dt*dx/(3._8*dy*Re)
	beta = 1._8*dx*dy + 4._8*dt*dx/(3._8*dy*Re)
	do i=1,(nx-1)*(ny-1)
		if ( mod(i,(ny-1)) .eq. 0 ) then
			diagv_u(i) = 0._8
			diagv_c(i) = beta
			diagv_d(i) = 0._8
		else
			diagv_u(i) = alfa
			diagv_c(i) = beta
			diagv_d(i) = alfa
		end if
	end do
!		segundo paso predictor -> resuelve dV_u
	call thomas((nx-1)*(ny-1),diagv_d,diagv_c,diagv_u,dV_v,dV_v0)

!	-----------------------------------------------

!			PREDICTOR DEL CAMPO DE VELOCIDADES
!				traspaso variables de vector a matriz

	contador = 0
	do j=3,nx+1
		do i=2,ny
			contador = contador + 1
			v_pred(i,j) = v_1(i,j) + dV_v(contador)
		end do
	end do
	

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end subroutine
