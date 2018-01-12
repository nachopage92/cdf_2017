!:::::::::::::::::::::::::::::::::::::::::::::

!	CALCULO DE RHS UTILIZANDO DISCRETIZACION
!	PROPUESTA POR PATANKAR

!:::::::::::::::::::::::::::::::::::::::::::::

subroutine  CALCULO_RHS(&
		& u_W,u_E,u_P,v_N,v_S,v_P, &
		& u0_W,u0_E,u0_P,v0_N,v0_S,v0_P, &
		& phi_W,phi_E,phi_S,phi_N,phi_P,&
		& phi0_W,phi0_E,phi0_S,phi0_N,phi0_P,&
		& P_1,P_2,dx,dy,dist,dt,Re,RHS)

!:::::::::::::::::::::::::::::::::::::::::::::

	implicit none

	!VELOCIDAD ( u , v )
	real(kind=8),intent(in) :: &
		& u_W,u_E,u_P,v_S,v_N,v_P, &				
		& u0_W,u0_E,u0_P,v0_S,v0_N,v0_P
	
	!ESCALAR PASIVO phi
	real(kind=8),intent(in) :: &
		& phi_W,phi_E,phi_S,phi_N,phi_P, &
		& phi0_W,phi0_E,phi0_S,phi0_N,phi0_P
		
	!PRESION ( P_E , P_W ) o ( P_N , P_S )
	real(kind=8),intent(in) :: P_1,P_2		
	
	!DISCRETIZACION ESPACIO-TIEMPO
	real(kind=8),intent(in) :: dx,dy,dist,dt

	!NUMERO DE REYNOLDS
	real(kind=8),intent(in)::Re
	
	!RHS (right hand side)
	real(kind=8),intent(out) :: RHS

	!VARIABLES LOCALES
	real(kind=8) :: & 		
		& F_w,F_e,F_n,F_s,F_p, &
		& F0_w,F0_e,F0_n,F0_s,F0_p, &
		& D_w,D_e,D_n,D_s,D_p, &
		& Pe_w,Pe_e,Pe_n,Pe_s, &
		& h_W,h_E,h_S,h_N,h_P, &
		& h0_W,h0_E,h0_S,h0_N,h0_P, &
		& g_W,g_E,g_S,g_N,g_P, &
		& H,H0,G,dP

!:::::::::::::::::::::::::::::::::::::::::::::

!			CALCULO DE LA CONVECCION 
		F_w = 0.5_8*( u_W + u_P )
		F_e = 0.5_8*( u_P + u_E )
		F_n = 0.5_8*( v_N + v_P )
		F_s = 0.5_8*( v_P + v_S )
		F_p = -( F_w + F_e + F_n + F_s )
		
		F0_w = 0.5_8*( u0_W + u0_P )
		F0_e = 0.5_8*( u0_P + u0_E )
		F0_n = 0.5_8*( v0_N + v0_P )
		F0_s = 0.5_8*( v0_P + v0_S )
		F0_p = -( F0_w + F0_e + F0_n + F0_s )
				
!:::::::::::::::::::::::::::::::::::::::::::::
		
!			CALCULO DE LA DIFUSION
		D_w = dy/Re
		D_e = dy/Re
		D_n = dx/Re
		D_s = dx/Re
		D_p = -( D_w + D_e + D_n + D_s )
					
!:::::::::::::::::::::::::::::::::::::::::::::
			
!			CALCULO DEL NUMERO DE PECLET
		Pe_w = F_w/D_w
		Pe_e = F_e/D_e
		Pe_n = F_n/D_n
		Pe_s = F_s/D_s
			
!:::::::::::::::::::::::::::::::::::::::::::::
		
!			COEFICIENTES DE LA CONVECCION	
		h_W = max( F_w,0._8)
		h_E = max(-F_e,0._8)
		h_S = max( F_s,0._8)
		h_N = max(-F_n,0._8)
		h_P = -( h_E + h_W + h_S + h_N )
		
		h0_W = max( F0_w,0._8)
		h0_E = max(-F0_e,0._8)
		h0_S = max( F0_s,0._8)
		h0_N = max(-F0_n,0._8)
		h0_P = -( h0_E + h0_W + h0_S + h0_N )
			
!:::::::::::::::::::::::::::::::::::::::::::::
			
!			COEFICIENTES DE LA DIFUSION
		g_W = max(0._8,( 1._8 - 0.1_8*abs(Pe_w) )**5._8)
		g_E = max(0._8,( 1._8 - 0.1_8*abs(Pe_w) )**5._8)
		g_S = max(0._8,( 1._8 - 0.1_8*abs(Pe_w) )**5._8)
		g_N = max(0._8,( 1._8 - 0.1_8*abs(Pe_w) )**5._8)
		g_p = -( g_W + g_E + g_S + g_N )
		
!:::::::::::::::::::::::::::::::::::::::::::::

!			DISCRETIZACION DE LA DIFUSION
		H = h_W*phi_W + h_E*phi_E + &
			& h_S*phi_S + h_N*phi_N + h_P*phi_P
		H0 = h0_W*phi0_W + h0_E*phi0_E + &
			&h0_S*phi0_S + h0_N*phi0_N + h0_P*phi0_P
!			DISCRETIZACION DE LA DIFUSION
		G = g_W*phi_W + g_W*phi_E + &
			&g_S*phi_S + g_N*phi_N + g_P*phi_P
!			DISCRETIZACION DEL GRADIENTE DE PRESION
		dP = (P_1-P_2)/dist
		
!			RIGHT HAND SIDE
		RHS = (phi_P-phi0_P)/3._8 - (2*dt)/3._8 * (2._8*H-H0+dP-G)

end subroutine
