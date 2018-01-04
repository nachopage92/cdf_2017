!:::::::::::::::::::::::::::::::::::::::::::::

! CALCULO DE RHS PARA CALCULO DE PHI

!:::::::::::::::::::::::::::::::::::::::::::::

subroutine  CALCULO_RHS(&
	&phi_W,phi_P,phi_E,phi_N,phi_S,&
	&phi0_W,phi0_P,phi0_E,phi0_N,phi0_S,&
	&F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S,&
	&F0_w,F0_e,F0_n,F0_s,&
	&P_U,P_D,dx,dy,dist1,dist2,dt,Re,RHS)

!:::::::::::::::::::::::::::::::::::::::::::::

	implicit none

	real(kind=8),intent(in):: &
	&phi_W,phi_P,phi_E,phi_N,phi_S,&
	&phi0_W,phi0_P,phi0_E,phi0_N,phi0_S,&
	&F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S,&
	&F0_w,F0_e,F0_n,F0_s,&
	&P_U,P_D,dx,dy,dist1,dist2,dt,Re


	real(kind=8),intent(out) :: RHS

	!VARIABLES LOCALES
	real(kind=8) :: & 		
		& Pe_w,Pe_e,Pe_n,Pe_s, &
		& h_W,h_E,h_S,h_N,h_P, &
		& h0_W,h0_E,h0_S,h0_N,h0_P, &
		& g_W,g_E,g_S,g_N,g_P, &
		& H,H0,G,dP
	
!:::::::::::::::::::::::::::::::::::::::::::::
						
!			COEFICIENTES DE LA CONVECCION	
		h_W = -0.5_8*dy*F_w
		h_E =  0.5_8*dy*F_e
		h_S = -0.5_8*dx*F_s
		h_N =  0.5_8*dx*F_n
		h_P = - h_W - h_E - h_S - h_N
		
		h0_W = -0.5_8*dy*F0_w
		h0_E =  0.5_8*dy*F0_e
		h0_S = -0.5_8*dx*F0_s
		h0_N =  0.5_8*dx*F0_n
		h0_P = - h0_W - h0_E - h0_S - h0_N 
		
		
!			CALCULO DE LA DIFUSION
		g_w = (-1._8/Re) * (dy/dx)
		g_e = (-1._8/Re) * (dy/dx)
		g_n = (-1._8/Re) * (dx/dy)
		g_s = (-1._8/Re) * (dx/dy)
		g_p = -( g_w + g_e + g_n + g_s )
	
				
!:::::::::::::::::::::::::::::::::::::::::::::

!			DISCRETIZACION DE LA CONVECCION
		H = h_W*phi_W + h_E*phi_E + &
			& h_S*phi_S + h_N*phi_N + h_P*phi_P
		H0 = h0_W*phi0_W + h0_E*phi0_E + &
			&h0_S*phi0_S + h0_N*phi0_N + h0_P*phi0_P

!			DISCRETIZACION DE LA DIFUSION
		G = g_W*phi_W + g_E*phi_E + &
			&g_S*phi_S + g_N*phi_N + g_P*phi_P

!			DISCRETIZACION DEL GRADIENTE DE PRESION
		dP = (P_U-P_D)*dist1
		
!			RIGHT HAND SIDE
		RHS = ( phi_P - phi0_P )*dy*dx/3._8 &
			& - 2._8*dt/3._8 * ( 2._8*H - H0 + dP - G )


!print*, phi0_w , phi0_p , phi0_e , phi0_n , phi0_s
!print*, phi_w , phi_p , phi_e , phi_n , phi_s
!print*, h0_w , h0_p , h0_e , h0_n , h0_s
!print*, h_w , h_p , h_e , h_n , h_s
!print*, g_w , g_p , g_e , g_n , g_s
!print*, dP
!print*, RHS
!print*, 

!!			CALCULO DEL NUMERO DE PECLET
!		Pe_w = abs(F_w)*Re*dx
!		Pe_e = abs(F_e)*Re*dx
!		Pe_n = abs(F_n)*Re*dy
!		Pe_s = abs(F_s)*Re*dy
		
		
!:::::::::::::::::::::::::::::::::::::::::::::		

end subroutine
