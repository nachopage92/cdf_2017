
subroutine  CALCULO_RHS(u,v,u0,v0,P,i,j,str,RHS)

!!:::::::::::::::::::::::::::::::::::::::::::::
	use variables
	implicit none
	!entrada
	integer,intent(in) :: i,j
	real(kind=8),dimension(ny+2,nx+2),intent(in) ::&
		& u,v,u0,v0,P
	character(len=1),intent(in) :: str
	!salida
	real(kind=8),intent(out)::RHS
	
	!locales
	real(kind=8)	::	H,H0,G,dP, &
	&	F_w,F_e,F_n,F_s, &
	&	F0_w,F0_e,F0_n,F0_s, &
	&	phi_W,phi_P,phi_E,phi_N,phi_S, &
	&	phi0_W,phi0_P,phi0_E,phi0_N,phi0_S 
	
!:::::::::::::::::::::::::::::::::::::::::::::
				
	call WPENS(str,u,v,i,j,F_w,F_e,F_n,F_s,&
		&phi_W,phi_P,phi_E,phi_N,phi_S)
	call WPENS(str,u0,v0,i,j,F0_w,F0_e,F0_n,F0_s,&
		&phi0_W,phi0_P,phi0_E,phi0_N,phi0_S)

	if ( str .eq. 'u' ) then
		dP = (P(i,j)-P(i,j-1))*dy
	else if ( str .eq. 'v' ) then
		dP = (P(i+1,j)-P(i,j))*dx
	else
		print*, 'ERROR'
		RETURN
	end if
				
!:::::::::::::::::::::::::::::::::::::::::::::

!			DISCRETIZACION DE LA CONVECCION

	H = ((phi_E+phi_P)*F_e-(phi_P+phi_W)*F_w)*0.5_8*dy + &
	&	((phi_N+phi_P)*F_n-(phi_P+phi_S)*F_s)*0.5_8*dx
		
	H0 = ((phi0_E+phi0_P)*F0_e-(phi0_P+phi0_W)*F0_w)*0.5_8*dy + &
	&	 ((phi0_N+phi0_P)*F0_n-(phi0_P+phi0_S)*F0_s)*0.5_8*dx	
		
	G = (phi_E-2._8*phi_P+phi_W)*(dy/dx)+(phi_N-2._8*phi_P+phi_S)*(dx/dy)

!			RIGHT HAND SIDE

	RHS = ( phi_P - phi0_P )*(dy*dx)/3._8 &
		& - (2._8*dt/3._8) * ( 2._8*H - H0 + dP - G )
	
!:::::::::::::::::::::::::::::::::::::::::::::	

	contains
	
	subroutine WPENS(str,u,v,i,j,F_w,F_e,F_n,F_s,&
	 &phi_W,phi_P,phi_E,phi_N,phi_S)
		implicit none
		character(len=1)::str
		integer,intent(in)::i,j
		real(kind=8),intent(in),dimension(:,:) :: u,v
		real(kind=8),intent(out):: &
			&F_w,F_e,F_n,F_s,phi_W,phi_E,phi_N,phi_S,phi_P
		if ( str .eq. 'u' ) then
			F_e = 0.5_8 * ( u(i,j+1) + u(i,j) )
			F_w = 0.5_8 * ( u(i,j) + u(i,j-1) )
			F_n = 0.5_8 * ( v(i,j) + v(i,j-1) ) 
			F_s = 0.5_8 * ( v(i-1,j) + v(i-1,j-1) )
			phi_W = u(i,j-1)
			phi_E = u(i,j+1)
			phi_N = u(i+1,j)
			phi_S = u(i-1,j)
			phi_P = u(i,j) 
		else if ( str .eq. 'v' ) then
			F_e = 0.5_8 * ( u(i,j+1) + u(i+1,j+1) )
			F_w = 0.5_8 * ( u(i,j) + u(i+1,j) )
			F_n = 0.5_8 * ( v(i+1,j) + v(i,j) ) 
			F_s = 0.5_8 * ( v(i,j) + v(i-1,j) )
			phi_W = v(i,j-1)
			phi_E = v(i,j+1)
			phi_N = v(i+1,j)
			phi_S = v(i-1,j)
			phi_P = v(i,j)
		end if
	end subroutine
		
!:::::::::::::::::::::::::::::::::::::::::::::	

end subroutine

