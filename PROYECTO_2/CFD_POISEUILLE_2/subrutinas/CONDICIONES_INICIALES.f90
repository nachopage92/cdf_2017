subroutine CONDICIONES_INICIALES(u,v)

	use variables
	implicit none	
	
	real(kind=8),dimension(ny+2,nx+2),intent(inout)::u,v
	
!		condicion inicial para u
	u = 0._8
	u(2:ny+1,nx+1) = u_init
	u(2:ny+1,nx+2) = u_init
	
!		condicion inicial para v
	v = 0._8
	
end subroutine
