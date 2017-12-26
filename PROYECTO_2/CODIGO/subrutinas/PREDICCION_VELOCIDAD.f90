subroutine PREDICCION_VELOCIDAD(nx,ny,dx,dy,dt,Re,RHS_u,RHS_v,u_0,u_1,v_0,v_1,u_pred,v_pred)

	implicit none
	
	INTERFACE

		subroutine thomas(n,d,a,c,x,b)
			implicit none
			integer,intent(in)::n
			double precision,dimension(n),intent(in) :: d,a,c,b
			double precision,dimension(n),intent(out)	:: x
		end subroutine

	END INTERFACE

	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy,dt,Re
	real(kind=8),dimension(nx*ny),intent(in) :: RHS_u,RHS_v
	real(kind=8),dimension(ny+2,nx+2),intent(in) :: &
		&u_0,u_1,v_0,v_1
	real(kind=8),dimension(ny+2,nx+2),intent(out) :: &
		& u_pred,v_pred
	
	!variables locales
	integer :: i,j,contador
	real(kind=8),dimension(nx*ny) :: diag_d,diag_c,diag_u
	real(kind=8),dimension(nx*ny) :: dV_u0,dV_v0,dV_u,dV_v
	

!	CALCULO DEL SEGUNDO PASO (EN DIRECCION X)
!		vectores de las diagonales de la matriz
	diag_d = 2._8*dt*dy/(3._8*Re)
	diag_c = dx*dy - 4._8*dt*dy/(3._8*Re)
	diag_u = 2._8*dt*dy/(3._8*Re)
!		primer paso predictor u-> resuelve dV_u0
	call thomas(nx*ny,diag_d,diag_c,diag_u,dV_u0,RHS_u)
!		primer paso predictor v-> resuelve dV_v0
	call thomas(nx*ny,diag_d,diag_c,diag_u,dV_v0,RHS_v)

!	CALCULO DEL SEGUNDO PASO (EN DIRECCION Y)
!		vectores de las diagonales de la matriz
	diag_d = 2._8*dt*dx/(3._8*Re)
	diag_c = dx*dy - 4._8*dt*dx/(3._8*Re)
	diag_u = 2._8*dt*dx/(3._8*Re)
!		segundo paso predictor -> resuelve dV_u1
	call thomas(nx*ny,diag_d,diag_c,diag_u,dV_u,dV_u0)
!		segundo paso predictor -> resuelve dV_v1
	call thomas(nx*ny,diag_d,diag_c,diag_u,dV_v,dV_v0)
		
!	PREDICTOR DEL CAMPO DE VELOCIDADES
	contador = 0
	do i=2,ny+1
		do j=2,nx+1
			contador = contador + 1
!			traspaso variables de vector a matriz
			u_pred(i,j) = u_1(i,j) + dV_u(contador)
			v_pred(i,j) = v_1(i,j) + dV_v(contador)
		end do
	end do
	
end subroutine
