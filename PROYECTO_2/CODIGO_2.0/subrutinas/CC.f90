subroutine CC(u,v)
	
	use variables

	implicit none
	real(kind=8),dimension(ny+2,nx+2),intent(inout) :: u,v

! CONDICIONES DE CONTORNO

!		cara inferior ( no-slip )
!			u(x,0) = 0
	u(1,:) = -u(2,:)
!			v(x,0) = 0	
	v(1,:) = 0._8

!		cara superior ( symmetry )
!			du/dy = 0
	u(ny+2,:) = u(ny+1,:)
!			v(x,y=h) = 0	
	v(ny+1,:) = 0._8
!			(no entrega informacion)
	v(ny+2,:) = 0._8

!		condicion de flujo periodico
!			u(0,y) = u(Lx,y)
!			du/dx(0,y) = du/dx(Lx,y)
	u(:,1) = u(:,nx+1)
	u(:,nx+2) = u(:,2)
!			v(0,y) = v(Lx,y)
!			dv/dx(0,y) = dv/dx(Lx,y)
	v(:,1) = v(:,nx+1)
	v(:,nx+2) = v(:,2)

	
end subroutine
