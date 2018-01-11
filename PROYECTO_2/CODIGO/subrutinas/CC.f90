subroutine CC(nx,ny,u,v,u_init)
	implicit none
	integer,intent(in)	:: nx,ny
	real(kind=8),intent(in)	:: u_init
	real(kind=8),dimension(:,:),intent(inout) :: u,v

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
	u(:,nx+2) = u_init
	u(:,ny+1) = u_init
	u(:,1) = u(:,ny+1)
	u(:,2) = u(:,nx+2)
!			v(0,y) = v(Lx,y)
!			dv/dx(0,y) = dv/dx(Lx,y)
	v(:,nx+2) = 0._8
	v(:,ny+1) = 0._8
	v(:,1) = v(:,ny+1)
	v(:,2) = v(:,nx+2)
	
	
end subroutine
