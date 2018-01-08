subroutine CC(nx,ny,u,v,u_init)
	implicit none
	integer,intent(in)	:: nx,ny
	real(kind=8),intent(in)	:: u_init
	real(kind=8),dimension(:,:),intent(inout) :: u,v

! CONDICIONES DE CONTORNO

!		cara lateral izquierda ( u = u_init , v = 0 )
!			v(0,y) = 0
	v(:,1) = 0._8
	v(:,2) = 0._8
!			u(0,y) = u_init
	u(:,1) = u_init
!			(no entrega informacion)
	u(:,2) = u_init
	
!		cara inferior ( u,v = 0 )
!			u(x,0) = 0
	u(1,:) = -u(2,:)
!			v(x,0) = 0	
	v(1,:) = 0._8

!		cara lateral derecha ( du/dx = 0 , dv/dx = 0 )
!			du/dx(x=L) = 0
	u(:,nx+2) = u(:,nx+1)
!			dv/dx(x=L) = 0
	v(:,nx+2) = v(:,nx+1)

!		cara superior ( u,v = 0 )
!			du/dy = 0
	u(ny+2,:) = -u(ny+1,:)
!			v(x,y=h/2) = 0	
	v(ny+1,:) = 0._8
!			(no entrega informacion)
	v(ny+2,:) = 0._8


end subroutine
