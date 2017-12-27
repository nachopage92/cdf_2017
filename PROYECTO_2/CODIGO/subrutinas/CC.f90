subroutine CC(nx,ny,u_0,u_1,v_0,v_1,phi,u_init)
	implicit none
	integer,intent(in)	:: nx,ny
	real(kind=8),intent(in)	:: u_init
	real(kind=8),dimension(:,:),intent(inout)	:: &
		&u_0,u_1,v_0,v_1,phi

! CONDICIONES DE CONTORNO

!		cara lateral izquierda ( u = u_init , v = 0 )
!			v(0,y) = 0
	v_0(:,1) = 0._8
	v_1(:,1) = 0._8
	v_0(:,2) = 0._8
	v_1(:,2) = 0._8
!			u(0,y) = u_init
	u_0(:,1) = u_init
	u_1(:,1) = u_init
!			(no entrega informacion)
	u_0(:,2) = u_init
	u_1(:,2) = u_init

	phi(:,1) = 0._8
	phi(:,2) = 0._8
	
	
!		cara inferior ( u,v = 0 )
!			u(x,0) = 0
	u_0(1,:) = -u_0(2,:)
	u_1(1,:) = -u_1(2,:)
!			v(x,0) = 0	
	v_0(1,:) = 0._8
	v_1(1,:) = 0._8
			
	phi(1,:) = 0._8
	phi(2,:) = 0._8


!		cara lateral derecha ( du/dx = 0 , dv/dx = 0 )
!			du/dx(x=L) = 0
	u_0(:,nx+2) = u_0(:,nx+1)
	u_1(:,nx+2) = u_1(:,nx+1)
!!!			dv/dx(x=L) = 0
!!	v_0(:,nx+2) = v_0(:,nx+1)
!!	v_1(:,nx+2) = v_1(:,nx+1)
!			v(x=L) = 0
	v_0(:,nx+2) = 0._8
	v_1(:,nx+2) = 0._8
	v_0(:,nx+1) = 0._8
	v_1(:,nx+1) = 0._8
	
	phi(:,nx+2) = 0._8
	phi(:,nx+1) = 0._8


!		cara superior ( simetria )
!			du/dy = 0
	u_0(ny+2,:) = u_0(ny+1,:)
	u_1(ny+2,:) = u_1(ny+1,:)
	
!			v(x,y=h/2) = 0	
	v_0(ny+1,:) = 0._8
	v_1(ny+1,:) = 0._8
!			(no entrega informacion)
	v_0(ny+2,:) = 0._8
	v_1(ny+2,:) = 0._8
	
	phi(ny+2,:) = 0._8
	phi(ny+1,:) = 0._8

end subroutine
