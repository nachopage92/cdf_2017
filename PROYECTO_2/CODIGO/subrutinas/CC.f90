subroutine CC(nx,ny,u_0,u_1,v_0,v_1,u_init)
	implicit none
	integer,intent(in)	:: nx,ny
	real(kind=8),intent(in)	:: u_init
	real(kind=8),dimension(:,:),intent(inout)	:: &
		&u_0,u_1,v_0,v_1

! CONDICIONES DE CONTORNO PARA u (NODOS FICTICIOS)

!		cara lateral izquierda ( u = u_init , v = 0 )
	u_0(:,2) = u_init
	u_1(:,2) = u_init
	u_0(:,1) = u_0(:,2)
	u_1(:,1) = u_1(:,2)	
	
	v_0(:,1) = -v_0(:,2)
	v_1(:,1) = -v_1(:,2)	

!		cara superior ( u,v = 0 )
	u_0(1,:) = -u_0(2,:)
	u_1(1,:) = -u_1(2,:)
	v_0(2,:) = 0._8
	v_1(2,:) = 0._8
	
	v_0(1,:) = v_1(2,:)
	v_1(1,:) = v_0(2,:)


!		cara lateral derecha ( du/dx = 0 , dv/dx = 0 )
	u_0(:,nx+2) = u_0(:,nx+1)
	u_1(:,nx+2) = u_1(:,nx+1)
	
	v_0(:,nx+2) = v_0(:,nx+1)
	v_1(:,nx+2) = v_1(:,nx+1)
	
!		cara inferior ( simetria )
	u_0(ny+2,:) = u_0(ny+1,:)
	u_1(ny+2,:) = u_1(ny+1,:)
	
	v_0(ny+2,:) = 0._8
	v_1(ny+2,:) = 0._8

end subroutine
