program T2_CDF_2017_2S


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	implicit none
	
	INTERFACE
			
		subroutine thomas(n,d,a,c,x,b)
			implicit none
			integer,intent(in)::n
			double precision,dimension(n),intent(in) :: d,a,c,b
			double precision,dimension(n),intent(out)	:: x
		end subroutine
		
		subroutine CC(nx,ny,u_0,u_1,v_0,v_1,u_init)
			implicit none
			integer,intent(in)	:: nx,ny
			real(kind=8),intent(in)::u_init
			real(kind=8),dimension(:,:),intent(inout)	:: &
				&u_0,u_1,v_0,v_1
		end subroutine
		
		subroutine RIGHT_HAND_SIDE(nx,ny,dx,dy,dt,Re,P,u_1,u_0,v_1,v_0,RHSy,RHSx)
			implicit none
			integer,intent(in)::nx,ny
			real(kind=8),intent(in)::dx,dy,dt,Re
			real(kind=8),dimension(nx+2,ny+2),intent(in):: &
				& P , u_1 , u_0 , v_1 , v_0
			real(kind=8),dimension(nx*ny),intent(out):: RHSx,RHSy
		end subroutine
		
		subroutine PREDICCION_VELOCIDAD(nx,ny,dx,dy,dt,Re,RHS_u,RHS_v,u_0,u_1,v_0,v_1)
			integer,intent(in) :: nx,ny
			real(kind=8),intent(in) :: dx,dy,dt,Re
			real(kind=8),dimension(nx*ny),intent(in) :: RHS_u,RHS_v
			real(kind=8),dimension(nx+2,ny+2),intent(inout) :: &
				&u_0,u_1,v_0,v_1
		end subroutine		
		
		subroutine TDMA_PHI(nx,ny,dx,dy,dt,u_0,u_1,v_0,v_1,phi,R)
			implicit none
			integer,intent(in) :: nx,ny
			real(kind=8),intent(in) :: dx,dy,dt
			real(kind=8),dimension((nx+2),(ny+2)),intent(in) :: &
				& u_0,u_1,v_0,v_1
			real(kind=8),intent(out) :: R
			real(kind=8),dimension((nx+2),(ny+2)),intent(out) :: & 
				& phi
		end subroutine

	END INTERFACE
	
	integer :: i,j,k,nx,ny,num_nodos,num_volumenes,contador,info,criterio2
	
	real(kind=8) :: Lx,Ly,dx,dy,u_init,rho,RHS,gama,dt,Re,criterio1,&
		& u_E,u_W,u_S,u_N,u_P,u0_E,u0_W,u0_S,u0_N,u0_P, &
		& v_E,v_W,v_S,v_N,v_P,v0_E,v0_W,v0_S,v0_N,v0_P, &
		& P_E,P_W,P_N,P_S,R

	
	real(kind=8),allocatable :: x_1(:),x_2(:),y_1(:),y_2(:),&
		& u_0(:,:),v_0(:,:),u_1(:,:),v_1(:,:),P_0(:,:),P_1(:,:),&
		& RHS_u(:),RHS_v(:),diag_d(:),phi(:,:)
	
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


!	PARAMETROS

!	DIMENSION TUBERIA
	Lx = 400._8
	Ly = 100._8
	
!	NUMERO DE VOLUMENES
!		(no considera nodos ficticios))
	nx = 40
	ny = 10
	num_volumenes = ny*nx
	
!	PASO DE TIEMPO
	dt = 1._8
	
!	NUMERO DE REYNOLDS
	Re = 100._8
	
!	VELOCIDAD ENTRADA
	u_init = 2._8
	
!	DENSIDAD
	rho = 1._8
	
!	VISCOSIDAD
!	gama = 1._8+exp(-5._8)		CUAL DE LOS DOS ES?????????????????????????????????????
	gama = 1d-5
	
!	CRITERIO RUTINA TDMA_PHI
	criterio1 = 0.5_8	! criterio1 -> max(RESTO)
	criterio2 = 10	! criterio2 -> iteraciones
	
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	
!	PREAMBULO
	
!------------------------------------------

!	VARIABLES EN LA MALLA DESPLAZADA
!		P en nodo ( x_2 , y_2 )
!		U en nodo ( x_1 , y_2 )
!		V en nodo ( x_2 , y_1 )

!	ALMACENAMIENTO DE LAS VARIABLES
!		ejemplo: componente velocidad u
!	u(1,1)	u(1,2)	u(1,3)	u(1,4)	...	u(1,n)
!	u(2,1)	u(2,2)	u(2,3)	u(2,4)	...	u(2,n)
!	u(3,1)	u(3,2)	u(3,3)	u(3,4)	...	u(3,n)
!	...		...		...		...			...

!------------------------------------------
	
!	ESPACIO INTERNODAL
	dx = Lx/dfloat(nx)
	dy = Ly/dfloat(ny)	

!	MALLA PRINCIPAL ( x_1 , y_1 )
	allocate(x_1(nx+2),y_1(ny+2))
	x_1 = (/ ( (dfloat(i)-1.0_8)*dx,i=0,nx+1 ) /)
	y_1 = (/ ( (dfloat(i)-1.0_8)*dy,i=0,ny+1 ) /) 

!	MALLA SECUNDARIA ( x_2 , y_2 )
	allocate(x_2(nx+2),y_2(ny+2))
	x_2 = (/ ( (dfloat(i)-0.5_8)*dx,i=0,nx+1 ) /)
	y_2 = (/ ( (dfloat(i)-0.5_8)*dy,i=0,ny+1 ) /)

!	VELOCIDAD INICIAL
!		en t_(n)
	allocate(u_0(nx+2,ny+2),v_0(nx+2,ny+2))
	u_0 = 0.1_8
	v_0 = 0.1_8
!		en t_(n+1)
	allocate(u_1(nx+2,ny+2),v_1(nx+2,ny+2))
	u_1 = 0.1_8
	v_1 = 0.1_8

!	CAMPO DE PRESION INICIAL
!		en t_(n)
	allocate(P_0(nx+2,ny+2))
	P_0 = 0.1_8
!		en t_(n+1)
	allocate(P_1(nx+2,ny+2))
	P_1 = 0.1_8	
	
!	FUNCION AUXILIAR PHI
	allocate(phi(nx+2,ny+2))
	phi = 0._8
	
!	IMPOSICION DE LAS CONDICIONES DE CONTORNO
	call CC(nx,ny,u_0,u_1,v_0,v_1,u_init)


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
		
				
!	PREDICCION DEL CAMPO DE VELOCIDAD NO SOLENOIDAL

	allocate(RHS_u(nx*ny),RHS_v(nx*ny))
		
	call RIGHT_HAND_SIDE(nx,ny,dx,dy,dt,Re,P_0,u_1,u_0,v_1,v_0,RHS_u,RHS_v)

	call PREDICCION_VELOCIDAD(nx,ny,dx,dy,dt,Re,RHS_u,RHS_v,u_0,u_1,v_0,v_1)
		
		
!!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	
!	RESOLUCION DE LA ECUACION DE POISSON PARA LA PRESION

	contador = 0
	info = 0
	
	do while ( info .eq. 0 ) 
		
		contador = contador + 1
		
		call TDMA_PHI(nx,ny,dx,dy,dt,u_0,u_1,v_0,v_1,phi,R)
		
		if ( R .le. criterio1 .or. contador .ge. criterio2 ) then
			info = 1
		end if
	
	end do
	

!!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!	CORRECCION DE LA PRESION
!	do i=1,nx+1
!		do j=1,ny+1
!			P_1(i,j) = (P_0(i,j)+phi(i,j))/(dx*dy) &
!				& (u_1(i,j)-u_1(i,j+1))*dx &
!				& (v_1(i-1,j)-v_1(i,j))*dy
!		end do
!	end do

!!	CORRECCION DE LA VELOCIDAD
!	u_0 = v_1
!	v_0 = v_1
!	do i=1,nx+1
!		do j=1,ny+1
!			u_1(i,j) = u_0 - dt*dx*(phi(i,j)-phi(i,j-1))
!			v_1(i,j) = v_0 - dt*dy*(phi(i+1,j)-phi(i,j))
!		end do
!	end do

!!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



end program
