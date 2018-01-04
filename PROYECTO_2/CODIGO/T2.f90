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
		
		subroutine CC(nx,ny,u,v,u_init)
			implicit none
			integer,intent(in)	:: nx,ny
			real(kind=8),intent(in)	:: u_init
			real(kind=8),dimension(:,:),intent(inout) :: u,v
		end subroutine
		
		subroutine PREDICCION_VELOCIDAD(nx,ny,dx,dy,dt,Re,P,u_1,u_0,v_1,v_0,u_pred,v_pred)
			implicit none
			integer,intent(in)::nx,ny
			real(kind=8),intent(in)::dx,dy,dt,Re
			real(kind=8),dimension(ny+2,nx+2),intent(in):: &
				& P , u_1 , u_0 , v_1 , v_0
			real(kind=8),dimension(ny+2,nx+2),intent(in):: &
				& u_pred , v_pred
		end subroutine
				
		subroutine TDMA_PHI(nx,ny,dx,dy,dt,u,v,phi,R)
			integer,intent(in) :: nx,ny
			real(kind=8),intent(in) :: dx,dy,dt
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: &
				& u,v
			real(kind=8),intent(out) :: R
			real(kind=8),dimension(ny+2,nx+2),intent(inout) :: & 
				& phi
		end subroutine
		
		subroutine CORRECCION_PRESION(nx,ny,dx,dy,Re,P,phi,u,v)
			implicit none
			integer,intent(in) :: nx,ny
			real(kind=8),intent(in) :: dx,dy,Re
			real(kind=8),dimension((nx+2),(ny+2)),intent(in) :: phi,u,v
			real(kind=8),dimension((nx+2),(ny+2)),intent(inout) :: P
		end subroutine

		subroutine CORRECCION_VELOCIDAD(nx,ny,dx,dy,dt,u_pred,v_pred,phi,u_0,u_1,v_0,v_1)
			implicit none
			integer,intent(in) :: nx,ny
			real(kind=8),intent(in) :: dx,dy,dt
			real(kind=8),dimension((nx+2),(ny+2)),intent(in) :: phi
			real(kind=8),dimension((nx+2),(ny+2)),intent(in) :: u_pred,v_pred
			real(kind=8),dimension((nx+2),(ny+2)),intent(inout) :: u_0,u_1,v_0,v_1
		end subroutine
		
		subroutine CALCULO_CONVERGENCIA(nx,ny,dx,dy,dt,u,v,phi,CDM)
			implicit none
			integer,intent(in) :: nx,ny
			real(kind=8),intent(in) :: dx,dy,dt
			real(kind=8),dimension((nx+2),(ny+2)),intent(in) :: u,v,phi
			real(kind=8),intent(out):: CDM
		end subroutine

	END INTERFACE
	
	integer :: i,j,k,nx,ny,nt,num_nodos,num_volumenes,contador,niter,&
		&info_TDMA,info_ITER,criterio2_TDMA,criterio2_convergencia
	
	real(kind=8) :: Lx,Ly,dx,dy,u_init,rho,RHS,gama,dt,Re,&
		& criterio1_TDMA,T,criterio1_convergencia,R,CFL,CDM,aux
	
	real(kind=8),allocatable :: x_1(:),x_2(:),y_1(:),y_2(:),&
		& u_0(:,:),v_0(:,:),u_1(:,:),v_1(:,:),P(:,:),&
		& RHS_u(:),RHS_v(:),phi(:,:),&
		& u_pred(:,:),v_pred(:,:)
		
	character(len=3) :: dummy = 'cfd' 
	
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


!	PARAMETROS

!	DIMENSION TUBERIA
	Lx = 2._8
	Ly = 1._8
	
!	NUMERO DE VOLUMENES
!		(no considera nodos ficticios))
	nx = 20
	ny = 10
	num_volumenes = ny*nx
	
!!	PASO DE TIEMPO
!	dt = 0.1_8
!	T = 10._8
!	nt = T/dt
	
!	NUMERO DE REYNOLDS
	Re = 2000._8
!	DENSIDAD
	rho = 1d3
!	VISCOSIDAD
	gama = 1.5_8
		
!	CRITERIO RUTINA TDMA_PHI
	criterio1_TDMA = 1d-3	! criterio1 -> max(RESTO)
	criterio2_TDMA = 100	! criterio2 -> iteraciones
	
!	CRITERIO CONVERGENCIA
	criterio1_convergencia = 1d-3 ! crit_conver -> max(eq_cdm_vol)
	criterio2_convergencia = 20
	
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

!	VELOCIDAD ENTRADA
	u_init = Re*gama/(rho*Ly)

!	PASO DE TIEMPO
	CFL = 1._8
	dt = CFL * dx / u_init
	T = 10._8
	nt = T/dt


!	MALLA PRINCIPAL ( x_1 , y_1 )
	allocate(x_1(ny+2),y_1(nx+2))
	x_1 = (/ ( (dfloat(i)-1.0_8)*dx,i=0,nx+1 ) /)
	y_1 = (/ ( (dfloat(i)-1.0_8)*dy,i=0,ny+1 ) /) 



!	MALLA SECUNDARIA ( x_2 , y_2 )
	allocate(x_2(ny+2),y_2(nx+2))
	x_2 = (/ ( (dfloat(i)-0.5_8)*dx,i=0,nx+1 ) /)
	y_2 = (/ ( (dfloat(i)-0.5_8)*dy,i=0,ny+1 ) /)



!	VELOCIDAD INICIAL

!		en t_(n)
	allocate(u_0(ny+2,nx+2),v_0(ny+2,nx+2))
!	do i=1,ny+2	
!		u_0(i,:) = (/ (0.1_8*u_init*(nx+2-j)*(i-1)/((ny+2)*(nx+2)),j=1,nx+2) /)
!	end do
	u_0 = 0.1_8
	v_0 = 0._8
	
!		en t_(n+1)
	allocate(u_1(ny+2,nx+2),v_1(ny+2,nx+2))
!	do i=1,ny+2	
!		u_1(i,:) = (/ (0.1_8*u_init*(nx+2-j)*(i-1)/((ny+2)*(nx+2)),j=1,nx+2) /)
!	end do
	u_1 = 0.1_8
	v_1 = 0._8
	
!	CAMPO DE PRESION INICIAL
!		en t_(n)
	allocate(P(ny+2,nx+2))
	P = 0._8

!	FUNCION AUXILIAR PHI
	allocate(phi(ny+2,nx+2))
	phi = 0.0_8
		
!	IMPOSICION DE LAS CONDICIONES DE CONTORNO
	call CC(nx,ny,u_0,v_0,u_init)
	call CC(nx,ny,u_1,v_1,u_init)


!:::::::::::::::::::::::::::::::::::::::::::::::::::::

!			INICIO INTEGRACION TEMPORAL

	allocate(u_pred(ny+2,nx+2),v_pred(ny+2,nx+2))
	u_pred = 0._8
	v_pred = 0._8


	info_ITER = 0

	do k = 1,5

		niter = 0

!:::::::::::::::::::::::::::::::::::::::::::::::::::::

		do while( info_ITER .eq. 0 )
		
			phi = 0._8
			niter = niter + 1 
		
		!	PREDICCION DEL CAMPO DE VELOCIDAD NO SOLENOIDAL

			call PREDICCION_VELOCIDAD(nx,ny,dx,dy,dt,Re,P,u_1,u_0,v_1,v_0,u_pred,v_pred)
			call CC(nx,ny,u_pred,v_pred,u_init)
			
			open(unit=10,file=dummy//'v.dat',access='SEQUENTIAL')
			do i=1,ny+2
				write(10,*) u_pred(i,:)
			end do	
			close(10)
			
									
		!	RESOLUCION DE LA ECUACION DE POISSON PARA LA PRESION

			contador = 0
			info_TDMA = 0
			do while ( info_TDMA .eq. 0 ) 
				contador = contador + 1		
				call TDMA_PHI(nx,ny,dx,dy,dt,u_pred,v_pred,phi,R)
				if ( R .lt. criterio1_TDMA .or. contador .gt. criterio2_TDMA ) then
					info_TDMA = 1
				end if	
			end do
			info_TDMA = 0

		!	CORRECCION DEL CAMPO DE VELOCIDAD Y PRESION

			call CORRECCION_PRESION(nx,ny,dx,dy,Re,P,phi,u_pred,v_pred)

			open(unit=10,file=dummy//'phi.dat',access='SEQUENTIAL')
			do i=1,ny+2
				write(10,*) phi(i,:)
			end do	
			close(10)
			
			open(unit=10,file=dummy//'p.dat',access='SEQUENTIAL')
			do i=1,ny+2
				write(10,*) P(i,:)
			end do	
			close(10)
		
			print*,CDM
		
			read(*,*)	
								
		!	VERIFICAR CONVERGENCIA DE LAS VARIABLES

			call CALCULO_CONVERGENCIA(nx,ny,dx,dy,dt,u_pred,v_pred,phi,CDM)
									
			if ( CDM .lt. criterio1_convergencia .or. & 
				& niter .gt. criterio2_convergencia ) then
				info_ITER = 1
			end if

!			print*, CDM	

		end do
		
		return
			
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
	
		info_ITER=0
			
!			CORRECCION DE LA VELOCIDAD
			
		call CORRECCION_VELOCIDAD(nx,ny,dx,dy,dt,u_pred,v_pred,phi,u_0,u_1,v_0,v_1)	

		
	end do
	
	

!!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end program
