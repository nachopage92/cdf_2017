program T2_CDF_2017_2S

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	use variables

	implicit none
	
	INTERFACE

		subroutine CC(u,v)
			use variables
			implicit none
			real(kind=8),dimension(ny+2,nx+2),intent(inout) :: u,v
		end subroutine
		
		subroutine PREDICCION_VELOCIDAD(P,u_1,u_0,v_1,v_0,u_pred,v_pred)
			use variables
			implicit none
			real(kind=8),dimension(ny+2,nx+2),intent(in):: &
				& P , u_1 , u_0 , v_1 , v_0
			real(kind=8),dimension(ny+2,nx+2),intent(out):: &
				& u_pred , v_pred
		end subroutine
		
		subroutine CORRECCION_PRESION(P,u,v,R,phi)
			use variables
			implicit none
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: u,v
			real(kind=8),dimension(ny+2,nx+2),intent(inout) :: P
			real(kind=8),intent(out) :: R
			real(kind=8),dimension(ny+2,nx+2),intent(out) :: phi
		end subroutine

		subroutine CORRECCION_VELOCIDAD(u_pred,v_pred,phi,u_0,u_1,v_0,v_1)
			use variables
			implicit none
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: phi
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: u_pred,v_pred
			real(kind=8),dimension(ny+2,nx+2),intent(inout) :: u_0,u_1,v_0,v_1
		end subroutine		

	END INTERFACE
	
	integer :: i,j,k,contador,niter
	
	real(kind=8) :: R
	
	real(kind=8),allocatable :: &
		& x_1(:),x_2(:),y_1(:),y_2(:), &
		& u_0(:,:),v_0(:,:), &
		& u_1(:,:),v_1(:,:), &
		P(:,:), phi(:,:), &
		& u_pred(:,:),v_pred(:,:)
		
	character(len=3) :: dummy = 'cfd' 
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
	
!	PREAMBULO 
!	predefinido en el modulo 'variables'
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::	


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
	u_0 = 0._8
	u_0(:,nx+1) = u_init
	u_0(:,nx+2) = u_init
	v_0 = 0._8
!		en t_(n+1)
	allocate(u_1(ny+2,nx+2),v_1(ny+2,nx+2))
	u_1 = 0._8
	u_1(:,nx+1) = u_init
	u_1(:,nx+2) = u_init
	v_1 = 0._8
	
	
!	CAMPO DE PRESION INICIAL
!		en t_(n)
	allocate(P(ny+2,nx+2))
	P = 0._8

	
!	FUNCION AUXILIAR PHI
	allocate(phi(ny+2,nx+2))


!	IMPOSICION DE LAS CONDICIONES DE CONTORNO
	call CC(u_0,v_0)
	call CC(u_1,v_1)
	
		
!:::::::::::::::::::::::::::::::::::::::::::::::::::::


!			INICIO INTEGRACION TEMPORAL

!	PREDICTOR DE VELOCIDAD
	allocate(u_pred(ny+2,nx+2),v_pred(ny+2,nx+2))

	do k = 1 , nt
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
	
		niter = 0
		
		do
			u_pred = 0._8
			v_pred = 0._8
			niter = niter + 1
				
				
		!	PREDICCION DEL CAMPO DE VELOCIDAD NO SOLENOIDAL
			call PREDICCION_VELOCIDAD(P,u_1,u_0,v_1,v_0,u_pred,v_pred)		
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if ( niter .eq. 400 ) then	
!	write(1,*) 'u_pred'		
!	do i=1,ny+2
!		write(1,*) u_pred(i,:)
!	end do
!	write(2,*) 'phi'	
!	do i=1,ny+2
!			write(2,*) phi(i,:)
!	end do
!	write(3,*) 'P'	
!	do i=1,ny+2
!			write(3,*) P(i,:)
!	end do
!	print*,niter
!	return
!end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
												
!			CORRECCION DE LA PRESION
			call CORRECCION_PRESION(P,u_pred,v_pred,R,phi)
			if ( R .lt. criterio1 .or. &
			 &	niter .gt. criterio2 ) then
				EXIT
			end if

		end do
		
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!do i=1,ny+2
!	write(1,*) u_pred(i,:)
!end do
!do i=1,ny+2
!	write(2,*) P(i,:)
!end do
!do i=1,ny+2
!	write(3,*) phi(i,:)
!end do
!print*,R,niter
!return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		
!:::::::::::::::::::::::::::::::::::::::::::::::::::::

!			CORRECCION DE LA VELOCIDAD
		
		call CORRECCION_VELOCIDAD(u_pred,v_pred,phi,u_0,u_1,v_0,v_1)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(unit=10,file=dummy//'u.dat',access='SEQUENTIAL')
	do i=2,ny+1
		write(10,*) u_1(i,:)
	end do	
 close(10)	
open(unit=10,file=dummy//'u0.dat',access='SEQUENTIAL')
	do i=2,ny+1
		write(10,*) u_0(i,:)
	end do	
 close(10)
open(unit=10,file=dummy//'phi.dat',access='SEQUENTIAL')
	do i=2,ny+1
		write(10,*) phi(i,1:nx+1)
	end do	
 close(10)
open(unit=10,file=dummy//'p.dat',access='SEQUENTIAL')
	do i=2,ny+1
		write(10,*) P(i,1:nx+1)
	end do	
 close(10)
call system( 'gnuplot plot' )
print*, 'resto',R
print*, 'n iter', niter
read(*,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
	call CC(u_1,v_1)
		
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
		
	end do
	
	

!!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end program
