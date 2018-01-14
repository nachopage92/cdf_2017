program T2_CDF_2017_2S

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	use variables

	implicit none
	
	INTERFACE

		subroutine CC(u,v)
			use variables
			real(kind=8),dimension(ny+2,nx+2),intent(inout) :: u,v
		end subroutine
		
		subroutine PREDICCION_VELOCIDAD(P,u_1,u_0,v_1,v_0,u_pred,v_pred)
			use variables
			real(kind=8),dimension(ny+2,nx+2),intent(in):: &
				& P , u_1 , u_0 , v_1 , v_0
			real(kind=8),dimension(ny+2,nx+2),intent(out):: &
				& u_pred , v_pred
		end subroutine
		
		subroutine CORRECCION_PRESION(P,u,v,R,phi)
			use variables
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: u,v
			real(kind=8),dimension(ny+2,nx+2),intent(inout) :: P
			real(kind=8),intent(out) :: R
			real(kind=8),dimension(ny+2,nx+2),intent(out) :: phi
		end subroutine

		subroutine CORRECCION_VELOCIDAD(u_pred,v_pred,phi,u_0,u_1,v_0,v_1)
			use variables
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: phi
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: u_pred,v_pred
			real(kind=8),dimension(ny+2,nx+2),intent(inout) :: u_0,u_1,v_0,v_1
		end subroutine		
		
		subroutine EXPORTAR_DATOS(x,y,u,v)
			use variables
			real(kind=8),dimension(nx+2),intent(in) :: x
			real(kind=8),dimension(ny+2),intent(in) :: y
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: u,v	
		end subroutine

	END INTERFACE
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	integer :: i,j,k,contador,niter
	
	real(kind=8) :: R
	
	real(kind=8),dimension(nx+2) :: x_1,x_2
	
	real(kind=8),dimension(ny+2) :: y_1,y_2
	
	real(kind=8),dimension(ny+2,nx+2) :: &
		& u_0,v_0,u_1,v_1,P,phi,u_pred,v_pred	 
		
	character(len=4) :: title = 'cfd_' 
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	
!	PREAMBULO 

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
	

!	MALLA PRINCIPAL
	x_1 = (/ ( (dfloat(i)-1.0_8)*dx,i=1,nx+2 ) /)
	y_1 = (/ ( (dfloat(i)-1.0_8)*dy,i=1,ny+2 ) /) 
	
	
!	MALLA SECUNDARIA
	x_2 = (/ ( (dfloat(i)-0.5_8)*dx,i=1,nx+2 ) /)
	y_2 = (/ ( (dfloat(i)-0.5_8)*dy,i=1,ny+2 ) /)


!	VELOCIDAD INICIAL
!		en t_(n)
	u_0 = 0._8
	u_0(2:ny+1,nx+1) = u_init
	u_0(2:ny+1,nx+2) = u_init
	v_0 = 0._8
!		en t_(n+1)
	u_1 = 0._8
	u_1(2:ny+1,nx+1) = u_init
	u_1(2:ny+1,nx+2) = u_init
	v_1 = 0._8
	
	
!	CAMPO DE PRESION INICIAL
!		en t_(n)
	P = 0._8
	

!	IMPOSICION DE LAS CONDICIONES DE CONTORNO
	call CC(u_0,v_0)
	call CC(u_1,v_1)
	
		
!:::::::::::::::::::::::::::::::::::::::::::::::::::::


!			INICIO INTEGRACION TEMPORAL

	do k = 1 , nt
	
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	
		niter = 0
		
		do		

			u_pred = 0._8
			v_pred = 0._8
			niter = niter + 1
				
		!	PREDICCION DEL CAMPO DE VELOCIDAD NO SOLENOIDAL
			call PREDICCION_VELOCIDAD(P,u_1,u_0,v_1,v_0,u_pred,v_pred)	
												
!			CORRECCION DE LA PRESION
			call CORRECCION_PRESION(P,u_pred,v_pred,R,phi)
			if ( R .lt. criterio1 .or. &
			 &	niter .gt. criterio2 ) then
				EXIT
			end if
			
		end do
	
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::


!			CORRECCION DE LA VELOCIDAD
		
		call CORRECCION_VELOCIDAD(u_pred,v_pred,phi,u_0,u_1,v_0,v_1)

		call EXPORTAR_DATOS(x_1,y_1,u_1,v_1)
		
		read(*,*)
		
		call CC(u_1,v_1)
		
	end do
			

!!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end program
