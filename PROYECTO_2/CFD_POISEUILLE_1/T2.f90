program T2_CDF_2017_2S

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	use variables

	implicit none
	
	INTERFACE

		subroutine CONDICIONES_INICIALES(u,v)
			use variables
			implicit none	
			real(kind=8),dimension(ny+2,nx+2),intent(inout)::u,v
		end subroutine

		subroutine CONDICIONES_CONTORNO(u,v)
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
		
		subroutine EXPORTAR_DATOS(x,y,u,v,P,k)
			use variables
			integer,intent(in)::k
			real(kind=8),dimension(nx+2),intent(in) :: x
			real(kind=8),dimension(ny+2),intent(in) :: y
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: u,v,P
		end subroutine
		
		subroutine EXPORTAR_DATOS_2(x,y,u,v,P)
			use variables
			real(kind=8),dimension(nx+2),intent(in) :: x
			real(kind=8),dimension(ny+2),intent(in) :: y
			real(kind=8),dimension(ny+2,nx+2),intent(in) :: u,v,P
		end subroutine

	END INTERFACE
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	integer :: i,j,k,contador,niter
	
	real(kind=8) :: R,start,finish
	
	real(kind=8),dimension(nx+2) :: x
	
	real(kind=8),dimension(ny+2) :: y
	
	real(kind=8),dimension(ny+2,nx+2) :: &
		& u_0,v_0,u_1,v_1,P,phi,u_pred,v_pred	 
		
	character(len=4) :: title = 'cfd_' 
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	call cpu_time(start)
	
	PRINT*, 'PROYECTO 2 - DINAMICA DE FLUIDOS COMPUTACIONAL'
	PRINT*, 'IGNACIO APABLAZA BENAVIDES'
	PRINT*, 'ROL. 201141007-6 ; ignacio.apablaza@alumnos.usm.cl'
	PRINT*, ''
	PRINT*, ''
	PRINT*, ''
	
	
!	PREAMBULO 

!	MALLA PRINCIPAL
	x = (/ ( (dfloat(i)-1.5_8)*dx,i=1,nx+2 ) /)
	y = (/ ( (dfloat(i)-1.5_8)*dy,i=1,ny+2 ) /) 


!	CAMPO DE VELOCIDAD INICIAL
	call CONDICIONES_INICIALES(u_0,v_0)
	call CONDICIONES_INICIALES(u_1,v_1)


!	CAMPO DE PRESION INICIAL
!		en t_(n)
	P = 0._8

	
!	IMPOSICION DE LAS CONDICIONES DE CONTORNO
	call CONDICIONES_CONTORNO(u_0,v_0)
	call CONDICIONES_CONTORNO(u_1,v_1)
	
		
!:::::::::::::::::::::::::::::::::::::::::::::::::::::


!			INICIO INTEGRACION TEMPORAL

	do k = 1 , nt
	
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	
		niter = 0
		
!			ITERACIÓN DE CORRECCIÓN DE PRESION		
		do		
		
			u_pred = 0._8
			v_pred = 0._8
			niter = niter + 1
				
!			PREDICCION DEL CAMPO DE VELOCIDAD NO SOLENOIDAL
			call PREDICCION_VELOCIDAD(P,u_1,u_0,v_1,v_0,u_pred,v_pred)	
												
!			CORRECCION DE LA PRESION
			call CORRECCION_PRESION(P,u_pred,v_pred,R,phi)
			
!				criterio1 : limite maximo convergencia variables
!				criterio2 : limite maximo de iteraciones
			if ( R .lt. criterio1 .or. &
			 &	niter .gt. criterio2 ) then
				print*, R,niter,k,nt
				call EXPORTAR_DATOS(x,y,u_1,v_1,P,k)
				EXIT
			end if
			
		end do
	
	
!:::::::::::::::::::::::::::::::::::::::::::::::::::::


!			CORRECCION DE LA VELOCIDAD
		call CORRECCION_VELOCIDAD(u_pred,v_pred,phi,u_0,u_1,v_0,v_1)
		
		
!			exporta los datos y los grafica una vez
!			alcanzada la estabilidad en el tiempo
		if ( (maxval(u_1-u_0) .lt. estabilidad) .and. &
		& maxval(v_1-v_0) .lt. estabilidad ) then
			exit
		end if
		
!			aplica las condiciones de contorno 
!			para una nueva iteracion
		call CONDICIONES_CONTORNO(u_1,v_1)
		
	end do
	
!	REQUIERE TENER INSTALADO imagemagick
!	sudo apt-get install imagemagick
	call system('cd ./animacion1/ && convert -delay 1'//&
	&' -loop 0 *.png animacion1.gif && rm *.png && rm *.dat')
	call system('cd ./animacion2/ && convert -delay 1'//&
	&' -loop 0 *.png animacion2.gif && rm *.png && rm *.dat')
	call system('cd ./animacion3/ && convert -delay 1'//&
	&' -loop 0 *.png animacion3.gif && rm *.png && rm *.dat')

	call EXPORTAR_DATOS_2(x,y,u_1,v_1,P)
	
	call cpu_time(finish)
	
	PRINT*, '----------------------------------'
	PRINT*, 'Tiempo de calculo' , finish - start , 'segundos'

!!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end program

