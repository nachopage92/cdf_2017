subroutine EXPORTAR_DATOS(x,y,u,v,P)
	use variables
	real(kind=8),dimension(nx+2),intent(in) :: x
	real(kind=8),dimension(ny+2),intent(in) :: y
	real(kind=8),dimension(ny+2,nx+2),intent(in) :: u,v,P
	real(kind=8),dimension(ny,nx,2) :: velocity_field
	integer :: i,j
	character(len=4) :: dummy='cfd_'
			 
!	CAMPO DE VELOCIDAD
	do i=2,ny+1
		do j=2,nx+1
			velocity_field(i-1,j-1,1) = (u(i,j+1)+u(i,j))*0.5_8
			velocity_field(i-1,j-1,2) = (v(i-1,j)+v(i,j))*0.5_8			 
		end do
	end do

!	EXPORTAR DATOS DE LA VELOCIDAD PARA GRAFICAR
	open(unit=10,file=dummy//'velocity_field.dat',access='SEQUENTIAL')
		do i=1,ny
			do j=1,nx
				write(10,*)  x(j) , y(i) , &
				& velocity_field(i,j,1)/50._8 , velocity_field(i,j,2)/50._8
			end do
		end do	
	 close(10)	
	 
!	EXPORTAR DATOS DE LA PRESION PARA GRAFICAR
	open(unit=10,file=dummy//'pressure_field.dat',access='SEQUENTIAL')
		do i=2,ny+1
			write(10,*) P(i,2:nx+1)
		end do	
	 close(10)	
	 
!	EXPORTAR DATOS DEL PERFIL
	open(unit=10,file=dummy//'velocity_profile.dat',access='SEQUENTIAL')
		do i=2,ny+1
			write(10,*) y(i) , u(i,(nx+2)/2) , POISEUILLE(P,Ly-y(i))
		end do	
	 close(10)	
	 
!	 GRAFICAR DATOS
	 call system( 'gnuplot plot' )	
	 
	contains
	
	real(kind=8) function POISEUILLE(P,y)
		use variables
		real(kind=8),dimension(ny+2,nx+2) :: P
		real(kind=8) :: y,Pmax,Pmin,maxx
		Pmax = sum( P( 2:ny+1 , 2 ) ) / ny
		Pmin = sum( P( 2:ny+1 , nx+1 ) ) / ny
		POISEUILLE = (Pmax-Pmin)/(4._8*Lx*gama) * (Ly**2._8 - y**2._8) - maxx
	end function
	

end subroutine
