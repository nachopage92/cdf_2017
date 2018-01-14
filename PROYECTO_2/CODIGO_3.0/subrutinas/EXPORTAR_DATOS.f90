subroutine EXPORTAR_DATOS(x,y,u,v)
	use variables
	real(kind=8),dimension(nx+2),intent(in) :: x
	real(kind=8),dimension(ny+2),intent(in) :: y
	real(kind=8),dimension(ny+2,nx+2),intent(in) :: u,v	
	real(kind=8),dimension(ny,nx,2) :: velocity_field
	integer :: i,j
	character(len=4) :: dummy='cfd_'
			 
	!CAMPO DE VELOCIDAD
	do i=2,ny+1
		do j=2,nx+1
			velocity_field(i-1,j-1,1) = (u(i,j+1)+u(i,j))*0.5_8
			velocity_field(i-1,j-1,2) = (v(i-1,j)+v(i,j))*0.5_8			 
		end do
	end do

	!EXPORTAR DATOS PARA GRAFICAR
	open(unit=10,file=dummy//'velocity_field.dat',access='SEQUENTIAL')
		do i=1,ny
			do j=1,nx
				write(10,*)  x(j) , y(i) , &
				& velocity_field(i,j,1)/50._8 , velocity_field(i,j,2)/50._8
			end do
		end do	
	 close(10)	
	 
	 call system( 'gnuplot plot' )

end subroutine
