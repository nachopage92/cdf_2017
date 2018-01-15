subroutine EXPORTAR_DATOS(x,y,u,v,P,k)
	use variables
	integer,intent(in)::k
	real(kind=8),dimension(nx+2),intent(in) :: x
	real(kind=8),dimension(ny+2),intent(in) :: y
	real(kind=8),dimension(ny+2,nx+2),intent(in) :: u,v,P
	real(kind=8),dimension(ny,nx,2) :: velocity_field
	integer :: i,j
	character(len=1) :: dummy='.'
	character(len=6) :: n='100000'
	
	write(n,'(I6)') 100000+k
	
!	CAMPO DE VELOCIDAD
	do i=2,ny+1
		do j=2,nx+1
			velocity_field(i-1,j-1,1) = (u(i,j+1)+u(i,j))*0.5_8
			velocity_field(i-1,j-1,2) = (v(i-1,j)+v(i,j))*0.5_8			 
		end do
	end do

!	EXPORTAR DATOS DE LA VELOCIDAD PARA GRAFICAR
	open(unit=10,file=dummy//'/animacion1/velocity_field_'//n//'.dat',access='SEQUENTIAL')
		do i=1,ny
			do j=1,nx
				write(10,*)  x(j) , y(i) , &
				& velocity_field(i,j,1)/50._8 , velocity_field(i,j,2)/50._8
			end do
		end do	
	 close(10)	
	 
!	EXPORTAR DATOS DE LA PRESION PARA GRAFICAR
	open(unit=10,file=dummy//'/animacion2/pressure_field_'//n//'.dat',access='SEQUENTIAL')
		do i=2,ny+1
			write(10,*) P(i,2:nx+1)
		end do	
	 close(10)	
	 
!	EXPORTAR DATOS DEL PERFIL
	open(unit=10,file=dummy//'/animacion3/velocity_profile_'//n//'.dat',access='SEQUENTIAL')
		do i=2,ny+1
			write(10,*) y(i) , u(i,(nx+2)/2) , POISEUILLE(P,Ly-y(i))
		end do	
	 close(10)	
	 
!	 GRAFICAR DATOS
	call GNUPLOT_SCRIPT()
	call system( 'gnuplot plot' )	
	 
	contains
	
	real(kind=8) function POISEUILLE(P,y)
		use variables
		real(kind=8),dimension(ny+2,nx+2) :: P
		real(kind=8) :: y,Pmax,Pmin
		Pmax = sum( P( 2:ny+1 , 2 ) ) / ny
		Pmin = sum( P( 2:ny+1 , nx+1 ) ) / ny
		POISEUILLE = (Pmax-Pmin)/(4._8*Lx*gama) * (Ly**2._8 - y**2._8)
	end function
	
	subroutine GNUPLOT_SCRIPT()
	
		open(unit=10,file=dummy//'/plot',access='SEQUENTIAL') 
		write(10,*) "set terminal png #size 600,400"
		write(10,*)
		write(10,*) "set pm3d map"
		write(10,*) "#set view map"
		write(10,*) "set pm3d interpolate 5,5"
		write(10,*)
		write(10,*) "set output './animacion3/velocity_profile_"//n//".png'"
		write(10,*) "set key left top"
		write(10,*) "set title 'Perfil de velocidad en x = Lx/2'"
		write(10,*) "set xlabel 'y'"
		write(10,*) "set ylabel 'u'"
		write(10,*) "plot './animacion3/velocity_profile_"//n//".dat' using 1:2 with points title 'Result. Vol. Finitos' , \"
		write(10,*) " '' using 1:3 with lines title 'Result. Analítico'"
		write(10,*)
		write(10,*) "set output './animacion2/pressure_field_"//n//".png'"
		write(10,*) "set title 'Campo de presión P en cada volumen'"
		write(10,*) "set xlabel 'j-ésimo volumen'"
		write(10,*) "set ylabel 'i-ésimo volumen'"
		write(10,*) "splot './animacion2/pressure_field_"//n//".dat' matrix notitle "
		write(10,*)
		write(10,*) "set output './animacion1/velocity_field_"//n//".png'"
		write(10,*) "set title 'Campo de velocidad V'"
!		write(10,*) "set yrange [0:"//Ly//"]"
!		write(10,*) "set xrange [0:"//Lx//"]"
		write(10,*) "set xlabel 'x'"
		write(10,*) "set ylabel 'y'"
		write(10,*) "plot './animacion1/velocity_field_"//n//".dat' using 1:2:3:4 with vectors head size 0.1,20,60 filled notitle"
		close(10)
	
	end subroutine

end subroutine
