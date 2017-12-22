subroutine thomas(n,d,a,c,x,b)

!::::::::::::::::::::::::::::::::::::::::::::::::::

	implicit none
	
	!input
	integer,intent(in)::n
	double precision,dimension(n),intent(in) :: d,a,c,b
	
	!output
	double precision,dimension(n),intent(out)	:: x
	
	!variables locales
	integer							:: i
	double precision,dimension(n)	:: y,alfa
	double precision,dimension(n)	:: beta

!::::::::::::::::::::::::::::::::::::::::::::::::::
	
	!calculo de los coeficientes de la matriz L y U 
	!(alfa, beta y gama) 
	
	!alfa 	= (alfa_1 , alfa_2 , ... , alfa_(n-1) , alfa_n)
	!beta 	= (0	  , beta_2 , ... , beta_(n-1) , beta_n)
	!c 		= (c_1 	  , c_2    , ... , c(n-1) 	  , 0)
	
	alfa(1) = a(1)
	do i=2,n
		beta(i)=d(i-1)/alfa(i-1)
		alfa(i)=a(i)-beta(i)*c(i-1)
	end do

	!resolver sistema Ly=b
	y(1)=b(1)
	do i=2,n
		y(i) = b(i) - beta(i) * y(i-1)
	end do
	
	!resolver sistema Ux=y
	x(n)=y(n)/alfa(n)
	do i=n-1,1,-1
		x(i) = ( y(i)-c(i)*x(i+1) ) / alfa(i)
	end do
	
end subroutine
