subroutine TDMA_PHI(nx,ny,dx,dy,dt,u,v,phi,R)

	implicit none
	
	INTERFACE
	
		subroutine thomas(n,d,a,c,x,b)
			implicit none
			integer,intent(in)::n
			double precision,dimension(n),intent(in) :: d,a,c,b
			double precision,dimension(n),intent(out)	:: x
		end subroutine
		
	END INTERFACE
	
		
	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy,dt
	real(kind=8),dimension(ny+2,nx+2),intent(in) :: &
		& u,v
	real(kind=8),intent(out) :: R
	real(kind=8),dimension(ny+2,nx+2),intent(inout) :: & 
		& phi
	
	
	integer :: i,j,k,contador
	real(kind=8) :: B
	real(kind=8),dimension(ny+2,nx+2) :: phi_
	real(kind=8),dimension((nx-2)*(ny-2)) :: R_vec
	real(kind=8),dimension(nx-2) :: RHS_x,diagx_1,diagx_2,diagx_3,phi_x
	real(kind=8),dimension(ny-2) :: RHS_y,diagy_1,diagy_2,diagy_3,phi_y
	


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	phi_ = phi

!		PREDICTOR PHI1_E, PHI1_P, PHI1_W (DIRECCION X)

!			recorriendo fila por fila
	do i=3,ny
	
		contador = 0
		
!			se crea el vector RHS_x
		do j=3,nx
			
			contador = contador + 1
			 
			B = 3._8/(2._8*dt)*(&
			& ( u(i,j+1)-u(i,j) )*dy + &
			& ( v(i,j)-v(i-1,j) )*dx )

			RHS_x(contador) = - B + &
			& (dx/dy)*phi_(i+1,j) + &
			& (dx/dy)*phi_(i-1,j)
	 
		end do 
		
		diagx_1 = -dy/dx
		diagx_2 = 2._8*( dy/dx + dx/dy )
		diagx_3 = -dy/dx	
		
!			resuelve phi en cada fila

		call thomas(nx-2,diagx_1,diagx_2,diagx_3,phi_x,RHS_x)

!			se reemplaza la solucion de phi (tmp_x) en su
!			fila correspondiente a la matriz phi
		phi_(i,3:nx) = phi_x(:)
				
	end do

!-----------------------------------------------------------		

!		PREDICTOR PHI1_N, PHI1_P, PHI1_S (DIRECCION Y)

!			recorriendo columna por columna
	do j=3,nx
	
		contador = 0
		
!			se crea el vector RHS_y
		do i=3,ny
			
			contador = contador + 1
			 
			B = 3._8/(2._8*dt)*(  &
			& ( u(i,j+1)-u(i,j) )*dy + &
			& ( v(i,j)-v(i-1,j) )*dx )

			RHS_y(contador) = - B + &
			& (dy/dx)*phi_(i,j+1) + &
			& (dy/dx)*phi_(i,j-1)
	 
		end do 
		
		diagy_1 = -dx/dy
		diagy_2 = 2._8*( dy/dx + dx/dy )
		diagy_3 = -dx/dy	
		
!			resuelve phi en cada columna
		call thomas(ny-2,diagy_1,diagy_2,diagy_3,phi_y,RHS_y)
		
!			se reemplaza la solucion de phi (tmp_x) en su
!			fila correspondiente a la matriz phi
		phi_(3:ny,j) = phi_y(:)
		
	end do
	
	phi = phi_

!-----------------------------------------------------------		

	!CALCULO DE RESTO R (CONVERGENCIA)
	
	contador = 0
	
	do i=3,ny
	
		do j=3,nx
	
			contador = contador + 1
			
			B = 3._8/(2._8*dt)*(  &
			& ( u(i,j+1)-u(i,j) )*dy+ &
			& ( v(i,j)-v(i-1,j) )*dx )
			
			R_vec(contador) = &
			& -dx/dy*(phi(i+1,j)+phi(i-1,j)) &
			& + 2._8*(dx/dy+dy/dx)*phi(i,j) &  
			& - dy/dx*(phi(i,j+1)+phi(i,j-1)) &
			& + B
			
		end do
		
	end do
			
	R = maxval(abs(R_vec))

end subroutine
