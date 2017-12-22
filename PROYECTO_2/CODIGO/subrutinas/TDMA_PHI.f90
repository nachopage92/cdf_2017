subroutine TDMA_PHI(nx,ny,dx,dy,dt,u_0,u_1,v_0,v_1,phi,R)

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
	real(kind=8),dimension((nx+2),(ny+2)),intent(in) :: &
		& u_0,u_1,v_0,v_1
	real(kind=8),intent(out) :: R
	real(kind=8),dimension((nx+2),(ny+2)),intent(inout) :: & 
		& phi
	
	
	integer :: i,j,contador
	real(kind=8) :: B
	real(kind=8),dimension((nx+2),(ny+2)) :: phi_
	real(kind=8),dimension(nx*ny) :: R_vec
	real(kind=8),dimension(nx) :: RHS_x,diagx_1,diagx_2,diagx_3,phi_x
	real(kind=8),dimension(ny) :: RHS_y,diagy_1,diagy_2,diagy_3,phi_y
	


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	phi_ = phi

!		PREDICTOR PHI1_E, PHI1_P, PHI1_W (DIRECCION X)

!			recorriendo fila por fila
	do j=2,ny+1
	
		contador = 0
		
!			se crea el vector RHS_x
		do i=2,nx+1
			
			contador = contador + 1
			 
			B = 3._8/(2._8*dt)*(  &
			&u_1(i,j+1)-u_1(i,j)+ &
			&v_1(i,j)-v_1(i-1,j) )

			RHS_x(contador) = B - &
			& (dx/dy)*phi_(i+1,j) - &
			& (dx/dy)*phi_(i-1,j)
	 
		end do 
		
		diagx_1 = dy/dx
		diagx_2 = -2._8*( dy/dx + dx/dy )
		diagx_3 = dy/dx	
		
!			resuelve phi en cada fila
		call thomas(nx,diagx_1,diagx_2,diagx_3,phi_x,RHS_x)
		
!			se reemplaza la solucion de phi (tmp_x) en su
!			fila correspondiente a la matriz phi
		phi_(:,j) = phi_x(:)
		
	end do

!-----------------------------------------------------------		

!		PREDICTOR PHI1_N, PHI1_P, PHI1_S (DIRECCION Y)

!			recorriendo columna por columna
	do i=2,nx+1
	
		contador = 0
		
!			se crea el vector RHS_y
		do j=2,ny+1
			
			contador = contador + 1
			 
			B = 3._8/(2._8*dt)*(  &
			&u_1(i,j+1)-u_1(i,j)+ &
			&v_1(i,j)-v_1(i-1,j) )

			RHS_y(contador) = B - &
			& (dx/dy)*phi_(i+1,j) - &
			& (dx/dy)*phi_(i-1,j)
	 
		end do 
		
		diagy_1 = dy/dx
		diagy_2 = -2._8*( dy/dx + dx/dy )
		diagy_3 = dy/dx	
		
!			resuelve phi en cada columna
		call thomas(ny,diagy_1,diagy_2,diagy_3,phi_y,RHS_y)
		
!			se reemplaza la solucion de phi (tmp_x) en su
!			fila correspondiente a la matriz phi
		phi_(i,:) = phi_y(:)
		
	end do

	!CALCULO DE RESTO R (CONVERGENCIA)
	contador = 0
	
	do i=2,nx+1
	
		do j=2,ny+1
	
			contador = contador + 1
			
			B = 3._8/(2._8*dt)*(  &
			&u_1(i,j+1)-u_1(i,j)+ &
			&v_1(i,j)-v_1(i-1,j) )
			
			R_vec(contador) = &
			& dx/dy*(phi(i+1,j)+phi(i-1,j)) &
			& -2._8*(dx/dy+dy/dx)*phi(i,j) &  
			& + dy/dx*(phi(i,j+1)+phi(i,j-1)) &
			& - B
			
		end do
		
	end do
	
	R = maxval(R_vec)

end subroutine
