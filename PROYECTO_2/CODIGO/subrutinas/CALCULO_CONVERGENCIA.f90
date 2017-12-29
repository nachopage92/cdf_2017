subroutine CALCULO_CONVERGENCIA(nx,ny,dx,dy,dt,u,v,phi,CDM)

	implicit none
	
	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy,dt
	real(kind=8),dimension((nx+2),(ny+2)),intent(in) :: u,v,phi
	real(kind=8),intent(out):: CDM
	
	integer :: i,j
	real(kind=8),dimension(nx,ny):: CDM_vol
	
	do i=2,nx+1
		do j=2,ny+1
			CDM_vol(i-1,j-1) = & 
				& (u(i+1,j)-u(i,j))*dy + (v(i,j)-v(i,j-1))*dx &
				& - dt*(phi(i,j+1)-2._8*phi(i,j)+phi(i,j-1))*(dy/dx) &
				& - dt*(phi(i+1,j)-2._8*phi(i,j)+phi(i-1,j))*(dx/dy)
		end do
	end do
	
	CDM = sum( abs(CDM_vol) )

end subroutine
