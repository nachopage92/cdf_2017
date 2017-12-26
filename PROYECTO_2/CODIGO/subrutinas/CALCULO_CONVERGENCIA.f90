subroutine CALCULO_CONVERGENCIA(nx,ny,dx,dy,u,v,CDM_vol)

	integer,intent(in) :: nx,ny
	real(kind=8),intent(in) :: dx,dy
	real(kind=8),dimension((nx+2),(ny+2)),intent(in) :: u,v
	real(kind=8),dimension(nx,ny),intent(out):: CDM_vol
	
	integer :: i,j
	
	do i=2,nx+1
		do j=2,ny+1
			CDM_vol(i-1,j-1) = ( u(i+1,j) - u(i,j) )*dy + ( v(i,j) - v(i,j-1) )*dx 
		end do
	end do

end subroutine
