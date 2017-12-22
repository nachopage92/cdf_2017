subroutine WPENS(phi,i,j,phi_W,phi_P,phi_E,phi_N,phi_S)
	integer,intent(in)::i,j
	real(kind=8),intent(in),dimension(:,:) :: phi
	real(kind=8),intent(out)::phi_W,phi_P,phi_E,phi_N,phi_S
	phi_W = phi(i,j-1)
	phi_P = phi(i,j)
	phi_E = phi(i,j+1)
	phi_S = phi(i-1,j)
	phi_N = phi(i+1,j)
end subroutine
