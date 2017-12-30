subroutine WPENS_U(u,v,P,i,j,F_w,F_e,F_n,F_s,u_W,u_P,u_E,u_N,u_S)
	implicit none
	integer,intent(in)::i,j
	real(kind=8),intent(in),dimension(:,:) :: u,v,P
	real(kind=8),intent(out):: &
		&F_w,F_e,F_n,F_s,u_W,u_E,u_N,u_S,u_P
		
	F_w = 0.5_8 * ( u(i,j) + u(i,j+1) )
	F_e = 0.5_8 * ( u(i,j-1) + u(i,j) )
	F_n = 0.5_8 * ( v(i,j) + v(i,j-1) ) 
	F_s = 0.5_8 * ( v(i-1,j) + v(i-1,j-1) )
	
	u_W = u(i,j-1)
	u_E = u(i,j+1)
	u_N = u(i+1,j)
	u_S = u(i-1,j)
	u_P = u(i,j)

	
end subroutine
