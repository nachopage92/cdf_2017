subroutine WPENS_V(u,v,P,i,j,F_w,F_e,F_n,F_s,v_W,v_P,v_E,v_N,v_S)
	implicit none
	integer,intent(in)::i,j
	real(kind=8),intent(in),dimension(:,:) :: u,v,P
	real(kind=8),intent(out):: &
		&F_w,F_e,F_n,F_s,v_W,v_E,v_N,v_S,v_P
		
	F_w = 0.5_8 * ( u(i+1,j) + u(i,j) )
	F_e = 0.5_8 * ( u(i+1,j+1) + u(i,j+1) )
	F_n = 0.5_8 * ( v(i+1,j) + v(i,j) ) 
	F_s = 0.5_8 * ( v(i,j) + v(i-1,j) )
	
	v_W = v(i,j-1)
	v_E = v(i,j+1)
	v_N = v(i+1,j)
	v_S = v(i-1,j)
	v_P = v(i,j)
	
end subroutine
