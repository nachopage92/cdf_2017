module variables

!::::::::::::::::::::::::::::::::::::::::::::::

	implicit none

!	PARAMETROS

!	NUMERO DE VOLUMENES
!	(no considera nodos ficticios)
	integer,parameter :: nx = 50
	integer,parameter :: ny = 10
	
!	DIMENSION DOMINIO
	real(kind=8),parameter :: Lx = 5._8
	real(kind=8),parameter :: Ly = 1._8
	 
!	NUMERO DE REYNOLDS
	real(kind=8),parameter :: Re = 1000._8
	
!	DENSIDAD
	real(kind=8),parameter :: rho = 1000._8
	
!	VISCOSIDAD
	real(kind=8),parameter :: gama = 1.5_8

!	CRITERIO CONVERGENCIA
	real(kind=8),parameter :: criterio1 = 1d-4
	integer,parameter :: criterio2 = 10000
	
!	ESPACIO INTERNODAL
	real(kind=8),parameter :: dx = Lx/dfloat(nx)
	real(kind=8),parameter :: dy = Ly/dfloat(ny)	

!	VELOCIDAD ENTRADA
	real(kind=8),parameter :: &
		& u_init = 2._8*Re*gama/(rho*Ly)

!	NUMERO DE Courant-Friedrich-Levy
	real(kind=8),parameter :: CFL = 0.15_8
	
!	PASO DE TIEMPO
	real(kind=8),parameter :: dt = CFL * dx / u_init
	
!	INTERVALO DE INTEGRACION TEMPORAL
	real(kind=8),parameter :: T = 10._8
	
!	NUMERO DE PASOS
	integer,parameter :: nt = T/dt
	
!	CRITERIO ESTABILIDAD
	real(kind=8),parameter :: estabilidad = 1d-3
	
! NUMERO PI
	real(kind=8),parameter :: pi = 3.14159265359_8
	
!::::::::::::::::::::::::::::::::::::::::::::::

end module
