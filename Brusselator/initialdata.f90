SUBROUTINE initialdata(Nx,Ny,Nz,x,y,z,u,v,decomp)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine gets initial data for Brusselator equation
	! in 3 dimensions
	! u_t = a + u^2v - (b+1)u + Du u_xx
    ! v_t = bu - u^2v + Dv v_xx
	!
	! The boundary conditions are u(x=-Lx*\pi,y,z)=u(x=Lx*\pi,y,z), 
	!	u(x,y=-Ly*\pi,z)=u(x,y=Ly*\pi,z),u(x,y,z=-Ly*\pi)=u(x,y,z=Ly*\pi)
	! The initial condition is 
	! u=0.5+exp(-1.0-(x**2+y**2+z**2))
    ! v=0.1+exp(-1.0-(x**2+y**2+z**2))
	!
	! INPUT
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  Nz				= number of modes in z - power of 2 for FFT
	! .. Vectors ..
	!  x				= x locations
	!  y				= y locations
	!  z				= z locations
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!					see http://www.2decomp.org/ for more information
	! OUTPUT
	!
	! .. Arrays ..
	!  u 				= initial solution
	!  uold 			= approximate solution based on time derivative of
	!					initial solution
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
	!
	! REFERENCES
	!
	! ACKNOWLEDGEMENTS
	!
	! ACCURACY
	!		
	! ERROR INDICATORS AND WARNINGS
	!
	! FURTHER COMMENTS
	! Check that the initial iterate is consistent with the 
	! boundary conditions for the domain specified
	!--------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	! 2DECOMP&FFT	 -- Domain decomposition and Fast Fourier Library
	!			(http://www.2decomp.org/index.html)
	! MPI library
	USE decomp_2d
	IMPLICIT NONE					 
	INCLUDE 'mpif.h'
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)							:: Nx,Ny,Nz
	TYPE(DECOMP_INFO), INTENT(IN)						::  decomp
	REAL(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1)), INTENT(IN) :: x
	REAL(KIND=8), DIMENSION(decomp%xst(2):decomp%xen(2)), INTENT(IN) :: y
	REAL(KIND=8), DIMENSION(decomp%xst(3):decomp%xen(3)), INTENT(IN) :: z
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   							decomp%xst(2):decomp%xen(2),&
   							decomp%xst(3):decomp%xen(3)),&
   											INTENT(OUT) :: u,v
	INTEGER(kind=4)										:: i,j,k
	
    DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	 u(i,j,k)=15.0d0+8.0d0*sin(2.0*x(i))*cos(2.0*y(j)) !*exp(-4.0*(x(i)**2+y(j)**2+0.0d0*z(k)**2))
     v(i,j,k)=0.10d0 +0.75d0*exp(-4.0d0*(x(i)**2+y(j)**2+0.0d0*z(k)**2))
    END DO; END DO; END DO	
	END SUBROUTINE initialdata
