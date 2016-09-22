SUBROUTINE nonlinear2(dt,coef,u,v,decomp)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine solves the second nonlinear equation in a splitting 
	! method for the Brusselator equation in 3 dimensions
	! u_t = a + u^2v - (b+1)u + Du u_xx
    ! v_t = bu - u^2v + Dv v_xx
	!
	!
	! INPUT
	!
	! .. Parameters ..
	!  dt				= timestep
	!  coef				= (complex) timestep coeffcient
	! .. Arrays ..
	!  u				= u field
	!  v				= v field
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!					see http://www.2decomp.org/ for more information
	! OUTPUT
	!
	! .. Arrays ..
	!  v 				= updated v field
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
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
	TYPE(DECOMP_INFO), INTENT(IN)						   :: decomp
	REAL(KIND=8), INTENT(IN)                               :: dt
	COMPLEX(KIND=8),  INTENT(IN)                           :: coef
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   							   decomp%xst(2):decomp%xen(2),&
   							   decomp%xst(3):decomp%xen(3)),&
   											   INTENT(IN) :: u
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   							   decomp%xst(2):decomp%xen(2),&
   							   decomp%xst(3):decomp%xen(3)),&
   											INTENT(INOUT)  :: v
	INTEGER(kind=4)										   :: i,j,k
    DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	 v(i,j,k)=v(i,j,k)*exp(-1.0d0*dt*coef*u(i,j,k)**2)
	END DO; END DO; END DO	
	END SUBROUTINE nonlinear2
