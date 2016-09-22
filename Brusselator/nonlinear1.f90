SUBROUTINE nonlinear1(dt,coef,u,v,decomp)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine solves the first nonlinear equation in a splitting 
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
	!  u 				= updated u field
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
	!  temp             = temporary storage variable
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
   											   INTENT(INOUT) :: u
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   							   decomp%xst(2):decomp%xen(2),&
   							   decomp%xst(3):decomp%xen(3)),&
   											INTENT(IN)     :: v
	INTEGER(kind=4)										   :: i,j,k
	COMPLEX(kind=8)                                        :: temp
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     temp=u(i,j,k)
	 u(i,j,k)=temp/(1-temp*v(i,j,k)*coef*dt)
	END DO; END DO; END DO
	END SUBROUTINE nonlinear1
