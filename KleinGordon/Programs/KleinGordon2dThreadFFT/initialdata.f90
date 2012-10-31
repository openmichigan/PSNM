	SUBROUTINE initialdata(Nx,Ny,x,y,u,uold)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine gets initial data for nonlinear Klein-Gordon equation
	! in 2 dimensions
	! u_{tt}-u_{xx}+u_{yy}+u=Es*u^3+
	!
	! The boundary conditions are u(x=-Lx*\pi,y)=u(x=Lx*\pi,y), 
	!	u(x,y=-Ly*\pi)=u(x,y=Ly*\pi)
	! The initial condition is u=0.5*exp(-x^2-y^2)*sin(10*x+12*y)
	!
	! INPUT
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	! .. Vectors ..
	!  x				= x locations
	!  y				= y locations
	!
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
	! OpenMP library
	USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)							:: Nx,Ny
	REAL(KIND=8), DIMENSION(1:NX), INTENT(IN) 			:: x
	REAL(KIND=8), DIMENSION(1:NY), INTENT(IN) 			:: y
	COMPLEX(KIND=8), DIMENSION(1:NX,1:NY), INTENT(OUT)	:: u,uold
	INTEGER(kind=4)										:: i,j
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		u(1:Nx,j)=0.5d0*exp(-1.0d0*(x(1:Nx)**2 +y(j)**2))*&
					sin(10.0d0*x(1:Nx)+12.0d0*y(j))
	END DO
	!$OMP END PARALLEL DO
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		uold(1:Nx,j)=0.5d0*exp(-1.0d0*(x(1:Nx)**2 +y(j)**2))*&
					sin(10.0d0*x(1:Nx)+12.0d0*y(j))
	END DO
	!$OMP END PARALLEL DO
	
	END SUBROUTINE initialdata