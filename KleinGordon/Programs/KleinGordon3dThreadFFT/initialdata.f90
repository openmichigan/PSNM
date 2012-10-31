	SUBROUTINE initialdata(Nx,Ny,Nz,x,y,z,u,uold)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine gets initial data for nonlinear Klein-Gordon equation
	! in 3 dimensions
	! u_{tt}-(u_{xx}+u_{yy}+u_{zz})+u=Es*u^3+
	!
	! The boundary conditions are u(x=-Lx*\pi,y,z)=u(x=Lx*\pi,y,z), 
	!	u(x,y=-Ly*\pi,z)=u(x,y=Ly*\pi,z),u(x,y,z=-Ly*\pi)=u(x,y,z=Ly*\pi)
	! The initial condition is u=0.5*exp(-x^2-y^2-z^2)*sin(10*x+12*y)
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
	! OpenMP library
	USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)							:: Nx,Ny,Nz
	REAL(KIND=8), DIMENSION(1:NX), INTENT(IN) 			:: x
	REAL(KIND=8), DIMENSION(1:NY), INTENT(IN) 			:: y
	REAL(KIND=8), DIMENSION(1:NZ), INTENT(IN) 			:: z
	COMPLEX(KIND=8), DIMENSION(1:NX,1:NY,1:NZ), INTENT(OUT)	:: u,uold
	INTEGER(kind=4)										:: i,j,k
	!$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(static)
	DO k=1,Nz
		DO j=1,Ny
			DO i=1,Nx
				u(i,j,k)=0.5d0*exp(-1.0d0*(x(i)**2 +y(j)**2+z(k)**2))!*&
					!sin(10.0d0*x(i)+12.0d0*y(j))
			END DO
		END DO
	END DO
	!$OMP END PARALLEL DO
	!$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(static)
	DO k=1,Nz
		DO j=1,Ny
			DO i=1,Nx
				uold(i,j,k)=0.5d0*exp(-1.0d0*(x(i)**2 +y(j)**2+z(k)**2))!*&
					!sin(10.0d0*x(i)+12.0d0*y(j))
			END DO
		END DO
	END DO
	!$OMP END PARALLEL DO
	
	END SUBROUTINE initialdata