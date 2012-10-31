	SUBROUTINE storeold(Nx,Ny,Nz,unew,u,uold,vnew,v,vold)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine copies arrays for a
	! pseudospectral simulation of the 2D nonlinear Klein-Gordon equation
	!
	! u_{tt}-(u_{xx}+u_{yy}+u_{zz})+u=Es*u^3
	!
	! INPUT
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  Nz				= number of modes in z - power of 2 for FFT
	!  .. Arrays ..
	!  unew 			= approximate solution
	!  vnew 			= Fourier transform of approximate solution
	!  u 				= approximate solution
	!  v 				= Fourier transform of approximate solution
	!  uold 			= approximate solution
	!  vold 			= Fourier transform of approximate solution
	!
	! OUTPUT
	!
	!  u 				= approximate solution
	!  v 				= Fourier transform of approximate solution
	!  uold 			= approximate solution
	!  vold 			= Fourier transform of approximate solution
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
	!--------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	! OpenMP library
	USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)								:: Nx,Ny,Nz
	COMPLEX(KIND=8), DIMENSION(1:NX,1:NY,1:NZ), INTENT(OUT)	:: vold,uold
	COMPLEX(KIND=8), DIMENSION(1:NX,1:NY,1:NZ), INTENT(INOUT):: u,v
	COMPLEX(KIND=8), DIMENSION(1:NX,1:NY,1:NZ), INTENT(IN)	:: unew,vnew
	INTEGER(kind=4)											:: i,j,k

	!$OMP PARALLEL PRIVATE(j,k)
	
	!$OMP DO SCHEDULE(static)
	DO k=1,Nz
		DO j=1,Ny
			DO i=1,Nx
				vold(i,j,k)=v(i,j,k)
			END DO
		END DO
	END DO
	!$OMP END DO NOWAIT

	!$OMP DO SCHEDULE(static)
	DO k=1,Nz
		DO j=1,Ny
			DO i=1,Nx
				uold(i,j,k)=u(i,j,k)
			END DO
		END DO
	END DO
	!$OMP END DO NOWAIT
	
	!$OMP DO SCHEDULE(static)
	DO k=1,Nz
		DO j=1,Ny
			DO i=1,Nx
				u(i,j,k)=unew(i,j,k)
			END DO
		END DO
	END DO
	!$OMP END DO NOWAIT
	
	!$OMP DO SCHEDULE(static)
	DO k=1,Nz
		DO j=1,Ny
			DO i=1,Nx
				v(i,j,k)=vnew(i,j,k)
			END DO
		END DO
	END DO
	!$OMP END DO NOWAIT
	
	!$OMP END PARALLEL 
	
	END SUBROUTINE storeold
