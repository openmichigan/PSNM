	SUBROUTINE backwardfft(Nx,Ny,Nz,planbz,planbxy,mythreadid,numthreads,&
							lowxv,lowzv,upxv,upzv,&
							arrayin,arrayout,vtemp1,vtemp2)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine does a forward FFT using 
	!
	!
	! INPUT
	!
	! .. Scalars ..
	!  Nx				= number of modes in x
	!  Ny				= number of modes in y
	!  Nz				= number of modes in z
	!  planfz			= plan for 1D forward FFT in z direction
	!  planfxy			= plan for 2D forward FFT in xy directions
	!  numthreads		= total number of OpenMP threads
	!  mythreadid		= thread number
	!  .. Arrays ..
	!  arrayin 			= input array
	!  vtemp1 			= temporary storage array
	!  vtemp2 			= temporary storage array
	!  .. Vectors ..
	!  lowxv			= array with entries of array to be used by each thread
	!  upxv				= array with entries of array to be used by each thread
	!  lowzv			= array with entries of array to be used by each thread
	!  upzv				= array with entries of array to be used by each thread
	!
	! OUTPUT
	!
	!  arrayout 		= output array
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in y direction
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
	INTEGER(KIND=4), INTENT(IN)		:: Nx,Ny,Nz,numthreads,mythreadid
	INTEGER(KIND=4), DIMENSION(1:numthreads), INTENT(IN) ::&
												lowxv,upxv,lowzv,upzv											
	INTEGER(KIND=8), INTENT(IN)								:: planbz,planbxy		
	COMPLEX(KIND=8), DIMENSION(1:Nz,1:Ny,1:Nx),INTENT(IN)	:: arrayin
	COMPLEX(KIND=8), DIMENSION(1:Nx,1:Ny,1:Nz),INTENT(OUT)	:: vtemp1
	COMPLEX(KIND=8), DIMENSION(1:Nx,1:Ny,1:Nz),INTENT(OUT)	:: arrayout
	COMPLEX(KIND=8), DIMENSION(1:Nz,1:Ny,1:Nx),INTENT(OUT)	:: vtemp2
	INTEGER(KIND=4)											:: i,j,k

	!$OMP DO SCHEDULE(static) 	
	DO i=1,Nx
		DO j=1,Ny
			CALL dfftw_execute_dft_(planbz,arrayin(1:Nz,j,i),vtemp2(1:Nz,j,i))
		END DO
	END DO
	!$OMP END DO
	!$OMP FLUSH
	!$OMP MASTER
	DO k=1,Nz
		DO j=1,Ny
			DO i=1,Nx
				vtemp1(i,j,k)=vtemp2(k,j,i)
			END DO
		END DO
	END DO
	!$OMP END MASTER
	!$OMP DO SCHEDULE(static) 
	DO k=1,Nz
		CALL dfftw_execute_dft_(planbxy,vtemp1(1:Nx,1:Ny,k),arrayout(1:Nx,1:Ny,k))
	END DO
	!$OMP END DO
	!$OMP FLUSH
	
	END SUBROUTINE backwardfft