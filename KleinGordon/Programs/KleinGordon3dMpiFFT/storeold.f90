	SUBROUTINE storeold(Nx,Ny,Nz,unew,u,uold,vnew,v,vold,decomp)
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
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!					see http://www.2decomp.org/ for more information
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
	! 2DECOMP&FFT	 -- Domain decomposition and Fast Fourier Library
	!			(http://www.2decomp.org/index.html)
	! MPI library
	IMPLICIT NONE					 
	USE decomp_2d
	USE decomp_2d_fft
	USE decomp_2d_io
	INCLUDE 'mpif.h'	
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)								:: Nx,Ny,Nz
	TYPE(DECOMP_INFO), INTENT(IN)							::  decomp
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   								decomp%xst(2):decomp%xen(2),&
   								decomp%xst(3):decomp%xen(3)), INTENT(OUT):: uold
	COMPLEX(KIND=8), DIMENSION(decomp%zst(1):decomp%zen(1),&
            					decomp%zst(2):decomp%zen(2),&
            					decomp%zst(3):decomp%zen(3)), INTENT(OUT):: vold
	COMPLEX(KIND=8), DIMENSION(decomp%zst(1):decomp%zen(1),&
            					decomp%zst(2):decomp%zen(2),&
            					decomp%zst(3):decomp%zen(3)), INTENT(INOUT):: v
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   								decomp%xst(2):decomp%xen(2),&
   								decomp%xst(3):decomp%xen(3)), INTENT(INOUT):: u
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   								decomp%xst(2):decomp%xen(2),&
   								decomp%xst(3):decomp%xen(3)), INTENT(IN):: unew
	COMPLEX(KIND=8), DIMENSION(decomp%zst(1):decomp%zen(1),&
            					decomp%zst(2):decomp%zen(2),&
            					decomp%zst(3):decomp%zen(3)), INTENT(IN):: vnew
	INTEGER(kind=4)														:: i,j,k

	DO k=decomp%zst(3),decomp%zen(3)
		DO j=decomp%zst(2),decomp%zen(2)
			DO i=decomp%zst(1),decomp%zen(1)
				vold(i,j,k)=v(i,j,k)
			END DO
		END DO
	END DO
	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				uold(i,j,k)=u(i,j,k)
			END DO
		END DO
	END DO
	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				u(i,j,k)=unew(i,j,k)
			END DO
		END DO
	END DO
	DO k=decomp%zst(3),decomp%zen(3)
		DO j=decomp%zst(2),decomp%zen(2)
			DO i=decomp%zst(1),decomp%zen(1)
				v(i,j,k)=vnew(i,j,k)
			END DO
		END DO
	END DO
	
	END SUBROUTINE storeold
