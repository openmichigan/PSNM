	SUBROUTINE getgrid(Nx,Ny,Lx,Ly,pi,name_config,x,y,kx,ky)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine gets grid points and fourier frequencies for a
	! pseudospectral simulation of the 2D nonlinear Klein-Gordon equation
	!
	! u_{tt}-u_{xx}+u_{yy}+u=Es*u^3
	!
	! The boundary conditions are u(x=0,y)=u(2*Lx*\pi,y), 
	!	u(x,y=0)=u(x,y=2*Ly*\pi)
	!
	! INPUT
	!
	! .. Scalars ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  pi				= 3.142....
	!  Lx				= width of box in x direction
	!  Ly				= width of box in y direction
	! OUTPUT
	!
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  x				= x locations
	!  y				= y locations
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
	IMPLICIT NONE					 
	USE omp_lib		 	   
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)							:: Nx,Ny
	REAL(kind=8), INTENT(IN)							:: Lx,Ly,pi
	REAL(KIND=8), DIMENSION(1:NX), INTENT(OUT) 			:: x
	REAL(KIND=8), DIMENSION(1:NY), INTENT(OUT) 			:: y
	COMPLEX(KIND=8), DIMENSION(1:NX), INTENT(OUT)		:: kx
	COMPLEX(KIND=8), DIMENSION(1:NY), INTENT(OUT)		:: ky
	CHARACTER*100, INTENT(OUT)							:: name_config
	INTEGER(kind=4)										:: i,j
		
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,1+Nx/2
		kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/Lx  			
	END DO
	!$OMP END PARALLEL DO
	kx(1+Nx/2)=0.0d0
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO
	!$OMP END PARALLEL DO
		
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,Nx
		x(i)=(-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*Lx
	END DO
	!$OMP END PARALLEL DO
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,1+Ny/2
		ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly  			
	END DO
	!$OMP END PARALLEL DO
	ky(1+Ny/2)=0.0d0
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j = 1,Ny/2 -1
		ky(j+1+Ny/2)=-ky(1-j+Ny/2)
	END DO
	!$OMP END PARALLEL DO
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		y(j)=(-1.0d0 + 2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)))*pi*Ly
	END DO
	!$OMP END PARALLEL DO
	! Save x grid points in text format
	name_config = 'xcoord.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO i=1,Nx
		WRITE(11,*) x(i)
	END DO
	CLOSE(11)
	! Save y grid points in text format
	name_config = 'ycoord.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,Ny
		WRITE(11,*) y(j)
	END DO
	CLOSE(11)

	
	END SUBROUTINE getgrid