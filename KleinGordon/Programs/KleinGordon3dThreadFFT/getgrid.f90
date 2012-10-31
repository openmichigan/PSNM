	SUBROUTINE getgrid(Nx,Ny,Nz,Lx,Ly,Lz,pi,name_config,x,y,z,kx,ky,kz)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine gets grid points and fourier frequencies for a
	! pseudospectral simulation of the 2D nonlinear Klein-Gordon equation
	!
	! u_{tt}-(u_{xx}+u_{yy}+u_{zz})+u=Es*u^3
	!
	! The boundary conditions are u(x=-Lx*pi,y,z)=u(x=Lx*\pi,y,z), 
	!	u(x,y=-Ly*pi,z)=u(x,y=Ly*pi,z),u(x,y,z=-Ly*pi)=u(x,y,z=Ly*pi),
	!
	! INPUT
	!
	! .. Scalars ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  Ny				= number of modes in z - power of 2 for FFT
	!  pi				= 3.142....
	!  Lx				= width of box in x direction
	!  Ly				= width of box in y direction
	!  Lz				= width of box in z direction
	! OUTPUT
	!
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  kz				= fourier frequencies in z direction
	!  x				= x locations
	!  y				= y locations
	!  z				= z locations
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
	! 	OpenMP library
	USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)							:: Nx,Ny,Nz
	REAL(kind=8), INTENT(IN)							:: Lx,Ly,Lz,pi
	REAL(KIND=8), DIMENSION(1:NX), INTENT(OUT) 			:: x
	REAL(KIND=8), DIMENSION(1:NY), INTENT(OUT) 			:: y
	REAL(KIND=8), DIMENSION(1:NZ), INTENT(OUT) 			:: z
	COMPLEX(KIND=8), DIMENSION(1:NX), INTENT(OUT)		:: kx
	COMPLEX(KIND=8), DIMENSION(1:NY), INTENT(OUT)		:: ky
	COMPLEX(KIND=8), DIMENSION(1:NZ), INTENT(OUT)		:: kz
	CHARACTER*100, INTENT(OUT)							:: name_config
	INTEGER(kind=4)										:: i,j,k
	
	!$OMP PARALLEL PRIVATE(i,j,k) 
	!$OMP DO SCHEDULE(static)
	DO i=1,1+Nx/2
		kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/Lx  			
	END DO
	!$OMP END DO
	kx(1+Nx/2)=0.0d0
	!$OMP DO SCHEDULE(static)
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO
	!$OMP END DO
	
	!$OMP DO SCHEDULE(static)
	DO i=1,Nx
		x(i)=(-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*Lx
	END DO
	!$OMP END DO

	!$OMP DO SCHEDULE(static)
	DO j=1,1+Ny/2
		ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly  			
	END DO
	!$OMP END DO
	ky(1+Ny/2)=0.0d0
	!$OMP DO SCHEDULE(static)
	DO j = 1,Ny/2 -1
		ky(j+1+Ny/2)=-ky(1-j+Ny/2)
	END DO
	!$OMP END DO
	!$OMP DO SCHEDULE(static)
	DO j=1,Ny
		y(j)=(-1.0d0 + 2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)))*pi*Ly
	END DO
	!$OMP END DO
	
	!$OMP DO SCHEDULE(static)
	DO k=1,1+Nz/2
		kz(k)= cmplx(0.0d0,1.0d0)*REAL(k-1,kind(0d0))/Lz  			
	END DO
	!$OMP END DO
	kz(1+Nz/2)=0.0d0
	!$OMP DO SCHEDULE(static)
	DO k = 1,Nz/2 -1
		kz(k+1+Nz/2)=-kz(1-k+Nz/2)
	END DO
	!$OMP END DO
	!$OMP DO SCHEDULE(static)
	DO k=1,Nz
		z(k)=(-1.0d0 + 2.0d0*REAL(k-1,kind(0d0))/REAL(Nz,kind(0d0)))*pi*Lz
	END DO
	!$OMP END DO
	!$OMP END PARALLEL
	
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
	! Save z grid points in text format
	name_config = 'zcoord.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO k=1,Nz
		WRITE(11,*) z(k)
	END DO
	CLOSE(11)
	
	END SUBROUTINE getgrid