	SUBROUTINE getgrid(myid,Nx,Ny,Nz,Lx,Ly,Lz,pi,name_config,x,y,z,kx,ky,kz,decomp)
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
	!  myid				= processor id
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!					see http://www.2decomp.org/ for more information
	!
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
	! 2DECOMP&FFT	 -- Domain decomposition and Fast Fourier Library
	!			(http://www.2decomp.org/index.html)
	! MPI library
	USE decomp_2d
	IMPLICIT NONE					 
	INCLUDE 'mpif.h'
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)							:: myid,Nx,Ny,Nz
	REAL(kind=8), INTENT(IN)							:: Lx,Ly,Lz,pi
	TYPE(DECOMP_INFO), INTENT(IN)						::  decomp
	REAL(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1)), INTENT(OUT) 	:: x
	REAL(KIND=8), DIMENSION(decomp%xst(2):decomp%xen(2)), INTENT(OUT) 	:: y
	REAL(KIND=8), DIMENSION(decomp%xst(3):decomp%xen(3)), INTENT(OUT) 	:: z
	COMPLEX(KIND=8), DIMENSION(decomp%zst(1):decomp%zen(1)), INTENT(OUT):: kx
	COMPLEX(KIND=8), DIMENSION(decomp%zst(2):decomp%zen(2)), INTENT(OUT):: ky
	COMPLEX(KIND=8), DIMENSION(decomp%zst(3):decomp%zen(3)), INTENT(OUT):: kz
	CHARACTER*100, INTENT(OUT)							:: name_config
	INTEGER(kind=4)										:: i,j,k
	
	
	DO i = 1,1+ Nx/2
		IF ((i.GE.decomp%zst(1)).AND.(i.LE.decomp%zen(1))) THEN
			kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/Lx
		END IF
	END DO
	IF ((Nx/2 + 1 .GE.decomp%zst(1)).AND.(Nx/2 + 1 .LE.decomp%zen(1))) THEN
		kx( Nx/2 + 1 ) = 0.0d0 
	ENDIF
	DO i = Nx/2+2, Nx  
		IF ((i.GE.decomp%zst(1)).AND.(i.LE.decomp%zen(1))) THEN
			Kx( i) = cmplx(0.0d0,-1.0d0)*REAL(1-i+Nx,KIND(0d0))/Lx  
		ENDIF
	END DO      
	DO i=decomp%xst(1),decomp%xen(1)
		x(i)=(-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*Lx
	END DO

	DO j = 1,1+ Ny/2
		IF ((j.GE.decomp%zst(2)).AND.(j.LE.decomp%zen(2))) THEN
			ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly
		END IF
	END DO
	IF ((Ny/2 + 1 .GE.decomp%zst(2)).AND.(Ny/2 + 1 .LE.decomp%zen(2))) THEN
		ky( Ny/2 + 1 ) = 0.0d0 
	ENDIF
	DO j = Ny/2+2, Ny  
		IF ((j.GE.decomp%zst(2)).AND.(j.LE.decomp%zen(2))) THEN
			ky(j) = cmplx(0.0d0,-1.0d0)*REAL(1-j+Ny,KIND(0d0))/Ly  
		ENDIF
	END DO      
	DO j=decomp%xst(2),decomp%xen(2)
		y(j)=(-1.0d0 + 2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)))*pi*Ly
	END DO
	
	DO k = 1,1+ Nz/2
		IF ((k.GE.decomp%zst(3)).AND.(k.LE.decomp%zen(3))) THEN
			kz(k)= cmplx(0.0d0,1.0d0)*REAL(k-1,kind(0d0))/Lz
		END IF
	END DO
	IF ((Nz/2 + 1 .GE.decomp%zst(3)).AND.(Nz/2 + 1 .LE.decomp%zen(3))) THEN
		kz( Nz/2 + 1 ) = 0.0d0 
	ENDIF
	DO k = Nz/2+2, Nz  
		IF ((k.GE.decomp%zst(3)).AND.(k.LE.decomp%zen(3))) THEN
			kz(k) = cmplx(0.0d0,-1.0d0)*REAL(1-k+Nz,KIND(0d0))/Lz  
		ENDIF
	END DO      
	DO k=decomp%xst(3),decomp%xen(3)
		z(k)=(-1.0d0 + 2.0d0*REAL(k-1,kind(0d0))/REAL(Nz,kind(0d0)))*pi*Lz
	END DO
	
	IF (myid.eq.0) THEN
		! Save x grid points in text format
		name_config = 'xcoord.dat' 
		OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
		REWIND(11)
		DO i=1,Nx
			WRITE(11,*) (-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*Lx
		END DO
		CLOSE(11)
		! Save y grid points in text format
		name_config = 'ycoord.dat' 
		OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
		REWIND(11)
		DO j=1,Ny
			WRITE(11,*) (-1.0d0 + 2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)))*pi*Ly
		END DO
		CLOSE(11)
		! Save z grid points in text format
		name_config = 'zcoord.dat' 
		OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
		REWIND(11)
		DO k=1,Nz
			WRITE(11,*) (-1.0d0 + 2.0d0*REAL(k-1,kind(0d0))/REAL(Nz,kind(0d0)))*pi*Lz
		END DO
		CLOSE(11)
	END IF
	
	END SUBROUTINE getgrid
