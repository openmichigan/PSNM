	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves nonlinear Schrodinger equation in 3 dimensions
	! i*u_t+Es*|u|^2u+u_{xx}+u_{yy}+u_{zz}=0
	! using a second order time spectral splitting scheme
	!
	! The boundary conditions are u(x=0,y,z)=u(2*Lx*\pi,y,z), 
	!	u(x,y=0,z)=u(x,y=2*Ly*\pi,z), u(x,y,z=0)=u(x,y,z=2*Lz*\pi)
	! The initial condition is u=exp(-x^2-y^2)
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  Nz				= number of modes in z - power of 2 for FFT
	!  Nt				= number of timesteps to take
	!  Tmax				= maximum simulation time
	!  plotgap			= number of timesteps between plots
	!  FFTW_IN_PLACE 	= value for FFTW input 
	!  FFTW_MEASURE 	= value for FFTW input
	!  FFTW_EXHAUSTIVE 	= value for FFTW input
	!  FFTW_PATIENT 	= value for FFTW input    
	!  FFTW_ESTIMATE 	= value for FFTW input
	!  FFTW_FORWARD     = value for FFTW input
	!  FFTW_BACKWARD	= value for FFTW input	
	!  pi = 3.14159265358979323846264338327950288419716939937510d0
	!  Lx				= width of box in x direction
	!  Ly				= width of box in y direction
	!  Lz				= width of box in z direction
	!  ES				= +1 for focusing and -1 for defocusing
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
	!  n				= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  start			= variable to record start time of program
	!  finish			= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  count			= keep track of information written to disk
	!  iol				= size of array to write to disk
	!  planfxyz			= Forward 3d fft plan 
	!  planbxyz			= Backward 3d fft plan
	!  dt				= timestep
	!  modescalereal	= Number to scale after backward FFT 
	!  ierr				= error code
	! .. Arrays ..
	!  unax 			= approximate solution
	!  vnax 			= Fourier transform of approximate solution
	!  potx 			= potential
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  x				= x locations
	!  y				= y locations
	!  time				= times at which save data
	!  name_config		= array to store filename for data to be saved    		
	!  fftfxy			= array to setup 2D Fourier transform
	!  fftbxy			= array to setup 2D Fourier transform
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
	! FFTW3	 -- Fast Fourier Transform in the West Library
	!			(http://www.fftw.org/)
	PROGRAM main
	! Declare variables
	IMPLICIT NONE					 
	INTEGER(kind=4), PARAMETER	::  Nx=2**5	
	INTEGER(kind=4), PARAMETER 	::  Ny=2**5	
	INTEGER(kind=4), PARAMETER 	::  Nz=2**5	
	INTEGER(kind=4), PARAMETER	::  Nt=50	
	INTEGER(kind=4), PARAMETER	::  plotgap=10 
	REAL(kind=8), PARAMETER	::&
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER	::  Lx=2.0d0,Ly=2.0d0,Lz=2.0d0	 
	REAL(kind=8), PARAMETER	::  Es=1.0d0	 
	REAL(kind=8)	::  dt=0.10d0/Nt		
	REAL(kind=8) :: modescalereal
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	::  kx,ky,kz 
	REAL(kind=8),  	 DIMENSION(:), ALLOCATABLE	::  x,y,z 
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	::  unax,vnax,potx 
	REAL(kind=8), 	 DIMENSION(:), ALLOCATABLE	::  time
	INTEGER(kind=4)	::  i,j,k,n,AllocateStatus,count,iol
		! timing 
	INTEGER(kind=4)	::  start, finish, count_rate
  		! fftw variables
	INTEGER(kind=8), PARAMETER	::  FFTW_IN_PLACE = 8, FFTW_MEASURE = 0, &
		FFTW_EXHAUSTIVE = 8, FFTW_PATIENT = 32, FFTW_ESTIMATE = 64
	INTEGER(kind=8),PARAMETER	::  FFTW_FORWARD = -1, FFTW_BACKWARD=1	
	INTEGER(kind=8)	::	planfxyz,planbxyz
	CHARACTER*100	::	name_config, number_file
	
	CALL system_clock(start,count_rate)	
	ALLOCATE(unax(1:Nx,1:Ny,1:Nz),vnax(1:Nx,1:Ny,1:Nz),potx(1:Nx,1:Ny,1:Nz),&
   			kx(1:Nx),ky(1:Ny),kz(1:Nz),x(1:Nx),y(1:Ny),z(1:Nz),&
   			time(1:1+Nt/plotgap),stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 
	PRINT *,'allocated space'
	modescalereal=1.0d0/REAL(Nx,KIND(0d0))
	modescalereal=modescalereal/REAL(Ny,KIND(0d0))
	modescalereal=modescalereal/REAL(Nz,KIND(0d0))

	! set up ffts
	CALL dfftw_plan_dft_3d_(planfxyz,Nx,Ny,Nz,unax(1:Nx,1:Ny,1:Nz),&
		vnax(1:Nx,1:Ny,1:Nz),FFTW_FORWARD,FFTW_ESTIMATE)
 	CALL dfftw_plan_dft_3d_(planbxyz,Nx,Ny,Nz,vnax(1:Nx,1:Ny,1:Nz),&
 		unax(1:Nx,1:Ny,1:Nz),FFTW_BACKWARD,FFTW_ESTIMATE)
		
	PRINT *,'Setup FFTs'
		
	! setup fourier frequencies and grid points
	DO i=1,1+Nx/2
		kx(i)= cmplx(0.0d0,1.0)*REAL(i-1,kind(0d0))/Lx  			
	END DO
	kx(1+Nx/2)=0.0d0
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO
  	DO i=1,Nx
		x(i)=(-1.0d0+2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)) )*pi*Lx
	END DO
	DO j=1,1+Ny/2
		ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly  			
	END DO
	ky(1+Ny/2)=0.0d0
	DO j = 1,Ny/2 -1
		ky(j+1+Ny/2)=-ky(1-j+Ny/2)
	END DO
  	DO j=1,Ny
		y(j)=(-1.0d0+2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)) )*pi*Ly
	END DO
	DO k=1,1+Nz/2
		kz(k)= cmplx(0.0d0,1.0d0)*REAL(k-1,kind(0d0))/Lz  			
	END DO
	kz(1+Nz/2)=0.0d0
	DO k = 1,Nz/2 -1
		kz(k+1+Nz/2)=-kz(1-k+Nz/2)
	END DO
  	DO k=1,Nz
		z(k)=(-1.0d0+2.0d0*REAL(k-1,kind(0d0))/REAL(Nz,kind(0d0)) )*pi*Lz
	END DO

	PRINT *,'Setup grid and fourier frequencies'

	DO k=1,Nz;	DO j=1,Ny; DO i=1,Nx
		unax(i,j,k)=exp(-1.0d0*(x(i)**2 +y(j)**2+z(k)**2))
	END DO; END DO; END DO
		
	name_config = 'uinitial.dat' 
	INQUIRE(iolength=iol) unax(1,1,1)
	OPEN(unit=11,FILE=name_config,form="unformatted", &
	     access="direct",recl=iol) 	
	count=1
	DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
		WRITE(11,rec=count) unax(i,j,k)
        count=count+1
	END DO; END DO; END DO
	CLOSE(11)

	CALL dfftw_execute_dft_(planfxyz,unax(1:Nx,1:Ny,1:Nz),vnax(1:Nx,1:Ny,1:Nz)) 
		
	PRINT *,'Got initial data, starting timestepping'
	time(1)=0
	DO n=1,Nt
		DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
				vnax(i,j,k)=exp(0.50d0*dt*&
					(kz(k)*kz(k) + kx(i)*kx(i) + ky(j)*ky(j))&
					*cmplx(0.0d0,1.0d0))*vnax(i,j,k)
		END DO; END DO; END DO
		CALL dfftw_execute_dft_(planbxyz,vnax(1:Nx,1:Ny,1:Nz),&
			unax(1:Nx,1:Ny,1:Nz)) 

		DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
			unax(i,j,k)=unax(i,j,k)*modescalereal
			potx(i,j,k)=Es*unax(i,j,k)*conjg(unax(i,j,k))
			unax(i,j,k)=exp(cmplx(0.0d0,-1.0d0)*dt*potx(i,j,k))&
					*unax(i,j,k)
		END DO; END DO; END DO	
		CALL dfftw_execute_dft_(planfxyz,unax(1:Nx,1:Ny,1:Nz),&
			vnax(1:Nx,1:Ny,1:Nz)) 

		DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
			vnax(i,j,k)=exp(0.5d0*dt*&
					(kx(i)*kx(i) + ky(j)*ky(j)+ kz(k)*kz(k))&
						*cmplx(0.0d0,1.0d0))*vnax(i,j,k)		
		END DO; END DO; END DO	
		IF (mod(n,plotgap)==0) THEN
			time(1+n/plotgap)=n*dt
			PRINT *,'time',n*dt
			CALL dfftw_execute_dft_(planbxyz,vnax(1:Nx,1:Ny,1:Nz),unax(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
				unax(i,j,k)=unax(i,j,k)*modescalereal
			END DO; END DO; END DO
			name_config='./data/u'
			WRITE(number_file,'(i0)') 10000000+1+n/plotgap
			ind=index(name_config,' ') -1
			name_config=name_config(1:ind)//numberfile
			ind=index(name_config,' ') -1
			name_config=name_config(1:ind)//'.dat'
			OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
			REWIND(11)
			DO j=1,Ny
				DO i=1,Nx
					WRITE(11,*) abs(unax(i,j))**2
				END DO
			END DO
			CLOSE(11)
	
		END IF
	END DO		
	PRINT *,'Finished time stepping'
		
		! transform back final data and do another half time step
	CALL system_clock(finish,count_rate)
	PRINT*,'Program took ',REAL(finish-start)/REAL(count_rate),'for execution'
		
	name_config = 'tdata.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,1+Nt/plotgap
		WRITE(11,*) time(j)
	END DO
	CLOSE(11)

	name_config = 'xcoord.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO i=1,Nx
		WRITE(11,*) x(i)
	END DO
	CLOSE(11)

	name_config = 'ycoord.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,Ny
		WRITE(11,*) y(j)
	END DO
	CLOSE(11)
		
	name_config = 'zcoord.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO k=1,Nz
		WRITE(11,*) z(k)
	END DO
	CLOSE(11)
	PRINT *,'Saved data'

	CALL dfftw_destroy_plan_(planbxyz)
	CALL dfftw_destroy_plan_(planfxyz)
	CALL dfftw_cleanup_()
		
	DEALLOCATE(unax,vnax,potx,&
   			kx,ky,kz,x,y,z,&
   			time,stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 
	
	PRINT *,'Program execution complete'
	END PROGRAM main