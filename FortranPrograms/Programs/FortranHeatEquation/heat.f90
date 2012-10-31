	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves heat equation in 1 dimension
	! u_t=\alpha*u_{xx}
	! using a the backward Euler method for x\in[0,2\pi]
	!
	! The boundary conditions are u(0)=u(2\pi)
	! The initial condition is u=sin(x)
	!
	! .. Parameters ..
	!  Nx	= number of modes in x - power of 2 for FFT
	!  Nt	= number of timesteps to take
	!  Tmax	= maximum simulation time
	!  plotgap			= number of timesteps between plots
	!  FFTW_IN_PLACE 	= value for FFTW input 
	!  FFTW_MEASURE 	= value for FFTW input
	!  FFTW_EXHAUSTIVE 	= value for FFTW input
	!  FFTW_PATIENT 	= value for FFTW input    
	!  FFTW_ESTIMATE 	= value for FFTW input
	!  FFTW_FORWARD     = value for FFTW input
	!  FFTW_BACKWARD	= value for FFTW input	
	!  pi = 3.14159265358979323846264338327950288419716939937510d0
	!  L				= width of box 
	!  alpha			= heat conductivity
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  n				= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  start			= variable to record start time of program
	!  finish			= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  planfx			= Forward 1d fft plan in x
	!  planbx			= Backward 1d fft plan in x
	!  dt				= timestep
	! .. Arrays ..
	!  u				= approximate REAL solution
	!  v				= Fourier transform of approximate solution
	!  vna 				= temporary field
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  x				= x locations
	!  time				= times at which save data
	!  name_config		= array to store filename for data to be saved    		
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
	INTEGER(kind=4),	PARAMETER 	::  Nx=64 
	INTEGER(kind=4),	PARAMETER	::  Nt=20  
	REAL(kind=8),	PARAMETER	&
		::  pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8),	PARAMETER	::  L=5.0d0			 
	REAL(kind=8),	PARAMETER	::  alpha=0.50d0	
	REAL(kind=8)	::  dt=0.2d0/REAL(Nt,kind(0d0))		
	COMPLEX(KIND=8), DIMENSION(:),ALLOCATABLE	::  kx	
	REAL(kind=8), DIMENSION(:),ALLOCATABLE	::  x	
	COMPLEX(KIND=8), DIMENSION(:,:),ALLOCATABLE	::  u,v	
	REAL(kind=8), DIMENSION(:),ALLOCATABLE	::  time
	COMPLEX(KIND=8), DIMENSION(:),ALLOCATABLE	::  vna	
	INTEGER(kind=4)	::  i,j,k,n
	INTEGER(kind=4)	:: start, finish, count_rate, AllocateStatus
	INTEGER(kind=4), PARAMETER	:: FFTW_IN_PLACE = 8, FFTW_MEASURE = 0, &
		FFTW_EXHAUSTIVE = 8, FFTW_PATIENT = 32, FFTW_ESTIMATE = 64
	INTEGER(kind=4), PARAMETER :: FFTW_FORWARD = -1, FFTW_BACKWARD=1	
	COMPLEX(KIND=8), DIMENSION(:),ALLOCATABLE :: fftfx,fftbx
	INTEGER(kind=8)	:: planfx,planbx
	CHARACTER*100	:: name_config
	    		
	CALL system_clock(start,count_rate)	
	ALLOCATE(kx(1:Nx),x(1:Nx),u(1:Nx,1:1+Nt),v(1:Nx,1:1+Nt),&
   			time(1:1+Nt),vna(1:Nx),fftfx(1:Nx),fftbx(1:Nx),&
   			stat=AllocateStatus)
   	IF (AllocateStatus .ne. 0) STOP 

		! set up ffts
	CALL dfftw_plan_dft_1d(planfx,Nx,fftfx(1:Nx),fftbx(1:Nx),&
			FFTW_FORWARD,FFTW_ESTIMATE)
	CALL dfftw_plan_dft_1d(planbx,Nx,fftbx(1:Nx),fftfx(1:Nx),&
			FFTW_BACKWARD,FFTW_ESTIMATE)
		
	PRINT *,'Setup FFTs'
		
	! setup fourier frequencies
	DO i=1,1+Nx/2
		kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/L  			
	END DO
	kx(1+Nx/2)=0.00d0
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO
	DO i=1,Nx
		x(i)=(-1.00d0 + 2.00d0*REAL(i-1,kind(0d0))/REAL(Nx,KIND(0d0)))*pi*L
	END DO
		
	PRINT *,'Setup grid and fourier frequencies and splitting coefficients'
	
	u(1:Nx,1)=sin(x(1:Nx)) 
		! transform initial data
	CALL dfftw_execute_dft_(planfx,u(1:Nx,1),v(1:Nx,1)) 
	PRINT *,'Got initial data, starting timestepping'
	time(1)=0.0d0

	vna(1:Nx)=v(1:Nx,1)
	PRINT *,'Starting timestepping'
	DO n=1,Nt
		DO i=1,Nx
			vna(i)=vna(i)/(1-dt*kx(i)*kx(i))
		END DO
		PRINT *,'storing plot data ',n
		time(n+1)=time(n)+dt
		v(1:Nx,n+1)=vna(1:Nx)
		CALL dfftw_execute_dft_(planbx,v(1:Nx,n+1),u(1:Nx,n+1))
		u(1:Nx,n+1)=u(1:Nx,n+1)/REAL(Nx,KIND(0d0))	! normalize
	END DO
	PRINT *,'Finished time stepping'
	CALL system_clock(finish,count_rate)
	PRINT*,'Program took ',REAL(finish-start)/REAL(count_rate),'for execution'

	! Write data out to disk	

	name_config = 'u.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,1+Nt
		DO i=1,Nx
			WRITE(11,*) REAL(u(i,j))
		END DO
	END DO
	CLOSE(11)
		
	name_config = 'tdata.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,1+Nt
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

	PRINT *,'Saved data'
	DEALLOCATE(kx,x,u,v,&
   			time,vna,fftfx,fftbx,&
   			stat=AllocateStatus)
   	IF (AllocateStatus .ne. 0) STOP 
	
	CALL dfftw_destroy_plan(planbx)
	CALL dfftw_destroy_plan(planfx)
	CALL dfftw_cleanup()
	PRINT *,'Program execution complete'
	END PROGRAM main

