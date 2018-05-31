	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves nonlinear Schrodinger equation in 1 dimension
	! i*u_t+Es*|u|^2u+u_{xx}=0
	! using a second order time spectral splitting scheme
	!
	! The boundary conditions are u(0)=u(2*L*\pi)
	! The initial condition is u=exp(-x^2)
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
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
	!  L				= width of box 
	!  ES				= +1 for focusing and -1 for defocusing
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
	!  u				= approximate solution
	!  v				= Fourier transform of approximate solution
	! .. Vectors ..
	!  una 				= temporary field
	!  unb 				= temporary field
	!  vna 				= temporary field
	!  pot 				= potential
	!  kx				= fourier frequencies in x direction
	!  x				= x locations
	!  time				= times at which save data
	!  name_config		= array to store filename for data to be saved    		
	!  fftfx			= array to setup x Fourier transform
	!  fftbx			= array to setup x Fourier transform
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
	INTEGER(kind=4), PARAMETER 	:: Nx=8*256		
	INTEGER(kind=4), PARAMETER	:: Nt=200	 
	REAL(kind=8), PARAMETER		&
		::  pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER	:: L=5.0d0	 
	REAL(kind=8), PARAMETER	:: Es=1.0d0	
	REAL(kind=8)	:: dt=2.0d0/Nt		
	COMPLEX(kind=8),  DIMENSION(:), ALLOCATABLE	:: kx	
	REAL(kind=8),  DIMENSION(:), ALLOCATABLE	:: x	
	COMPLEX(kind=8),  DIMENSION(:,:), ALLOCATABLE	:: u	
	COMPLEX(kind=8),  DIMENSION(:,:), ALLOCATABLE	:: v 	
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	:: una,vn 
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	:: unb,pot 
	REAL(kind=8), DIMENSION(:), ALLOCATABLE		:: time
	INTEGER(kind=4)	:: i,j,k,n,modes,AllocateStatus
	INTEGER(kind=4)	:: start, finish, count_rate
	INTEGER(kind=4), PARAMETER	:: FFTW_IN_PLACE = 8, FFTW_MEASURE = 0, &
		FFTW_EXHAUSTIVE = 8, FFTW_PATIENT = 32, FFTW_ESTIMATE = 64                                   
	INTEGER(kind=4), PARAMETER  :: FFTW_FORWARD = -1, FFTW_BACKWARD=1	
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	:: fftfx,fftbx
	INTEGER(kind=8)	:: planfx,planbx
  	CHARACTER*100	:: name_config
	    		
	CALL system_clock(start,count_rate)
	ALLOCATE(kx(1:Nx),x(1:Nx),u(1:Nx,1:Nt+1),v(1:Nx,1:Nt+1),&
		una(1:Nx),vn(1:Nx),unb(1:Nx),pot(1:Nx),time(1:Nt+1),&
		fftfx(1:Nx),fftbx(1:Nx),stat=AllocateStatus)
	IF (allocatestatus .ne. 0) STOP		
	! set up ffts
	CALL dfftw_plan_dft_1d_(planfx,Nx,fftfx(1:Nx),fftbx(1:Nx),&
		FFTW_FORWARD,FFTW_PATIENT)
	CALL dfftw_plan_dft_1d_(planbx,Nx,fftbx(1:Nx),fftfx(1:Nx),&
		FFTW_BACKWARD,FFTW_PATIENT)
	PRINT *,'Setup FFTs'
		! setup fourier frequencies
	DO i=1,1+Nx/2
		kx(i)= cmplx(0.0d0,1.0d0)*(i-1.0d0)/L  			
	END DO
	kx(1+Nx/2)=0.0d0
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO
	DO i=1,Nx
		x(i)=(-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*L
	END DO
	PRINT *,'Setup grid and fourier frequencies'
	
	DO i=1,Nx
		u(i,1)=exp(-1.0d0*(x(i)**2))
	END DO 
	! transform initial data
	CALL dfftw_execute_dft_(planfx,u(1:Nx,1),v(1:Nx,1)) 
	PRINT *,'Got initial data, starting timestepping'
	time(1)=0.0d0
	DO n=1,Nt
		time(n+1)=n*dt
		DO i=1,Nx
			vn(i)=exp(0.5d0*dt*kx(i)*kx(i)*cmplx(0.0d0,1.0d0))*v(i,n)
		END DO
		CALL dfftw_execute_dft_(planbx,vn(1:Nx),una(1:Nx))
		! normalize
		DO i=1,Nx
			una(i)=una(i)/REAL(Nx,kind(0d0))	
			pot(i)=Es*una(i)*conjg(una(i))
			unb(i)=exp(cmplx(0.0d0,-1.0d0)*dt*pot(i))*una(i)
		END DO
		CALL dfftw_execute_dft_(planfx,unb(1:Nx),vn(1:Nx))
		DO i=1,Nx
			v(i,n+1)=exp(0.50d0*dt*kx(i)*kx(i)*cmplx(0.0d0,1.0d0))*vn(i)
		END DO
		CALL dfftw_execute_dft_(planbx,v(1:Nx,n+1),u(1:Nx,n+1))
		! normalize
		DO i=1,Nx
			u(i,n+1)=u(i,n+1)/REAL(Nx,kind(0d0))
		END DO	
	END DO
	PRINT *,'Finished time stepping'
	CALL system_clock(finish,count_rate)
	PRINT*,'Program took ',&
		REAL(finish-start,kind(0d0))/REAL(count_rate,kind(0d0)),'for execution'
		
	name_config = 'u.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,Nt
		DO i=1,Nx
			WRITE(11,*) abs(u(i,j))**2
		END DO
	END DO
	CLOSE(11)
		
	name_config = 'tdata.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,Nt
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
		
	CALL dfftw_destroy_plan_(planbx)
	CALL dfftw_destroy_plan_(planfx)
	CALL dfftw_cleanup_()
		
	DEALLOCATE(kx,x,u,v,una,vn,unb,&
				pot,time,fftfx,fftbx,&
	   			stat=AllocateStatus)
	IF (allocatestatus .ne. 0) STOP 
	PRINT *,'deallocated memory'
	PRINT *,'Program execution complete'
	END PROGRAM main
