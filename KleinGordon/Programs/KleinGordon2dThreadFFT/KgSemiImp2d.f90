	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves nonlinear Klein-Gordon equation in 2 dimensions
	! u_{tt}-u_{xx}+u_{yy}+u=Es*|u|^2u
	! using a second order implicit-explicit time stepping scheme.
	!
	! The boundary conditions are u(x=0,y)=u(2*Lx*\pi,y), 
	!	u(x,y=0)=u(x,y=2*Ly*\pi)
	! The initial condition is u=0.5*exp(-x^2-y^2)*sin(10*x+12*y)
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
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
	!  ES				= +1 for focusing and -1 for defocusing
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  n				= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  start			= variable to record start time of program
	!  finish			= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  planfxy			= Forward 2d fft plan 
	!  planbxy			= Backward 2d fft plan
	!  dt				= timestep
	!  ierr				= error code
	!  plotnum			= number of plot
	! .. Arrays ..
	!  unew 			= approximate solution
	!  vnew 			= Fourier transform of approximate solution
	!  u 			= approximate solution
	!  v 			= Fourier transform of approximate solution
	!  uold 			= approximate solution
	!  vold 			= Fourier transform of approximate solution
	!  nonlin 			= nonlinear term, u^3
	!  nonlinhat 		= Fourier transform of nonlinear term, u^3
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  x				= x locations
	!  y				= y locations
	!  time				= times at which save data
	!  en				= total energy	
	!  enstr			= strain energy
	!  enpot			= potential energy
	!  enkin			= kinetic energy
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
	!	getgrid.f90	-- Get initial grid of points
	! 	initialdata.f90 -- Get initial data
	!	enercalc.f90 -- Subroutine to calculate the energy
	!	savedata.f90 -- Save initial data
	!	storeold.f90 -- Store old data
	! External libraries required
	! 	FFTW3	 -- Fast Fourier Transform in the West Library
	!			(http://www.fftw.org/)
	! 	OpenMP library
		
	PROGRAM Kg
	USE omp_lib		 	   
	! Declare variables
	IMPLICIT NONE					 
	INTEGER(kind=4), PARAMETER	:: Nx=128
	INTEGER(kind=4), PARAMETER 	:: Ny=128
	INTEGER(kind=4), PARAMETER	:: Nt=20 
	INTEGER(kind=4), PARAMETER	:: plotgap=5	
	REAL(kind=8), PARAMETER		:: &
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER		::  Lx=3.0d0
	REAL(kind=8), PARAMETER		::  Ly=3.0d0 
	REAL(kind=8), PARAMETER		::  Es=1.0d0	
	REAL(kind=8)				::  dt=0.10d0/REAL(Nt,kind(0d0))	
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	::  kx,ky 
	REAL(kind=8),  	 DIMENSION(:), ALLOCATABLE	::  x,y	
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  u,nonlin
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  v,nonlinhat
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  uold
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  vold
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  unew
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  vnew
	REAL(kind=8), DIMENSION(:,:), ALLOCATABLE :: savearray
	REAL(kind=8), DIMENSION(:), ALLOCATABLE	::  time,enkin,enstr,enpot,en
	INTEGER(kind=4)	::  ierr,i,j,n,allocatestatus,numthreads
	INTEGER(kind=4)	::  start, finish, count_rate, plotnum
	INTEGER(kind=4), PARAMETER	::  FFTW_IN_PLACE = 8, FFTW_MEASURE = 0, &
		FFTW_EXHAUSTIVE = 8, FFTW_PATIENT = 32, FFTW_ESTIMATE = 64
	INTEGER(kind=4),PARAMETER	::  FFTW_FORWARD = -1, FFTW_BACKWARD=1	
	INTEGER(kind=8)	::  planfxy,planbxy
	CHARACTER*100	::  name_config
	! Start short parallel region to count threads 
	numthreads=omp_get_max_threads()
	PRINT *,'There are ',numthreads,' threads.'
	ALLOCATE(kx(1:Nx),ky(1:Ny),x(1:Nx),y(1:Ny),u(1:Nx,1:Ny),&
			v(1:Nx,1:Ny),nonlin(1:Nx,1:Ny),nonlinhat(1:Nx,1:Ny),&
			uold(1:Nx,1:Ny),vold(1:Nx,1:Ny),&
			unew(1:Nx,1:Ny),vnew(1:Nx,1:Ny),savearray(1:Nx,1:Ny),&
			time(1:1+Nt/plotgap),enkin(1:1+Nt/plotgap),&
			enstr(1:1+Nt/plotgap),enpot(1:1+Nt/plotgap),&
			en(1:1+Nt/plotgap),stat=allocatestatus)	
	IF (allocatestatus .ne. 0) stop 
	PRINT *,'allocated arrays'
		
	! set up multithreaded ffts
	CALL dfftw_init_threads_(ierr)
	PRINT *,'Initiated threaded FFTW'
	CALL dfftw_plan_with_nthreads_(numthreads)
	PRINT *,'Inidicated number of threads to be used in planning'
	CALL dfftw_plan_dft_2d_(planfxy,Nx,Ny,u,v,&
							FFTW_FORWARD,FFTW_ESTIMATE)
 	CALL dfftw_plan_dft_2d_(planbxy,Nx,Ny,v,u,&
 							FFTW_BACKWARD,FFTW_ESTIMATE)
	PRINT *,'Setup FFTs'
	! setup fourier frequencies
	CALL getgrid(Nx,Ny,Lx,Ly,pi,name_config,x,y,kx,ky)
	PRINT *,'Setup grid and fourier frequencies'
	CALL initialdata(Nx,Ny,x,y,u,uold)
	plotnum=1	
	name_config = 'data/u' 
	savearray=REAL(u)
	CALL savedata(Nx,Ny,plotnum,name_config,savearray)

	CALL dfftw_execute_dft_(planfxy,u,v) 
	CALL dfftw_execute_dft_(planfxy,uold,vold) 

	CALL enercalc(Nx,Ny,planfxy,planbxy,dt,Es,&
				enkin(plotnum),enstr(plotnum),&
				enpot(plotnum),en(plotnum),&
				kx,ky,nonlin,nonlinhat,&
				v,vold,u,uold)
			
	PRINT *,'Got initial data, starting timestepping'
	time(plotnum)=0.0d0
  	CALL system_clock(start,count_rate)	
	DO n=1,Nt					
		!$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(static)
		DO j=1,Ny
			DO i=1,Nx
				nonlin(i,j)=(abs(u(i,j))*2)*u(i,j)
			END DO
		END DO
		!$OMP END PARALLEL DO
		CALL dfftw_execute_dft_(planfxy,nonlin,nonlinhat)
		!$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(static)
		DO j=1,Ny
			DO i=1,Nx
				vnew(i,j)=( 0.25*(kx(i)*kx(i) + ky(j)*ky(j)-1.0d0)&
					*(2.0d0*v(i,j)+vold(i,j))&
					+(2.0d0*v(i,j)-vold(i,j))/(dt*dt)&
					+Es*nonlinhat(i,j) )&
					/(1/(dt*dt)-0.25*(kx(i)*kx(i) + ky(j)*ky(j)-1.0d0))
			END DO
		END DO
		!$OMP END PARALLEL DO
		CALL dfftw_execute_dft_(planbxy,vnew,unew)
		! normalize result
		!$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(static)
		DO j=1,Ny
			DO i=1,Nx
				unew(i,j)=unew(i,j)/REAL(Nx*Ny,kind(0d0))
			END DO
		END DO
		!$OMP END PARALLEL DO		
		IF (mod(n,plotgap)==0) then
			plotnum=plotnum+1
			time(plotnum)=n*dt
			PRINT *,'time',n*dt
			CALL enercalc(Nx,Ny,planfxy,planbxy,dt,Es,&
				enkin(plotnum),enstr(plotnum),&
				enpot(plotnum),en(plotnum),&
				kx,ky,&
				nonlin,nonlinhat,&
				vnew,v,unew,u)
			savearray=REAL(unew,kind(0d0))
			CALL savedata(Nx,Ny,plotnum,name_config,savearray)
		END IF
			! .. Update old values ..
		CALL storeold(Nx,Ny,&
				unew,u,uold,&
				vnew,v,vold)
	END DO		
	PRINT *,'Finished time stepping'
	CALL system_clock(finish,count_rate)
	PRINT*,'Program took ',&
		REAL(finish-start,kind(0d0))/REAL(count_rate,kind(0d0)),&
		'for Time stepping'
	CALL saveresults(Nt,plotgap,time(1:1+n/plotgap),en(1:1+n/plotgap),&
			enstr(1:1+n/plotgap),enkin(1:1+n/plotgap),enpot(1:1+n/plotgap))	
			
	! Save times at which output was made in text format
	PRINT *,'Saved data'

	CALL dfftw_destroy_plan_(planbxy)
	CALL dfftw_destroy_plan_(planfxy)
	CALL dfftw_cleanup_threads_()
	DEALLOCATE(kx,ky,x,y,u,v,nonlin,nonlinhat,savearray,&
		uold,vold,unew,vnew,time,enkin,enstr,enpot,en,&
		stat=allocatestatus)	
	IF (allocatestatus .ne. 0) STOP 
	PRINT *,'Deallocated arrays'
	PRINT *,'Program execution complete'
	END PROGRAM Kg
