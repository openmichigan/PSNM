	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves nonlinear Schrodinger equation in 2 dimensions
	! i*u_t+Es*|u|^2u+u_{xx}+u_{yy}=0
	! using a second order time spectral splitting scheme
	!
	! The boundary conditions are u(x=0,y)=u(2*Lx*\pi,y), 
	!	u(x,y=0)=u(x,y=2*Ly*\pi)
	! The initial condition is u=exp(-x^2-y^2)
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
	!  numthreads		= number of openmp threads
	!  ierr				= error return code
	!  start			= variable to record start time of program
	!  finish			= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  planfx			= Forward 1d fft plan in x
	!  planbx			= Backward 1d fft plan in x
	!  planfy			= Forward 1d fft plan in y
	!  planby			= Backward 1d fft plan in y
	!  dt				= timestep
	! .. Arrays ..
	!  u				= approximate solution
	!  v				= Fourier transform of approximate solution
	!  unax 			= temporary field
	!  vnax 			= temporary field
	!  vnbx 			= temporary field
	!  vnay 			= temporary field
	!  vnby 			= temporary field
	!  potx 			= potential
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  x				= x locations
	!  y				= y locations
	!  time				= times at which save data
	!  name_config		= array to store filename for data to be saved    		
	!  fftfx			= array to setup x Fourier transform
	!  fftbx			= array to setup x Fourier transform
	!  fftfy			= array to setup y Fourier transform
	!  fftby			= array to setup y Fourier transform
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
	! OpenMP library		
	PROGRAM main
	USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(kind=4), PARAMETER 	::  Nx=1024
	INTEGER(kind=4), PARAMETER 	::  Ny=1024	
	INTEGER(kind=4), PARAMETER	::  Nt=20		
	INTEGER(kind=4), PARAMETER	::  plotgap=5	
	REAL(kind=8), PARAMETER		::  &
	pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER		::  Lx=2.0d0		
	REAL(kind=8), PARAMETER		::  Ly=2.0d0		
	REAL(kind=8), PARAMETER		::  Es=1.0d0		
	REAL(kind=8)				::  dt=0.10d0/Nt	
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	::  kx			
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	::  ky			
	REAL(kind=8),  	 DIMENSION(:), ALLOCATABLE	::  x			
	REAL(kind=8),  	 DIMENSION(:), ALLOCATABLE	::  y			
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  unax,vnax,vnbx,potx 
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  vnay,vnby
	REAL(kind=8), 	 DIMENSION(:), ALLOCATABLE	::  time
	INTEGER(kind=4)				::  i,j,k,n,allocatestatus,ierr,ind
	INTEGER(kind=4)			 	::  start, finish, count_rate, numthreads
	INTEGER(kind=8), PARAMETER  ::  FFTW_IN_PLACE=8, FFTW_MEASURE=0,&
									FFTW_EXHAUSTIVE=8, FFTW_PATIENT=32,&
                   	 				FFTW_ESTIMATE=64                                   
	INTEGER(kind=8),PARAMETER   ::  FFTW_FORWARD=-1, FFTW_BACKWARD=1	
	INTEGER(kind=8)				::  planfxy,planbxy
   	CHARACTER*100			 	::  name_config,number_file
	
	numthreads=omp_get_max_threads() 
	PRINT *,'There are ',numthreads,' threads.'
	    			   	
	ALLOCATE(kx(1:Nx),ky(1:Nx),x(1:Nx),y(1:Nx),unax(1:Nx,1:Ny),&
			vnax(1:Nx,1:Ny),potx(1:Nx,1:Ny),time(1:1+Nt/plotgap),&
			stat=allocatestatus)	
	IF (allocatestatus .ne. 0) stop 
	PRINT *,'allocated memory'

	! set up multithreaded ffts
	CALL dfftw_init_threads_(ierr)
	PRINT *,'Initiated threaded FFTW'
	CALL dfftw_plan_with_nthreads_(numthreads)
	PRINT *,'Inidicated number of threads to be used in planning'
	CALL dfftw_plan_dft_2d_(planfxy,Nx,Ny,unax(1:Nx,1:Ny),vnax(1:Nx,1:Ny),&
							FFTW_FORWARD,FFTW_ESTIMATE)
 	CALL dfftw_plan_dft_2d_(planbxy,Nx,Ny,vnax(1:Nx,1:Ny),unax(1:Nx,1:Ny),&
 							FFTW_BACKWARD,FFTW_ESTIMATE)
	PRINT *,'Setup FFTs'
		
	! setup fourier frequencies
	!$OMP PARALLEL PRIVATE(i,j) 
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
		x(i)=(-1.0d0+2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)) )*pi*Lx
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
		y(j)=(-1.0d0+2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)) )*pi*Ly
	END DO
	!$OMP END DO
	PRINT *,'Setup grid and fourier frequencies'
	!$OMP DO SCHEDULE(static)
	DO j=1,Ny
		unax(1:Nx,j)=exp(-1.0d0*(x(1:Nx)**2 +y(j)**2)) 
	END DO
	!$OMP END DO
	!$OMP END PARALLEL
	name_config = 'uinitial.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,Ny
		DO i=1,Nx
			WRITE(11,*) abs(unax(i,j))**2
		END DO
	END DO
	CLOSE(11)
	! transform initial data and do first half time step
	CALL dfftw_execute_dft_(planfxy,unax(1:Nx,1:Ny),vnax(1:Nx,1:Ny)) 

	PRINT *,'Got initial data, starting timestepping'
	time(1)=0.0d0
	CALL system_clock(start,count_rate)	
	DO n=1,Nt					
		!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
		DO j=1,Ny
			DO i=1,Nx
				vnax(i,j)=exp(0.5d0*dt*(kx(i)*kx(i) + ky(j)*ky(j))&
					*cmplx(0.0d0,1.0d0))*vnax(i,j)
			END DO
		END DO
		!$OMP END PARALLEL DO
		CALL dfftw_execute_dft_(planbxy,vnax(1:Nx,1:Ny),unax(1:Nx,1:Ny))
		!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
		DO j=1,Ny
			DO i=1,Nx
				unax(i,j)=unax(i,j)/REAL(Nx*Ny,kind(0d0))
				potx(i,j)=Es*unax(i,j)*conjg(unax(i,j))
				unax(i,j)=exp(cmplx(0.0d0,-1.0d0)*dt*potx(i,j))&
						*unax(i,j)
			END DO
		END DO
		!$OMP END PARALLEL DO 
		CALL dfftw_execute_dft_(planfxy,unax(1:Nx,1:Ny),vnax(1:Nx,1:Ny))
		!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
		DO j=1,Ny
			DO i=1,Nx
				vnax(i,j)=exp(0.5d0*dt*(kx(i)*kx(i) + ky(j)*ky(j))&	
						*cmplx(0.0d0,1.0d0))*vnax(i,j)
			END DO	
		END DO
		!$OMP END PARALLEL DO
		IF (mod(n,plotgap)==0) then
			time(1+n/plotgap)=n*dt
			PRINT *,'time',n*dt
			CALL dfftw_execute_dft_(planbxy,vnax(1:Nx,1:Ny),unax(1:Nx,1:Ny))
			!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
			DO j=1,Ny
				DO i=1,Nx
					unax(i,j)=unax(i,j)/REAL(Nx*Ny,kind(0d0))
				END DO
			END DO
			!$OMP END PARALLEL DO 
			name_config='./data/u'
			WRITE(number_file,'(i0)') 10000000+1+n/plotgap
			ind=index(name_config,' ') -1
			name_config=name_config(1:ind)//number_file
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
	CALL system_clock(finish,count_rate)
	PRINT*,'Program took ',REAL(finish-start)/REAL(count_rate),&
		'for Time stepping'
		
		
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
	PRINT *,'Saved data'

	CALL dfftw_destroy_plan_(planbxy)
	CALL dfftw_destroy_plan_(planfxy)
	CALL dfftw_cleanup_threads_()

	DEALLOCATE(unax,vnax,potx,stat=allocatestatus)	
	IF (allocatestatus .ne. 0) STOP 
	PRINT *,'Deallocated memory'
		
   	PRINT *,'Program execution complete'
	END PROGRAM main
