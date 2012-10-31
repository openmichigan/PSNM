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
	INTEGER(kind=4), PARAMETER 	::  Nx=2**8	
	INTEGER(kind=4), PARAMETER 	::  Ny=2**8	
	INTEGER(kind=4), PARAMETER	::  Nt=20	
	INTEGER(kind=4), PARAMETER	::  plotgap=5	
	REAL(kind=8), PARAMETER		:: &
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER		::  Lx=2.0d0	 
	REAL(kind=8), PARAMETER		::  Ly=2.0d0	 
	REAL(kind=8), PARAMETER		::  Es=0.0d0	 
	REAL(kind=8)				::  dt=0.10d0/Nt		
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	::  kx,ky 
	REAL(kind=8),  	 DIMENSION(:), ALLOCATABLE	::  x,y 
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  unax,vnax,vnbx,potx 
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE::  vnay,vnby
	REAL(kind=8), 	 DIMENSION(:), ALLOCATABLE	::  time
	INTEGER(kind=4)				::  i,j,k,n,allocatestatus
	INTEGER(kind=4)			 	::  start, finish, count_rate
	INTEGER(kind=8), PARAMETER  ::  FFTW_IN_PLACE=8, FFTW_MEASURE=0,&
									FFTW_EXHAUSTIVE=8, FFTW_PATIENT=32,&
                   	 				FFTW_ESTIMATE=64                                   
	INTEGER(kind=8),PARAMETER   ::  FFTW_FORWARD=-1, FFTW_BACKWARD=1	
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE :: fftfx,fftbx,fftfy,fftby
	INTEGER(kind=8)				::  planfx,planbx,planfy,planby
   	CHARACTER*100			 	::  name_config
	    			   	
	ALLOCATE(kx(1:Nx),ky(1:Nx),x(1:Nx),y(1:Nx),unax(1:Nx,1:Ny),&
			vnax(1:Nx,1:Ny),vnbx(1:Nx,1:Ny),potx(1:Nx,1:Ny),fftfx(1:Nx),&
			fftbx(1:Nx),fftfy(1:Nx),fftby(1:Nx),vnay(1:Ny,1:Nx),&
			vnby(1:Ny,1:Nx),time(1:1+Nt/plotgap),stat=allocatestatus)	
	IF (allocatestatus .ne. 0) stop 
	PRINT *,'allocated memory'
		! set up ffts
	CALL dfftw_plan_dft_1d_(planfx,Nx,fftfx(1:Nx),fftbx(1:Nx),&
			FFTW_FORWARD,FFTW_ESTIMATE)
 	CALL dfftw_plan_dft_1d_(planbx,Nx,fftbx(1:Nx),fftfx(1:Nx),&
 			FFTW_BACKWARD,FFTW_ESTIMATE)
	CALL dfftw_plan_dft_1d_(planfy,Ny,fftfy(1:Ny),fftby(1:Ny),&
			FFTW_FORWARD,FFTW_ESTIMATE)
 	CALL dfftw_plan_dft_1d_(planby,Ny,fftby(1:Ny),fftfy(1:Ny),&
 			FFTW_BACKWARD,FFTW_ESTIMATE)	
	PRINT *,'Setup FFTs'
		
	! setup fourier frequencies
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
		x(i)=(-1.0d0+2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)) )*pi*Lx
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
		y(j)=(-1.0d0+2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)) )*pi*Ly
	END DO
	!$OMP END PARALLEL DO
	PRINT *,'Setup grid and fourier frequencies'
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		DO i=1,Nx
			unax(i,j)=exp(-1.0d0*(x(i)**2 +y(j)**2)) 
		END DO
	END DO
	!$OMP END PARALLEL DO
	name_config = 'uinitial.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,Ny
		DO i=1,Nx
			WRITE(11,*) abs(unax(i,j))**2
		END DO
	END DO
	CLOSE(11)
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		DO i=1,Nx
			CALL dfftw_execute_dft_(planfx,unax(i,j),vnax(i,j)) 
		END DO
	END DO
	!$OMP END PARALLEL DO
	vnay(1:Ny,1:Nx)=TRANSPOSE(vnax(1:Nx,1:Ny))
	! transform initial data and do first half time step
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,Nx
		CALL dfftw_execute_dft_(planfy,vnay(1:Ny,i),vnby(1:Ny,i)) 
		DO j=1,Ny
			vnby(j,i)=exp(0.5d0*dt*(kx(i)*kx(i) + ky(j)*ky(j))&
					*cmplx(0.0d0,1.0d0))*vnby(j,i)
		END DO
		CALL dfftw_execute_dft_(planby,vnby(j,i),vnay(j,i))
	END DO
	!$OMP END PARALLEL DO
	PRINT *,'Got initial data, starting timestepping'
	time(1)=0.0d0
	CALL system_clock(start,count_rate)	
	DO n=1,Nt					
		vnbx(1:Nx,1:Ny)=TRANSPOSE(vnay(1:Ny,1:Nx))/REAL(Ny,kind(0d0))		
		!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
		DO j=1,Ny
			CALL dfftw_execute_dft_(planbx,vnbx(1:Nx,j),unax(1:Nx,j))
			DO i=1,Nx
				unax(i,j)=unax(1:Nx,j)/REAL(Nx,kind(0d0))
				potx(i,j)=Es*unax(i,j)*conjg(unax(i,j))
				unax(i,j)=exp(cmplx(0.0d0,-1.0d0)*dt*potx(i,j))&
							*unax(i,j)
			END DO
			CALL dfftw_execute_dft_(planfx,unax(1:Nx,j),vnax(1:Nx,j))
		END DO
		!$OMP END PARALLEL DO 
		vnby(1:Ny,1:Nx)=TRANSPOSE(vnax(1:Nx,1:Ny))		
		!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
		DO i=1,Nx
			CALL dfftw_execute_dft_(planfy,vnby(1:Ny,i),vnay(1:Ny,i))
			DO j=1,Ny
				vnby(j,i)=exp(dt*(kx(i)*kx(i) + ky(j)*ky(j))&	
						*cmplx(0.0d0,1.0d0))*vnay(j,i)	
			END DO
			CALL dfftw_execute_dft_(planby,vnby(1:Ny,i),vnay(1:Ny,i))
		END DO
		!$OMP END PARALLEL DO
		IF (mod(n,plotgap)==0) then
			time(1+n/plotgap)=n*dt
			PRINT *,'time',n*dt
		END IF
	END DO		
	PRINT *,'Finished time stepping'
	CALL system_clock(finish,count_rate)
	PRINT*,'Program took ',REAL(finish-start)/REAL(count_rate),&
		'for Time stepping'
		
	! transform back final data and do another half time step
	vnbx(1:Nx,1:Ny)=transpose(vnay(1:Ny,1:Nx))/REAL(Ny,kind(0d0))				
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		CALL dfftw_execute_dft_(planbx,vnbx(1:Nx,j),unax(1:Nx,j))
		unax(1:Nx,j)=unax(1:Nx,j)/REAL(Nx,kind(0d0))
		potx(1:Nx,j)=Es*unax(1:Nx,j)*conjg(unax(1:Nx,j))
		unax(1:Nx,j)=exp(cmplx(0,-1)*dt*potx(1:Nx,j))*unax(1:Nx,j)
		CALL dfftw_execute_dft_(planfx,unax(1:Nx,j),vnax(1:Nx,j))
	END DO
	!$OMP END PARALLEL DO
	vnby(1:Ny,1:Nx)=TRANSPOSE(vnax(1:Nx,1:Ny))				
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(static)
	DO i=1,Nx
		CALL dfftw_execute_dft_(planfy,vnby(1:Ny,i),vnay(1:Ny,i))
		vnby(1:Ny,i)=exp(0.5d0*dt*(kx(i)*kx(i) + ky(1:Ny)*ky(1:Ny))&
					*cmplx(0,1))*vnay(1:Ny,i)		
		CALL dfftw_execute_dft_(planby,vnby(1:Ny,i),vnay(1:Ny,i))
	END DO
	!$OMP END PARALLEL DO
	vnbx(1:Nx,1:Ny)=TRANSPOSE(vnay(1:Ny,1:Nx))/REAL(Ny,kind(0d0))				
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		CALL dfftw_execute_dft_(planbx,vnbx(1:Nx,j),unax(1:Nx,j))
		unax(1:Nx,j)=unax(1:Nx,j)/REAL(Nx,kind(0d0))
	END DO	
	!$OMP END PARALLEL DO
	name_config = 'ufinal.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,Ny
		DO i=1,Nx
			WRITE(11,*) abs(unax(i,j))**2
		END DO
	END DO
	CLOSE(11)
		
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

	CALL dfftw_destroy_plan_(planbx)
	CALL dfftw_destroy_plan_(planfx)
	CALL dfftw_destroy_plan_(planby)
	CALL dfftw_destroy_plan_(planfy)
	CALL dfftw_cleanup_()

	DEALLOCATE(unax,vnax,vnbx,potx, vnay,vnby,stat=allocatestatus)	
	IF (allocatestatus .ne. 0) STOP 
	PRINT *,'Deallocated memory'
		
   	PRINT *,'Program execution complete'
	END PROGRAM main
