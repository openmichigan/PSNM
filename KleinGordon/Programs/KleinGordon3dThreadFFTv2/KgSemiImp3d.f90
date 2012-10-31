	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves nonlinear Klein-Gordon equation in 3 dimensions
	! u_{tt}-(u_{xx}+u_{yy}+u_{zz})+u=Es*|u|^2u
	! using a second order implicit-explicit time stepping scheme.
	!
	! The boundary conditions are u(x=-Lx*pi,y,z)=u(x=Lx*\pi,y,z), 
	!	u(x,y=-Ly*pi,z)=u(x,y=Ly*pi,z),u(x,y,z=-Ly*pi)=u(x,y,z=Ly*pi),
	! The initial condition is u=0.5*exp(-x^2-y^2-z^2)*sin(10*x+12*y)
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
	!  planfz			= Forward 1d fft plan 
	!  planbz			= Backward 1d fft plan
	!  planfxy			= Forward 2d fft plan 
	!  planbxy			= Backward 2d fft plan
	!  planfxyz			= Forward 3d fft plan 
	!  planbxyz			= Backward 3d fft plan
	!  dt				= timestep
	!  modescalereal	= Number to scale after backward FFT 
	!  ierr				= error code
	!  plotnum			= number of plot
	!  numthreads		= total number of OpenMP threads
	!  mythreadid		= thread number
	! .. Arrays ..
	!  unew 			= approximate solution
	!  vnew 			= Fourier transform of approximate solution
	!  u 				= approximate solution
	!  v 				= Fourier transform of approximate solution
	!  uold 			= approximate solution
	!  vold 			= Fourier transform of approximate solution
	!  nonlin 			= nonlinear term, u^3
	!  nonlinhat 		= Fourier transform of nonlinear term, u^3
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  kz				= fourier frequencies in z direction
	!  x				= x locations
	!  y				= y locations
	!  z				= z locations
	!  time				= times at which save data
	!  en				= total energy	
	!  enstr			= strain energy
	!  enpot			= potential energy
	!  enkin			= kinetic energy
	!  name_config		= array to store filename for data to be saved 
	!  lowxv			= array with entries of array to be used by each thread
	!  upxv				= array with entries of array to be used by each thread
	!  lowzv			= array with entries of array to be used by each thread
	!  upzv				= array with entries of array to be used by each thread
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
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(kind=4), PARAMETER	:: Nx=64
	INTEGER(kind=4), PARAMETER 	:: Ny=64
	INTEGER(kind=4), PARAMETER 	:: Nz=64
	INTEGER(kind=4), PARAMETER	:: Nt=20 
	INTEGER(kind=4), PARAMETER	:: plotgap=10	
	REAL(kind=8), PARAMETER		:: &
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER		::  Lx=64.0d0
	REAL(kind=8), PARAMETER		::  Ly=64.0d0 
	REAL(kind=8), PARAMETER		::  Lz=64.0d0 
	REAL(kind=8), PARAMETER		::  Es=1.0d0	
	REAL(kind=8)				::  dt=0.050d0/REAL(Nt,kind(0d0))	
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	::  kx,ky,kz 
	REAL(kind=8),  	 DIMENSION(:), ALLOCATABLE	::  x,y,z	
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  u,nonlin,vtemp1
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  v,nonlinhat,vtemp2
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  uold
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  vold
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  unew
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  vnew
	REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: savearray
	REAL(kind=8), DIMENSION(:), ALLOCATABLE	::  time,enkin,enstr,enpot,en
	INTEGER(kind=4), DIMENSION(:), ALLOCATABLE :: lowxv, upxv, lowzv, upzv
	REAL(kind=8) :: modescalereal
	INTEGER(kind=4)	::  ierr,i,j,k,n,allocatestatus,numthreads,mythreadid
	INTEGER(kind=4)	::  start, finish, count_rate, plotnum
	INTEGER(kind=4), PARAMETER	::  FFTW_IN_PLACE = 8, FFTW_MEASURE = 0, &
		FFTW_EXHAUSTIVE = 8, FFTW_PATIENT = 32, FFTW_ESTIMATE = 64
	INTEGER(kind=4),PARAMETER	::  FFTW_FORWARD = -1, FFTW_BACKWARD=1	
	INTEGER(kind=8)	::  planfxy,planbxy,planfz,planbz
	CHARACTER*100	::  name_config
	! Start short parallel region to count threads 
	!$OMP PARALLEL
    numthreads=omp_get_num_threads()
    !$OMP END PARALLEL
    PRINT *,'There are ',numthreads,' threads.'
	ALLOCATE(kx(1:Nx),ky(1:Ny),kz(1:Nz),x(1:Nx),y(1:Ny),z(1:Nz),u(1:Nx,1:Ny,1:Nz),&
			v(1:Nz,1:Ny,1:Nx),nonlin(1:Nx,1:Ny,1:Nz),nonlinhat(1:Nz,1:Ny,1:Nx),&
			uold(1:Nx,1:Ny,1:Nz),vold(1:Nz,1:Ny,1:Nx),vtemp1(1:Nx,1:Ny,1:Nz),&
			vtemp2(1:Nz,1:Ny,1:Nx),unew(1:Nx,1:Ny,1:Nz),vnew(1:Nz,1:Ny,1:Nx),&
			savearray(1:Nx,1:Ny,1:Nz),time(1:1+Nt/plotgap),enkin(1:1+Nt/plotgap),&
			enstr(1:1+Nt/plotgap),enpot(1:1+Nt/plotgap),en(1:1+Nt/plotgap),&
			lowxv(1:numthreads),upxv(1:numthreads),lowzv(1:numthreads),upzv(1:numthreads),&
			stat=allocatestatus)	
	IF (allocatestatus .ne. 0) stop 
	PRINT *,'allocated arrays'
	
	! set up ffts
	CALL dfftw_plan_dft_2d_(planfxy,Nx,Ny,u(1:Nx,1:Ny,1),vtemp1(1:Nx,1:Ny,1),&
							FFTW_FORWARD,FFTW_ESTIMATE)
 	CALL dfftw_plan_dft_2d_(planbxy,Nx,Ny,vtemp1(1:Nx,1:Ny,1),u(1:Nx,1:Ny,1),&
 							FFTW_BACKWARD,FFTW_ESTIMATE)
	CALL dfftw_plan_dft_(planfxy,Nx,Ny,v(1:Nz,1:Ny,1),vold(1:Nx,1:Ny,1),&
							FFTW_FORWARD,FFTW_ESTIMATE)
 	CALL dfftw_plan_dft_(planbz,Nz,vold(1:Nz,1,1),v(1:Nz,1,1),&
 							FFTW_BACKWARD,FFTW_ESTIMATE)
	PRINT *,'Setup FFTs'
	DO i=1,numthreads
		lowxv(i)=1+(i-1)*Nx/numthreads
	END DO
	DO i=1,numthreads
		upxv(i)=i*Nx/numthreads
	END DO
	DO k=1,numthreads
		lowzv(k)=1+(k-1)*Nz/numthreads
	END DO
	DO k=1,numthreads
		upxv(k)=k*Nz/numthreads
	END DO
	upxv(numthreads)=Nx
	upzv(numthreads)=Nz
	
	! setup fourier frequencies
	CALL getgrid(Nx,Ny,Nz,Lx,Ly,Lz,pi,name_config,x,y,z,kx,ky,kz)
	PRINT *,'Setup grid and fourier frequencies'
	!$OMP PARALLEL PRIVATE(i,j,k)
	! setup array ranges
	mythreadid=omp_get_thread_num()
	CALL initialdata(Nx,Ny,Nz,mythreadid,numthreads,lowxv,lowzv,upxv,upzv,&
				x,y,z,u,uold)
	!$OMP MASTER			
	plotnum=1	
	name_config = 'data/u' 
	!$OMP END MASTER
	DO k=lowzv(mythreadid),upzv(mythreadid)
		DO j=1,Ny
			DO i=1,Nx
				savearray(i,j,k)=REAL(unew(i,j,k),kind(0d0))
			END DO
		END DO
	END DO	
	!$OMP MASTER
	CALL savedata(Nx,Ny,Nz,plotnum,name_config,savearray)
	!$OMP END MASTER
	! forward transform
	CALL forwardfft(Nx,Ny,Nz,planfz,planfxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,u,v,vtemp1,vtemp2)					
	CALL forwardfft(Nx,Ny,Nz,planfz,planfxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,uold,vold,vtemp1,vtemp2)
	!$OMP MASTER
	modescalereal=1.0d0/REAL(Nx,KIND(0d0))
	modescalereal=modescalereal/REAL(Ny,KIND(0d0))
	modescalereal=modescalereal/REAL(Nz,KIND(0d0))
	!$OMP END MASTER

	CALL enercalc(Nx,Ny,Nz,planbxy,planbz,planfxy,planfz,dt,Es,&
			modescalereal,enkin(plotnum),enstr(plotnum),enpot(plotnum),&
			en(plotnum),mythreadid,numthreads,lowxv,lowzv,upxv,upzv,&
			kx,ky,kz,nonlin,nonlinhat,vtemp1,vtemp2,v,vold,u,uold)
	!$OMP MASTER
	PRINT *,'Got initial data, starting timestepping'
	time(plotnum)=0.0d0
  	CALL system_clock(start,count_rate)	
	!$OMP END MASTER
	DO n=1,Nt					
		DO k=lowzv(mythreadid),upzv(mythreadid)
			DO j=1,Ny
				DO i=1,Nx
					nonlin(i,j,k)=(abs(u(i,j,k))*2)*u(i,j,k)
				END DO
			END DO
		END DO
		CALL forwardfft(Nx,Ny,Nz,planfz,planfxy,mythreadid,numthreads,&
					   lowxv,lowzv,upxv,upzv,nonlin,nonlinhat,vtemp1,vtemp2)

		DO i=lowxv(mythreadid),upxv(mythreadid)
			DO j=1,Ny
				DO k=1,Nz
					vnew(k,j,i)=&
					( 0.25d0*(kx(i)*kx(i) + ky(j)*ky(j)+ kz(k)*kz(k)-1.0d0)&
					*(2.0d0*v(k,j,i)+vold(k,j,i))&
					+(2.0d0*v(k,j,i)-vold(k,j,i))/(dt*dt)&
					+Es*nonlinhat(k,j,i) )&
					/(1/(dt*dt)-0.25*(kx(i)*kx(i)+ ky(j)*ky(j)+ kz(k)*kz(k)-1.0d0))
				END DO	
			END DO
		END DO
		CALL backwardfft(Nx,Ny,Nz,planbz,planbxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,vnew,unew,vtemp1,vtemp2)
										
		! normalize result
		DO k=lowzv(mythreadid),upzv(mythreadid)
			DO j=1,Ny
				DO i=1,Nx
					unew(i,j,k)=unew(i,j,k)*modescalereal
				END DO
			END DO
		END DO
		IF (mod(n,plotgap)==0) then
			!$OMP MASTER
			plotnum=plotnum+1
			time(plotnum)=n*dt
			PRINT *,'time',n*dt
			!$OMP END MASTER
			CALL enercalc(Nx,Ny,Nz,planbxy,planbz,planfxy,planfz,dt,Es,&
				modescalereal,enkin(plotnum),enstr(plotnum),enpot(plotnum),&
				en(plotnum),mythreadid,numthreads,lowxv,lowzv,upxv,upzv,kx,ky,kz,&
				nonlin,nonlinhat,vtemp1,vtemp2,v,vold,u,uold)
			DO k=lowzv(mythreadid),upzv(mythreadid)
				DO j=1,Ny
					DO i=1,Nx
						savearray(i,j,k)=REAL(unew(i,j,k),kind(0d0))
					END DO
				END DO
			END DO	
			!$OMP BARRIER
			!$OMP MASTER
			CALL savedata(Nx,Ny,Nz,plotnum,name_config,savearray)
			!$OMP END MASTER
		END IF
			! .. Update old values ..
		!$OMP BARRIER	
		CALL storeold(Nx,Ny,Nz,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,unew,u,uold,vnew,v,vold)
		!$OMP BARRIER			
	END DO	
	!$OMP END PARALLEL
	PRINT *,'Finished time stepping'
	CALL system_clock(finish,count_rate)
	PRINT*,'Program took ',&
		REAL(finish-start,kind(0d0))/REAL(count_rate,kind(0d0)),&
		'for Time stepping'
	CALL saveresults(Nt,plotgap,time,en,enstr,enkin,enpot)	
			
	! Save times at which output was made in text format
	PRINT *,'Saved data'

	CALL dfftw_destroy_plan_(planbxy)
	CALL dfftw_destroy_plan_(planfxy)
	CALL dfftw_destroy_plan_(planbz)
	CALL dfftw_destroy_plan_(planfz)
	CALL dfftw_cleanup()
	DEALLOCATE(kx,ky,kz,x,y,z,u,&
			v,nonlin,nonlinhat,&
			uold,vold,vtemp1,&
			vtemp2,unew,vnew,&
			savearray,time,enkin,&
			enstr,enpot,en,&
			lowxv,upxv,lowzv,upzv,&
			stat=allocatestatus)	
	IF (allocatestatus .ne. 0) STOP 
	PRINT *,'Deallocated arrays'
	PRINT *,'Program execution complete'
	END PROGRAM Kg
