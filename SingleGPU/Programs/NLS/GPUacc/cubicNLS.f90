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
	! AUTHORS
	!
	! B. Cloutier, B.K. Muite, P. Rigge
	! 4 June 2012
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
	!  planfxy			= Forward 2d fft plan 
	!  planbxy			= Backward 2d fft plan
	!  dt				= timestep
	!  InMass			= initial mass
	!  FiMass			= final mass
	!  InEner			= initial energy
	!  FiEner			= final energy
	! .. Arrays ..
	!  u				= approximate solution
	!  v				= Fourier transform of approximate solution
	!  temp1 			= temporary field
	!  temp2 			= temporary field
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  x				= x locations
	!  y				= y locations
	!  time				= times at which save data
	!  name_config		= array to store filename for data to be saved    		
	!
	! REFERENCES
	!
	! This program is based on example code to demonstrate usage of Fortran and 
	! CUDA FFT routines taken from 
	! http://cudamusing.blogspot.com/2010/05/CALLing-cufft-from-cuda-fortran.html
	! 
	! and
	!
	! http://cudamusing.blogspot.com/search?q=cublas
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
	! precision
	! cufft
	!
	! External libraries required
	! CuFFT	 -- Cuda FFT Library
	! OpenACC
	
	!
	! Define the INTERFACE to the NVIDIA CUFFT routines
	!

	module precision
	! Precision control

	integer, parameter, public :: Single = kind(0.0) ! Single precision
	integer, parameter, public :: Double = kind(0.0d0) ! Double precision

	integer, parameter, public :: fp_kind = Double
	!integer, parameter, public :: fp_kind = Single

	end module precision

	module cufft

	integer, public :: CUFFT_FORWARD = -1
	integer, public :: CUFFT_INVERSE = 1
	integer, public :: CUFFT_R2C = Z'2a' ! Real to Complex (interleaved)
	integer, public :: CUFFT_C2R = Z'2c' ! Complex (interleaved) to Real
	integer, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
	integer, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex
	integer, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double
	integer, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! cufftPlan2d(cufftHandle *plan, int nx,int ny, cufftType type)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	interface cufftPlan2d
	subroutine cufftPlan2d(plan, nx, ny, type) bind(C,name='cufftPlan2d')
	use iso_c_binding
	integer(c_int):: plan
	integer(c_int),value:: nx, ny, type
	end subroutine cufftPlan2d
	end interface cufftPlan2d

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! cufftDestroy(cufftHandle plan)
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	interface cufftDestroy
	subroutine cufftDestroy(plan) bind(C,name='cufftDestroy')
	use iso_c_binding
	integer(c_int),value:: plan
	end subroutine cufftDestroy
	end interface cufftDestroy

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!
	! cufftExecZ2Z(cufftHandle plan,
	! cufftDoubleComplex *idata,
	! cufftDoubleComplex *odata,
	! int direction;
	!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	interface cufftExecZ2Z
	subroutine cufftExecZ2Z(plan, idata, odata, direction) &
	& bind(C,name='cufftExecZ2Z')
	use iso_c_binding
	use precision
	integer(c_int),value:: direction
	integer(c_int),value:: plan
	complex(fp_kind),device,dimension(1:nx,1:ny):: idata,odata
	end subroutine cufftExecZ2Z
	end interface cufftExecZ2Z	
	end module cufft
		
	PROGRAM main
	USE precision
	USE cufft	
	USE openacc

	! Declare variables
	IMPLICIT NONE					 
	INTEGER(kind=4), PARAMETER 	::  Nx=128
	INTEGER(kind=4), PARAMETER 	::  Ny=128	
	INTEGER(kind=4), PARAMETER	::  Nt=20		
	INTEGER(kind=4), PARAMETER	::  plotgap=20	
	REAL(fp_kind), PARAMETER		::  &
	pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(fp_kind), PARAMETER		::  Lx=5.0d0		
	REAL(fp_kind), PARAMETER		::  Ly=5.0d0		
	REAL(fp_kind), PARAMETER		::  Es=1.0d0		
	REAL(fp_kind)					::  dt=0.10d0**5
	REAL(fp_kind)					:: scalemodes
	COMPLEX(fp_kind)				:: InMass,FiMass,InEner,FiEner 		
	COMPLEX(fp_kind), DIMENSION(:), ALLOCATABLE		::  kx			
	COMPLEX(fp_kind), DIMENSION(:), ALLOCATABLE		::  ky			
	REAL(fp_kind),  	 DIMENSION(:), ALLOCATABLE	::  x			
	REAL(fp_kind),  	 DIMENSION(:), ALLOCATABLE	::  y			
	COMPLEX(fp_kind), DIMENSION(:,:), ALLOCATABLE	::  u,v,temp1,temp2 
	REAL(fp_kind), 	 DIMENSION(:), ALLOCATABLE		::  time
	INTEGER(kind=4)				::  i,j,k,n,allocatestatus,ierr, vecsize,gangsize
	REAL(fp_kind)			 	::  start_time,stop_time
	INTEGER(kind=4)				::  plan
   	CHARACTER*100			 	::  name_config
	
	vecsize=32
	gangsize=16
	PRINT *,'Program starting'
	PRINT *,'Grid: ',Nx,'X',Ny	
	    			   	
	ALLOCATE(kx(1:Nx),ky(1:Nx),x(1:Nx),y(1:Nx),u(1:Nx,1:Ny),&
			v(1:Nx,1:Ny),temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny),&
			time(1:1+Nt/plotgap),stat=allocatestatus)	
	IF (allocatestatus .ne. 0) stop 
	PRINT *,'allocated memory'

	!$acc data copy(InMass,FiMass,InEner,FiEner,kx,ky,x,y,u,v,temp1,temp2,time)
	
	! set up ffts
	CALL cufftPlan2D(plan,nx,ny,CUFFT_Z2Z)
	PRINT *,'Setup FFTs'
		
	! setup fourier frequencies
	!$acc kernels loop	
	DO i=1,1+Nx/2
		kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/Lx  			
	END DO
	!$acc end kernels
	kx(1+Nx/2)=0.0d0
	!$acc kernels loop	
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO
	!$acc end kernels
	!$acc kernels loop	
  	DO i=1,Nx
		x(i)=(-1.0d0+2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)) )*pi*Lx
	END DO
	!$acc end kernels
	!$acc kernels loop	
	DO j=1,1+Ny/2
		ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly  			
	END DO
	!$acc end kernels
	ky(1+Ny/2)=0.0d0
	!$acc kernels loop	
	DO j = 1,Ny/2 -1
		ky(j+1+Ny/2)=-ky(1-j+Ny/2)
	END DO
	!$acc end kernels
	!$acc kernels loop	
  	DO j=1,Ny
		y(j)=(-1.0d0+2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)) )*pi*Ly
	END DO
	!$acc end kernels
	scalemodes=1.0d0/REAL(Nx*Ny,kind(0d0))
	PRINT *,'Setup grid and fourier frequencies'
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			u(i,j)=exp(-1.0d0*(x(i)**2 +y(j)**2)) 
		END DO
	END DO
	!$acc end kernels
	! transform initial data 
	CALL cufftExecZ2Z(plan,u,v,CUFFT_FORWARD)

	PRINT *,'Got initial data'
	! get initial mass
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp1(i,j)=abs(u(i,j))**2
		END DO 
	END DO
	!$acc end kernels
	! Use FFT to get initial mass
	CALL cufftExecZ2Z(plan,temp1,temp2,CUFFT_FORWARD)
	!$acc end data
	InMass=temp2(1,1)
	! Get initial energy
	!$acc data copy(InMass,FiMass,InEner,FiEner,kx,ky,x,y,u,v,temp1,temp2,time)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp1(i,j)=-ES*0.25d0*abs(u(i,j))**4
		END DO 
	END DO
	!$acc end kernels
	! Use FFT to find mean
	CALL cufftExecZ2Z(plan,temp1,temp2,CUFFT_FORWARD)
	!$acc end data
	InEner=temp2(1,1)
	!$acc data copy(InMass,FiMass,InEner,FiEner,kx,ky,x,y,u,v,temp1,temp2,time)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp2(i,j)=kx(i)*v(i,j)*scalemodes
		END DO 
	END DO
	!$acc end kernels
	CALL cufftExecZ2Z(plan,temp2,temp1,CUFFT_INVERSE)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp2(i,j)=0.5d0*abs(temp1(i,j))**2
		END DO 
	END DO
	!$acc end kernels
	! Use FFT to find mean
	CALL cufftExecZ2Z(plan,temp2,temp1,CUFFT_FORWARD)
	!$acc end data
	InEner=InEner+temp1(1,1)
	!$acc data copy(InMass,FiMass,InEner,FiEner,kx,ky,x,y,u,v,temp1,temp2,time)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp2(i,j)=ky(j)*v(i,j)*scalemodes
		END DO 
	END DO
	!$acc end kernels
	CALL cufftExecZ2Z(plan,temp2,temp1,CUFFT_INVERSE)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp2(i,j)=0.5d0*abs(temp1(i,j))**2
		END DO 
	END DO
	!$acc end kernels
	! Use FFT to find mean
	CALL cufftExecZ2Z(plan,temp2,temp1,CUFFT_FORWARD)
	!$acc end data
	InEner=InEner+temp1(1,1)
	!$acc data copy(InMass,FiMass,InEner,FiEner,kx,ky,x,y,u,v,temp1,temp2,time)
	CALL cpu_time(start_time)


	! transform initial data and do first half time step
	!$acc kernels loop gang(gangsize), vector(vecsize)	
	DO j=1,Ny
		DO i=1,Nx
			v(i,j)=exp(0.5d0*dt*(kx(i)*kx(i) + ky(j)*ky(j))&
					*cmplx(0.0d0,1.0d0))*v(i,j)
		END DO
	END DO
	!$acc end kernels
	PRINT *,'Got initial data, starting timestepping'
	time(1)=0.0d0
	DO n=1,Nt					
		CALL cufftExecZ2Z(plan,v,u,CUFFT_INVERSE)
	    !$acc kernels loop gang(gangsize), vector(vecsize)		
		DO j=1,Ny
			DO i=1,Nx
				v(i,j)=Es*u(i,j)*conjg(u(i,j))*scalemodes**2
			END DO
		END DO
	    !$acc end kernels
	    !$acc kernels loop gang(gangsize), vector(vecsize)	
		DO j=1,Ny
			DO i=1,Nx
				u(i,j)=exp(cmplx(0.0d0,-1.0d0)*dt*v(i,j))&
						*u(i,j)*scalemodes
			END DO
		END DO
	    !$acc end kernels
		CALL cufftExecZ2Z(plan,u,v,CUFFT_FORWARD)
	    !$acc kernels loop gang(gangsize), vector(vecsize)		
		DO j=1,Ny
			DO i=1,Nx
				v(i,j)=exp(dt*(kx(i)*kx(i) + ky(j)*ky(j))&	
						*cmplx(0.0d0,1.0d0))*v(i,j)	
			END DO
		END DO
	    !$acc end kernels
		IF (mod(n,plotgap)==0) then
			time(1+n/plotgap)=n*dt
			PRINT *,'time',n*dt
		END IF
	END DO				
	! transform back final data and do another half time step
	CALL cufftExecZ2Z(plan,v,u,CUFFT_INVERSE)
	!$acc kernels loop gang(gangsize), vector(vecsize)		
	DO j=1,Ny
		DO i=1,Nx
			v(i,j)=Es*u(i,j)*conjg(u(i,j))*scalemodes**2
		END DO
	END DO
	!$acc end kernels
	!$acc kernels loop gang(gangsize), vector(vecsize)	
	DO j=1,Ny
		DO i=1,Nx
			u(i,j)=exp(cmplx(0,-1)*dt*v(i,j))*u(i,j)*scalemodes
		END DO
	END DO
	!$acc end kernels
	CALL cufftExecZ2Z(plan,u,v,CUFFT_FORWARD)
	!$acc kernels loop gang(gangsize), vector(vecsize)	
	DO j=1,Ny
		DO i=1,Nx
			v(i,j)=exp(0.5d0*dt*(kx(i)*kx(i) + ky(j)*ky(j))&
					*cmplx(0.0d0,1.0d0))*v(i,j)
		END DO		
	END DO
	!$acc end kernels
	CALL cufftExecZ2Z(plan,v,u,CUFFT_INVERSE)
	!$acc kernels loop gang(gangsize), vector(vecsize)		
	DO j=1,Ny
		DO i=1,Nx
			u(i,j)=u(i,j)*scalemodes
		END DO
	END DO	
	!$acc end kernels
	PRINT *,'Finished time stepping'
	CALL cpu_time(stop_time)
	!$acc end data
	PRINT*,'Program took ',stop_time-start_time,&
		'for Time stepping'
	!$acc data copy(InMass,FiMass,InEner,FiEner,kx,ky,x,y,u,v,temp1,temp2,time)

	! calculate final mass
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp1(i,j)=abs(u(i,j))**2
		END DO 
	END DO
	!$acc end kernels
	! Use FFT to get initial mass
	CALL cufftExecZ2Z(plan,temp1,temp2,CUFFT_FORWARD)
	!$acc end data
	FiMass=temp2(1,1)


	! Get final energy
	!$acc data copy(InMass,FiMass,InEner,FiEner,kx,ky,x,y,u,v,temp1,temp2,time)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp1(i,j)=-ES*0.25d0*abs(u(i,j))**4
		END DO 
	END DO
	!$acc end kernels
	! Use FFT to find mean
	CALL cufftExecZ2Z(plan,temp1,temp2,CUFFT_FORWARD)
	!$acc end data
	FiEner=temp2(1,1)
	!$acc data copy(InMass,FiMass,InEner,FiEner,kx,ky,x,y,u,v,temp1,temp2,time)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp2(i,j)=kx(i)*v(i,j)*scalemodes
		END DO 
	END DO
	!$acc end kernels
	CALL cufftExecZ2Z(plan,temp2,temp1,CUFFT_INVERSE)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp2(i,j)=0.5d0*abs(temp1(i,j))**2
		END DO 
	END DO
	!$acc end kernels
	! Use FFT to find mean
	CALL cufftExecZ2Z(plan,temp2,temp1,CUFFT_FORWARD)
	!$acc end data
	FiEner=FiEner+temp1(1,1)
	!$acc data copy(InMass,FiMass,InEner,FiEner,kx,ky,x,y,u,v,temp1,temp2,time)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp2(i,j)=ky(j)*v(i,j)*scalemodes
		END DO 
	END DO
	!$acc end kernels
	CALL cufftExecZ2Z(plan,temp2,temp1,CUFFT_INVERSE)
	!$acc kernels loop	
	DO j=1,Ny
		DO i=1,Nx
			temp2(i,j)=0.5d0*abs(temp1(i,j))**2
		END DO 
	END DO
	!$acc end kernels
	! Use FFT to find mean
	CALL cufftExecZ2Z(plan,temp2,temp1,CUFFT_FORWARD)
	!$acc end data
	FiEner=FiEner+temp1(1,1)
	
 	PRINT *,'Results copied back to host'
	PRINT*,'Initial mass',InMass
	PRINT*,'Final mass',FiMass
	PRINT*,'Final Mass/Initial Mass', &
	   ABS(REAL(FiMass,kind(0d0))/REAL(InMass,kind(0d0)))
	PRINT*,'Initial energy',InEner
	PRINT*,'Final energy',FiEner
	PRINT*,'Final Energy/Initial Energy', &
	   ABS(REAL(FiEner,kind(0d0))/REAL(InEner,kind(0d0)))

	name_config = 'ufinal.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,Ny
		DO i=1,Nx
			WRITE(11,*) abs(u(i,j))**2
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

	CALL cufftDestroy(plan)

	DEALLOCATE(u,v,temp1,temp2,time,kx,ky,x,y,stat=allocatestatus)	
	IF (allocatestatus .ne. 0) STOP 
	PRINT *,'Deallocated memory'
		
   	PRINT *,'Program execution complete'
	END PROGRAM main
