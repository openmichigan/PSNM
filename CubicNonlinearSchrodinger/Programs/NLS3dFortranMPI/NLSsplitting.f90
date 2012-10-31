




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
	!  dt				= timestep
	!  modescalereal	= Number to scale after backward FFT 
	!  myid				= Process id
	!  ierr				= error code
	!  p_row			= number of rows for domain decomposition
	!  p_col			= number of columns for domain decomposition
	!  filesize			= total filesize
	!  disp				= displacement to start writing data from
	!  ind				= index in array to write
	!  plotnum			= number of plot to save
	!  numberfile		= number of the file to be saved to disk
	!  stat				= error indicator when reading inputfile
	! .. Arrays ..
	!  u 				= approximate solution
	!  v				= Fourier transform of approximate solution
	!  pot				= potential
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  kz				= fourier frequencies in z direction
	!  x				= x locations
	!  y				= y locations
	!  z				= z locations
	!  time				= times at which save data
	!  nameconfig		= array to store filename for data to be saved    		
	!  InputFileName	= name of the Input File
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!					see http://www.2decomp.org/ for more information
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
		
	PROGRAM main
	USE decomp_2d
	USE decomp_2d_fft
	USE decomp_2d_io
	USE MPI
	! Declare variables
	IMPLICIT NONE					 
	INTEGER(kind=4)		::  Nx=2**5
	INTEGER(kind=4) 	::  Ny=2**5
	INTEGER(kind=4) 	::  Nz=2**5
	INTEGER(kind=4)		::  Nt=50	
	INTEGER(kind=4)		::  plotgap=10 
	REAL(kind=8), PARAMETER	::&
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8)		::  Lx=2.0d0,Ly=2.0d0,Lz=2.0d0	 
	REAL(kind=8)		::  Es=1.0d0	 
	REAL(kind=8)		::  dt=0.0010d0		
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE	::  kx,ky,kz
	REAL(kind=8),  	 DIMENSION(:), ALLOCATABLE	::  x,y,z 
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	::  u,v,pot 
	REAL(kind=8), 	 DIMENSION(:), ALLOCATABLE	::  time
	INTEGER(KIND=4), DIMENSION(1:5)	::  intcomm
	REAL(KIND=8), DIMENSION(1:5)	::  dpcomm
	REAL(kind=8) :: modescalereal
	INTEGER(kind=4)	::  i,j,k,n,AllocateStatus,stat
	INTEGER(kind=4)	:: myid,numprocs,ierr
	TYPE(DECOMP_INFO)	::  decomp
	INTEGER(kind=MPI_OFFSET_KIND) :: filesize, disp
	INTEGER(kind=4)	::  p_row=0, p_col=0	
 	INTEGER(kind=4)	::  start, finish, count_rate, ind, plotnum
	CHARACTER*50	::	nameconfig
	CHARACTER*20	::	numberfile, InputFileName
	    ! initialisation of MPI
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
	
	IF(myid.eq.0) THEN	
		CALL GET_ENVIRONMENT_VARIABLE(NAME='inputfile',VALUE=InputFileName, STATUS=stat)
	END IF
	CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	IF(stat.NE.0) THEN
		IF(myid.eq.0) THEN
				PRINT*,"Need to set environment variable inputfile to the name of the &
						file where the simulation parameters are set"
		END IF
		STOP
	END IF
	IF(myid.eq.0) THEN
		InputFileName='./INPUTFILE'
		OPEN(unit=11,FILE=trim(InputFileName),status="OLD") 
		REWIND(11)
		READ(11,*) intcomm(1), intcomm(2), intcomm(3), intcomm(4), intcomm(5), &
			 		dpcomm(1), dpcomm(2),  dpcomm(3),  dpcomm(4),  dpcomm(5)
		CLOSE(11)
		PRINT *,"NX ",intcomm(1)
		PRINT *,"NY ",intcomm(2)
		PRINT *,"NZ ",intcomm(3)
		PRINT *,"NT ",intcomm(4)
		PRINT *,"plotgap ",intcomm(5)
		PRINT *,"Lx ",dpcomm(1)
		PRINT *,"Ly ",dpcomm(2)
		PRINT *,"Lz ",dpcomm(3)
		PRINT *,"Es ",dpcomm(4)		
		PRINT *,"Dt ",dpcomm(5)
		PRINT *,"Read inputfile"
	END IF
	CALL MPI_BCAST(dpcomm,5,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(intcomm,5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	Nx=intcomm(1)
	Ny=intcomm(2)
	Nz=intcomm(3)
	Nt=intcomm(4)
	plotgap=intcomm(5)
	Lx=dpcomm(1)
	Ly=dpcomm(2)
	Lz=dpcomm(3)
	Es=dpcomm(4)
	DT=dpcomm(5)

	! initialisation of 2decomp
	! do automatic domain decomposition
	CALL decomp_2d_init(Nx,Ny,Nz,p_row,p_col)
	! get information about domain decomposition choosen
	CALL decomp_info_init(Nx,Ny,Nz,decomp)
	! initialise FFT library
	CALL decomp_2d_fft_init
	ALLOCATE(u(decomp%xst(1):decomp%xen(1),&
   				decomp%xst(2):decomp%xen(2),&
   				decomp%xst(3):decomp%xen(3)),&
            v(decomp%zst(1):decomp%zen(1),&
            	decomp%zst(2):decomp%zen(2),&
            	decomp%zst(3):decomp%zen(3)),&
            pot(decomp%xst(1):decomp%xen(1),&
            	decomp%xst(2):decomp%xen(2),&
            	decomp%xst(3):decomp%xen(3)),&
   			kx(1:Nx),ky(1:Ny),kz(1:Nz),&
   			x(1:Nx),y(1:Ny),z(1:Nz),&
   			time(1:1+Nt/plotgap),stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 

	IF (myid.eq.0) THEN
		PRINT *,'allocated space'
	END IF	
	
	modescalereal=1.0d0/REAL(Nx,KIND(0d0))
	modescalereal=modescalereal/REAL(Ny,KIND(0d0))
	modescalereal=modescalereal/REAL(Nz,KIND(0d0))
	
	! setup fourier frequencies and grid points
	DO i=1,1+Nx/2
		kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/Lx  			
	END DO
	kx(1+Nx/2)=0.0d0
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO
  	DO i=1,Nx
		x(i)=(-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*Lx
	END DO
	DO j=1,1+Ny/2
		ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly  			
	END DO
	ky(1+Ny/2)=0.0d0
	DO j = 1,Ny/2 -1
		ky(j+1+Ny/2)=-ky(1-j+Ny/2)
	END DO
  	DO j=1,Ny
		y(j)=(-1.0d0 + 2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)))*pi*Ly
	END DO
	DO k=1,1+Nz/2
		kz(k)= cmplx(0.0d0,1.0d0)*REAL(k-1,kind(0d0))/Lz  			
	END DO
	kz(1+Nz/2)=0.0d0
	DO k = 1,Nz/2 -1
		kz(k+1+Nz/2)=-kz(1-k+Nz/2)
	END DO
  	DO k=1,Nz
		z(k)=(-1.0d0 + 2.0d0*REAL(k-1,kind(0d0))/REAL(Nz,kind(0d0)))*pi*Lz
	END DO

	IF (myid.eq.0) THEN
		PRINT *,'Setup grid and fourier frequencies'
	END IF
	
	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				u(i,j,k)=exp(-1.0d0*(x(i)**2 +y(j)**2+z(k)**2))
			END DO
		END DO
	END DO
		
	! write out using 2DECOMP&FFT MPI-IO routines
	nameconfig='./data/u'
	plotnum=0
	WRITE(numberfile,'(i0)') 10000000+plotnum
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//numberfile
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//'.datbin'
	CALL decomp_2d_write_one(1,u,nameconfig)
  
	CALL decomp_2d_fft_3d(u,v,DECOMP_2D_FFT_FORWARD)
	IF (myid.eq.0) THEN
		PRINT *,'Got initial data, starting timestepping'
	END IF
	CALL system_clock(start,count_rate)	
	time(1)=0
	DO n=1,Nt
		! Use Strang splitting
		DO k=decomp%zst(3),decomp%zen(3)
			DO j=decomp%zst(2),decomp%zen(2)
				DO i=decomp%zst(1),decomp%zen(1)
					v(i,j,k)=exp(0.50d0*dt*&
						(kz(k)*kz(k) + kx(i)*kx(i) + ky(j)*ky(j))&
						*cmplx(0.0d0,1.0d0))*v(i,j,k)
				END DO
			END DO
		END DO

		CALL decomp_2d_fft_3d(v,u,DECOMP_2D_FFT_BACKWARD)

		DO k=decomp%xst(3),decomp%xen(3)
			DO j=decomp%xst(2),decomp%xen(2)
				DO i=decomp%xst(1),decomp%xen(1)
					u(i,j,k)=u(i,j,k)*modescalereal
					pot(i,j,k)=Es*u(i,j,k)*conjg(u(i,j,k))
					u(i,j,k)=exp(cmplx(0.0d0,-1.0d0)*dt*pot(i,j,k))*u(i,j,k)
				END DO
			END DO
		END DO	
		CALL decomp_2d_fft_3d(u,v,DECOMP_2D_FFT_FORWARD)

		DO k=decomp%zst(3),decomp%zen(3)
			DO j=decomp%zst(2),decomp%zen(2)
				DO i=decomp%zst(1),decomp%zen(1)
					v(i,j,k)=exp(dt*0.5d0*&
						(kx(i)*kx(i) +ky(j)*ky(j) +kz(k)*kz(k))&
						*cmplx(0.0d0,1.0d0))*v(i,j,k)		
				END DO
			END DO
		END DO	
		IF (mod(n,plotgap)==0) THEN
			time(1+n/plotgap)=n*dt
			IF (myid.eq.0) THEN
				PRINT *,'time',n*dt
			END IF
			CALL decomp_2d_fft_3d(v,u,DECOMP_2D_FFT_BACKWARD)
			u=u*modescalereal
			nameconfig='./data/u'
			plotnum=plotnum+1
			WRITE(numberfile,'(i0)') 10000000+plotnum
			ind=index(nameconfig,' ') -1
			nameconfig=nameconfig(1:ind)//numberfile
			ind=index(nameconfig,' ') -1
			nameconfig=nameconfig(1:ind)//'.datbin'
			! write out using 2DECOMP&FFT MPI-IO routines
			CALL decomp_2d_write_one(1,u,nameconfig)
		END IF
	END DO		
	IF (myid.eq.0) THEN
		PRINT *,'Finished time stepping'
	END IF

	CALL system_clock(finish,count_rate)

	IF (myid.eq.0) THEN
		PRINT*,'Program took ',REAL(finish-start)/REAL(count_rate),'for execution'
	END IF		
		
	IF (myid.eq.0) THEN	
		! Save times at which output was made in text format
		nameconfig = './data/tdata.dat' 
		OPEN(unit=11,FILE=nameconfig,status="UNKNOWN") 	
		REWIND(11)
		DO j=1,1+Nt/plotgap
			WRITE(11,*) time(j)
		END DO
		CLOSE(11)
		! Save x grid points in text format
		nameconfig = './data/xcoord.dat' 
		OPEN(unit=11,FILE=nameconfig,status="UNKNOWN") 	
		REWIND(11)
		DO i=1,Nx
			WRITE(11,*) x(i)
		END DO
		CLOSE(11)
		! Save y grid points in text format
		nameconfig = './data/ycoord.dat' 
		OPEN(unit=11,FILE=nameconfig,status="UNKNOWN") 	
		REWIND(11)
		DO j=1,Ny
			WRITE(11,*) y(j)
		END DO
		CLOSE(11)
		! Save z grid points in text format	
		nameconfig = './data/zcoord.dat' 
		OPEN(unit=11,FILE=nameconfig,status="UNKNOWN") 	
		REWIND(11)
		DO k=1,Nz
			WRITE(11,*) z(k)
		END DO
		CLOSE(11)
		PRINT *,'Saved data'
	END IF
	
	! clean up 
  	CALL decomp_2d_fft_finalize
  	CALL decomp_2d_finalize
 	DEALLOCATE(u,v,pot,&
   			kx,ky,kz,x,y,z,&
   			time,stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 
	IF (myid.eq.0) THEN
	   	PRINT *,'Program execution complete'
	END IF
	CALL MPI_FINALIZE(ierr)		
	END PROGRAM main