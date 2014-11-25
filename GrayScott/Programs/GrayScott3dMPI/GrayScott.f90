    !--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves GrayScott equation in 3 dimensions
	! u_t=D_u*(u_{xx}+u_{yy}+u_{zz})^2-u*v^2-F*(1-u)
	! v_t=D_v*(v_{xx}+v_{yy}+v_{zz})^2+u*v^2-(F+k)*v
	! 
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
	!  A				=F
	!  B				=F+k
	! .. Scalars ..
	!  Ahat 			= Fourier transform of A
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
	!  umean			= temporary for fix point iteration
	!  vmean			= temporary for fix point iteration
	! .. Arrays ..
	!  u 				= approximate solution u
	!  v				= approximate solution v
	!  uold 				= previous timestep of u
	!  vold				= previous timestep of v
	!  utemp 			= temporary for fix point iteration u
	!  vtemp			= temporary for fix point iteration v
	!  uhat 			= Fourier transform of u
	!  vhat				= Fourier transform of v
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction squared
	!  ky				= fourier frequencies in y direction squared
	!  kz				= fourier frequencies in z direction squared
	!  x				= x locations
	!  y				= y locations
	!  z				= z locations
	!  time				= times at which save data
	!  nameconfig		= array to store filename for data to be saved    		
	!  InputFileName	= name of the Input File
	!  dpcomm
	!  intcomm
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!  sp				see http://www.2decomp.org/ for more information
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
	INTEGER(kind=4)		::  Nx=0
	INTEGER(kind=4) 	::  Ny=0
	INTEGER(kind=4) 	::  Nz=0
	INTEGER(kind=4)		::  Nt=0	
	INTEGER(kind=4)		::  plotgap=0 
	REAL(kind=8), PARAMETER	::&
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8)		::  Lx=0.0,Ly=0.0,Lz=0.0
	REAL(kind=8)		::  dt=0.0		
	REAL(kind=8)		::  tol=0.1**12
	REAL(kind=8)		::  A=0.0
	REAL(kind=8)		::  Ahat=0
	REAL(kind=8)		::  B=0.0
	REAL(kind=8)		::  Du=0.0
	REAL(kind=8)		::  Dv=0.0
	REAL(kind=8)		::  chg=1.0
	REAL(kind=8), DIMENSION(:), ALLOCATABLE	::  kx,ky,kz
	REAL(kind=8),  	 DIMENSION(:), ALLOCATABLE	::  x,y,z 
	REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE	::  u,v
	REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE	::  uold,vold
	REAL(kind=8) 		:: utemp,vtemp
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	::  uhat,vhat
	COMPLEX(kind=8) :: uhatmp
	REAL(kind=8) :: umean,vmean
	REAL(kind=8), 	 DIMENSION(:), ALLOCATABLE	::  time
	INTEGER(KIND=4), DIMENSION(1:5)	::  intcomm
	REAL(KIND=8), DIMENSION(1:8)	::  dpcomm
	REAL(kind=8) :: modescalereal
	INTEGER(kind=4)	::  i,j,k,n,AllocateStatus,stat
	INTEGER(kind=4)	:: myid,numprocs,ierr
	TYPE(DECOMP_INFO)	::  decomp,sp
	INTEGER(kind=MPI_OFFSET_KIND) :: filesize, disp
	INTEGER(kind=4)	::  p_row=0, p_col=0	
 	INTEGER(kind=4)	::  start, finish, count_rate, ind, plotnum
	CHARACTER*50	::	nameconfig
	CHARACTER*20	::	numberfile, InputFileName
	    ! initialisation of MPI
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
	

	CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	IF(myid.eq.0) THEN
		InputFileName='INPUTFILE'
		OPEN(unit=11,FILE=trim(InputFileName),status="OLD") 
		REWIND(11)
		READ(11,*) intcomm(1), intcomm(2), intcomm(3), intcomm(4),intcomm(5),&
			 		dpcomm(1), dpcomm(2),  dpcomm(3),  dpcomm(4),  dpcomm(5), &
					dpcomm(6), dpcomm(7),  dpcomm(8)
		CLOSE(11)
		PRINT *,"NX ",intcomm(1)
		PRINT *,"NY ",intcomm(2)
		PRINT *,"NZ ",intcomm(3)
		PRINT *,"NT ",intcomm(4)
		PRINT *,"plotgap ",intcomm(5)
		PRINT *,"Lx",dpcomm(1)
		PRINT *,"Ly",dpcomm(2)
		PRINT *,"Lz",dpcomm(3)
		PRINT *,"dt",dpcomm(4)		
		PRINT *,"Du",dpcomm(5)
		PRINT *,"Dv",dpcomm(6)
		PRINT *,"F",dpcomm(7)
		PRINT *,"k",dpcomm(8)		
		PRINT *,"Read inputfile"
	END IF
	CALL MPI_BCAST(dpcomm,8,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(intcomm,5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	Nx=intcomm(1)
	Ny=intcomm(2)
	Nz=intcomm(3)
	Nt=intcomm(4)
	plotgap=intcomm(5)
	Lx=dpcomm(1)
	Ly=dpcomm(2)
	Lz=dpcomm(3)
	dt=dpcomm(4)
	Du=dpcomm(5)
	Dv=dpcomm(6)
	A=dpcomm(7)
	B=dpcomm(8)+dpcomm(7)
	! initialisation of 2decomp
	! do automatic domain decomposition
	CALL decomp_2d_init(Nx,Ny,Nz,p_row,p_col)
	! get information about domain decomposition choosen
	CALL get_decomp_info(decomp) ! physical domain
	CALL decomp_info_init(Nx/2+1,Ny,Nz,sp) ! spectral domain
	! initialise FFT library
	CALL decomp_2d_fft_init
	ALLOCATE(u(decomp%xst(1):decomp%xen(1),&
   				decomp%xst(2):decomp%xen(2),&
   				decomp%xst(3):decomp%xen(3)),&
            v(decomp%xst(1):decomp%xen(1),&
            	decomp%xst(2):decomp%xen(2),&
            	decomp%xst(3):decomp%xen(3)),&
			uold(decomp%xst(1):decomp%xen(1),&
   				decomp%xst(2):decomp%xen(2),&
   				decomp%xst(3):decomp%xen(3)),&
            vold(decomp%xst(1):decomp%xen(1),&
            	decomp%xst(2):decomp%xen(2),&
            	decomp%xst(3):decomp%xen(3)),&
			uhat(sp%zst(1):sp%zen(1),&
   				sp%zst(2):sp%zen(2),&
   				sp%zst(3):sp%zen(3)),&
            vhat(sp%zst(1):sp%zen(1),&
            	sp%zst(2):sp%zen(2),&
            	sp%zst(3):sp%zen(3)),&
		kx(sp%zst(1):sp%zen(1)), &
		ky(sp%zst(2):sp%zen(2)), &
		kz(sp%zst(3):sp%zen(3)), &
		x(decomp%xst(1):decomp%xen(1)),&
		y(decomp%xst(2):decomp%xen(2)),&
		z(decomp%xst(3):decomp%xen(3)),&
   		time(1:1+Nt/plotgap),&
		stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 

	IF (myid.eq.0) THEN
		PRINT *,'allocated space'
	END IF	
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	modescalereal=1.0d0/REAL(Nx,KIND(0d0))
	modescalereal=modescalereal/REAL(Ny,KIND(0d0))
	modescalereal=modescalereal/REAL(Nz,KIND(0d0))
	
	! setup fourier frequencies and grid points
DO i = 1,1+ Nx/2
		IF ((i.GE.sp%zst(1)).AND.(i.LE.sp%zen(1))) THEN
			kx(i)=-(REAL(i-1,kind(0d0))/Lx)**2
		END IF
	END DO
	IF ((Nx/2 + 1 .GE.sp%zst(1)).AND.(Nx/2 + 1 .LE.sp%zen(1))) THEN
		kx( Nx/2 + 1 ) = 0.0d0 
	END IF
	DO i = Nx/2+2, Nx  
		IF ((i.GE.sp%zst(1)).AND.(i.LE.sp%zen(1))) THEN
			Kx( i) = -(REAL(1-i+Nx,KIND(0d0))/Lx)**2  
		END IF
	END DO      
	DO i=decomp%xst(1),decomp%xen(1)
		x(i)=(-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*Lx
	END DO

	DO j = 1,1+ Ny/2
		IF ((j.GE.sp%zst(2)).AND.(j.LE.sp%zen(2))) THEN
			ky(j)= -(REAL(j-1,kind(0d0))/Ly)**2
		END IF
	END DO
	IF ((Ny/2 + 1 .GE.sp%zst(2)).AND.(Ny/2 + 1 .LE.sp%zen(2))) THEN
		ky( Ny/2 + 1 ) = 0.0d0 
	ENDIF
	DO j = Ny/2+2, Ny  
		IF ((j.GE.sp%zst(2)).AND.(j.LE.sp%zen(2))) THEN
			ky(j) = -(REAL(1-j+Ny,KIND(0d0))/Ly)**2  
		END IF
	END DO      
	DO j=decomp%xst(2),decomp%xen(2)
		y(j)=(-1.0d0 + 2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)))*pi*Ly
	END DO
	
	DO k = 1,1+ Nz/2
		IF ((k.GE.sp%zst(3)).AND.(k.LE.sp%zen(3))) THEN
			kz(k)= -(REAL(k-1,kind(0d0))/Lz)**2
		END IF
	END DO
	IF ((Nz/2 + 1 .GE.sp%zst(3)).AND.(Nz/2 + 1 .LE.sp%zen(3))) THEN
		kz( Nz/2 + 1 ) = 0.0d0 
	ENDIF
	DO k = Nz/2+2, Nz  
		IF ((k.GE.sp%zst(3)).AND.(k.LE.sp%zen(3))) THEN
			kz(k) = -(REAL(1-k+Nz,KIND(0d0))/Lz)**2  
		END IF
	END DO      
	DO k=decomp%xst(3),decomp%xen(3)
		z(k)=(-1.0d0 + 2.0d0*REAL(k-1,kind(0d0))/REAL(Nz,kind(0d0)))*pi*Lz
	END DO

	IF (myid.eq.0) THEN
		PRINT *,'Setup grid and fourier frequencies'
	END IF
	!compute Ahat
	Ahat=A*Nx*Ny*Nz
	!initial condition
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				u(i,j,k)=1+(exp(-1*(x(i)**2+y(j)**2+z(k)**2)-1))			
				v(i,j,k)=0+(exp(-1*(x(i)**2+y(j)**2+z(k)**2)-1))
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
	!same for v
	nameconfig='./data/v'
	plotnum=0
	WRITE(numberfile,'(i0)') 10000000+plotnum
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//numberfile
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//'.datbin'
	CALL decomp_2d_write_one(1,v,nameconfig)
	
	IF (myid.eq.0) THEN
		PRINT *,'Got initial data, starting timestepping'
	END IF
	CALL system_clock(start,count_rate)	
	time(1)=0
	DO n=1,Nt
	!use fixed point iteration to solve nonlinear system in real space
		DO k=decomp%xst(3),decomp%xen(3)
			DO j=decomp%xst(2),decomp%xen(2)
				DO i=decomp%xst(1),decomp%xen(1)
					uold(i,j,k)=u(i,j,k)
					vold(i,j,k)=v(i,j,k)
				END DO
			END DO
		END DO
		DO k=decomp%xst(3),decomp%xen(3)
			DO j=decomp%xst(2),decomp%xen(2)
				DO i=decomp%xst(1),decomp%xen(1)
					chg=1
					DO WHILE (chg>tol)
						utemp=u(i,j,k)
						vtemp=v(i,j,k)
						umean=0.5*(u(i,j,k)+uold(i,j,k))
						vmean=0.5*(v(i,j,k)+vold(i,j,k))
						u(i,j,k)=uold(i,j,k)+0.5*dt*(-umean*vmean**2)
						v(i,j,k)=vold(i,j,k)+0.5*dt*(umean*vmean**2)
						utemp=u(i,j,k)-utemp
						vtemp=v(i,j,k)-vtemp
						chg=abs(utemp)+abs(vtemp)
				END DO
			END DO
		END DO			
	END DO
! solve linear part exactly in Fourier space
		CALL decomp_2d_fft_3d(u,uhat)
		CALL decomp_2d_fft_3d(v,vhat)
		IF (myid.eq.0) THEN
			uhatmp=uhat(1,1,1)
		END IF
		DO k=sp%zst(3),sp%zen(3)
			DO j=sp%zst(2),sp%zen(2)
				DO i=sp%zst(1),sp%zen(1)
					uhat(i,j,k)=exp(dt*(-A+Du*(kx(i)+ky(j)+kz(k))))*&
								uhat(i,j,k)
					vhat(i,j,k)=exp(dt*(-B+Dv*(kx(i)+ky(j)+kz(k))))*&
								vhat(i,j,k)
				END DO
			END DO
		END DO
		IF (myid.eq.0) THEN
			uhat(1,1,1)=exp(dt*(-A+Du*(kx(1)+ky(1)+kz(1))))*&
						(uhatmp-Ahat/(A+Du*(kx(1)+ky(1)+kz(1))))+&
						(Ahat/(A+Du*(kx(1)+ky(1)+kz(1))))
		END IF
		!dont forget to scale
		CALL decomp_2d_fft_3d(uhat,u)
		CALL decomp_2d_fft_3d(vhat,v)
		DO k=decomp%xst(3),decomp%xen(3)
			DO j=decomp%xst(2),decomp%xen(2)
				DO i=decomp%xst(1),decomp%xen(1)
					u(i,j,k)=u(i,j,k)*modescalereal
					v(i,j,k)=v(i,j,k)*modescalereal
				END DO
			END DO
		END DO
		!use fixed point iteration to solve nonlinear system in real space
		DO k=decomp%xst(3),decomp%xen(3)
			DO j=decomp%xst(2),decomp%xen(2)
				DO i=decomp%xst(1),decomp%xen(1)
					uold(i,j,k)=u(i,j,k)
					vold(i,j,k)=v(i,j,k)
				END DO
			END DO
		END DO
		DO k=decomp%xst(3),decomp%xen(3)
			DO j=decomp%xst(2),decomp%xen(2)
				DO i=decomp%xst(1),decomp%xen(1)
					chg=1
					DO WHILE (chg>tol)
						utemp=u(i,j,k)
						vtemp=v(i,j,k)
						umean=0.5*(u(i,j,k)+uold(i,j,k))
						vmean=0.5*(v(i,j,k)+vold(i,j,k))
						u(i,j,k)=uold(i,j,k)+0.5*dt*(-umean*vmean**2)
						v(i,j,k)=vold(i,j,k)+0.5*dt*(umean*vmean**2)
						utemp=u(i,j,k)-utemp
						vtemp=v(i,j,k)-vtemp
						chg=abs(utemp)+abs(vtemp)
					END DO
				END DO
			END DO
		END DO			
!write data		
		IF (mod(n,plotgap)==0) THEN
			time(1+n/plotgap)=n*dt
			IF (myid.eq.0) THEN
				PRINT *,'time',n*dt
			END IF
			nameconfig='./data/u'
			plotnum=plotnum+1
			WRITE(numberfile,'(i0)') 10000000+plotnum
			ind=index(nameconfig,' ') -1
			nameconfig=nameconfig(1:ind)//numberfile
			ind=index(nameconfig,' ') -1
			nameconfig=nameconfig(1:ind)//'.datbin'
			!write out using 2DECOMP&FFT MPI-IO routines
			CALL decomp_2d_write_one(1,u,nameconfig)
			nameconfig='./data/v'
			WRITE(numberfile,'(i0)') 10000000+plotnum
			ind=index(nameconfig,' ') -1
			nameconfig=nameconfig(1:ind)//numberfile
			ind=index(nameconfig,' ') -1
			nameconfig=nameconfig(1:ind)//'.datbin'
			!write out using 2DECOMP&FFT MPI-IO routines
			CALL decomp_2d_write_one(1,v,nameconfig)
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
 	DEALLOCATE(u,v,uold,vold,uhat,vhat,&
   		kx,ky,kz,x,y,z,&
   		time,stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 
	IF (myid.eq.0) THEN
	   	PRINT *,'Program execution complete'
	END IF
	CALL MPI_FINALIZE(ierr)		
	END PROGRAM main
