!-----------------------------------------------------------------------------------
!
!
! PURPOSE
!
! This program numerically solves the 3D incompressible Navier-Stokes
! on a Cubic Domain [0,2pi]x[0,2pi]x[0,2pi] using pseudo-spectral methods and
! fourth order Carpenter Kennedy timestepping. The numerical solution is compared to
! an exact solution reported by Shapiro
!
! Analytical Solution:
! u(x,y,z,t)=-0.25*(cos(x)sin(y)sin(z)+sin(x)cos(y)cos(z))exp(-t/Re)
! v(x,y,z,t)= 0.25*(sin(x)cos(y)sin(z)-cos(x)sin(y)cos(z))exp(-t/Re)
! w(x,y,z,t)= 0.5*cos(x)cos(y)sin(z)exp(-t/Re)
!
! .. Parameters ..
! Nx = number of modes in x - power of 2 for FFT
! Ny = number of modes in y - power of 2 for FFT
! Nz = number of modes in z - power of 2 for FFT
! Nt = number of timesteps to take
! Tmax = maximum simulation time
! pi = 3.14159265358979323846264338327950288419716939937510d0
! Re = Reynolds number
! .. Scalars ..
! i = loop counter in x direction
! j = loop counter in y direction
! k = loop counter in z direction
! n = loop counter for timesteps direction
! nn = loop conuter for timesteps between plots
! numplots = number of plots to make
! allocatestatus = error indicator during allocation
! count = keep track of information written to disk
! iol = size of array to write to disk
! start = variable to record start time of program
! finish = variable to record end time of program
! count_rate = variable for clock count rate
! theta = variable to modify exact solution used
! dt = timestep
! filesize = MPI variable to indicate file size to be written to disk
! disp = MPI variable to indicate where in file local portion to be written
! p_row = variable used by 2decomp to store number of rows
! p_col= variable used by 2decomp to store number of columns
! numprocs = number of processros
! myid = number of local process
! ierr = error codes returned in this variable
! t = time
! plotInt= plot number
! ind = index for file number
! start = timing variable
! finish = fiming variable
! count_rate = clock rate to convert to seconds

! .. Arrays ..
! u = velocity in x direction
! v = velocity in y direction
! w = velocity in z direction
! uexact = array for exact solution
! vexact = array for exact solution
! wexact = array for exact solution
! uold = velocity in x direction at previous timestep
! vold = velocity in y direction at previous timestep
! wold = velocity in z direction at previous timestep
! ux = x derivative of velocity in x direction
! uy = y derivative of velocity in x direction
! uz = z derivative of velocity in x direction
! vx = x derivative of velocity in y direction
! vy = y derivative of velocity in y direction
! vz = z derivative of velocity in y direction
! wx = x derivative of velocity in z direction
! wy = y derivative of velocity in z direction
! wz = z derivative of velocity in z direction
! uhat = Fourier transform of u
! vhat = Fourier transform of v
! what = Fourier transform of w
! nonlinuhat = Fourier transform of nonlinear term for u
! nonlinvhat = Fourier transform of nonlinear term for u
! nonlinwhat = Fourier transform of nonlinear term for u
! phat = Fourier transform of nonlinear term for pressure, p
! temp_c = temporary storage for Fourier transforms
!
! .. Vectors ..
! kx = fourier frequencies in x direction
! ky = fourier frequencies in y direction
! kz = fourier frequencies in z direction
! x = x locations
! y = y locations
! z = y locations
! time = times at which save data
! name_config = array to store filename for data to be saved
! number_file = array to store number of filename being saved
! temp1 = array to collect local sums for MPI reductions
! temp2 = array to perform MPI reductions in
! mychg = array to collect local differences between fixed point iterates
! allchg = array to collect global differences between fixed point iterates
! alpha = array of Runge Kutta Coefficients
! beta = array of Runge Kutta Coefficients
! gamma = array of Runge Kutta Coefficients
! KE = array to record kinetic energy
! Enstrophy = array to record enstrophy
! KEdissipationRate = array to record kinetic energy dissipation rate
!
! .. Special types ..
!
! decomp = structure from 2decomp with information on complex array sizes
! sp = structure from 2decomp with information on real array sizes
!
! REFERENCES
!
! A. Shapiro " The use of an exact solution of the Navier-Stokes equations
! in a validation test of a three-dimensional nonhydrostatic numerical model"
! Monthly Weather Review vol. 121, 2420-2425, (1993).
!
! ACKNOWLEDGEMENTS
!
! ACCURACY
!
! ERROR INDICATORS AND WARNINGS
!
! FURTHER COMMENTS
!
! This program has not been optimized to use the least amount of memory
! or be as efficeint as possible, but is intended to be portable and produc
! correcte results
!
!--------------------------------------------------------------------------------
! External routines required
!
! External libraries required
! 2DECOMP&FFT -- Fast Fourier Transform in the West Library
! (http://2decomp.org/)

PROGRAM main
	USE decomp_2d
	USE decomp_2d_fft
	USE decomp_2d_io
	USE MPI	
	IMPLICIT NONE	
	!Grid & Stepping
   	INTEGER(kind=4), PARAMETER              :: Nx=128, Ny=128, Nz=128
   	INTEGER(kind=4), PARAMETER              :: Lx=1, Ly=1, Lz=1
	INTEGER(kind=4), PARAMETER              :: Nt=20,plotgap=10
	INTEGER					:: numplots=Nt/plotgap
	!choose initial conditions
	logical                                 ::readinput=.false.,&
                                            TGvortex=.false.,&
                                            exactsoln=.true.
	!runime parameters
	INTEGER(kind=4), PARAMETER              :: savedata1=0
	REAL(kind=8), PARAMETER                 :: dt=0.0001d0
	REAL(kind=8), PARAMETER                 :: Re=1600.0d0
	REAL(kind=8)                            :: theta=0.0d0
	!additional parameters & constants
	REAL(kind=8), PARAMETER	&
		::  pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER                 :: dx=2.0d0*pi*Lx/REAL(Nx,kind(0d0))
	REAL(kind=8), PARAMETER                 :: dy=2.0d0*pi*Ly/REAL(Ny,kind(0d0))
	REAL(kind=8), PARAMETER                 :: dz=2.0d0*pi*Lz/REAL(Nz,kind(0d0))
	REAL(kind=8), PARAMETER                 :: ReInv=1.0d0/REAL(Re,kind(0d0))
	REAL(kind=8), PARAMETER                 :: dtInv=1.0d0/REAL(dt,kind(0d0))
	!Carpenter-Kennedy coefficients
	REAL(kind=8), DIMENSION(6), PARAMETER	:: alpha=(/0.0,0.1496590219993,0.3704009573644,0.6222557631345,0.9582821306748,1.0/)
	REAL(kind=8), DIMENSION(5), PARAMETER	:: beta=(/0.0,-0.4178904745,-1.192151694643,-1.697784692471,-1.514183444257/)
	REAL(kind=8), DIMENSION(5), PARAMETER	:: gamma1=(/0.1496590219993,0.3792103129999,0.8229550293869,0.6994504559488,0.1530572479681/)		
	!computational arrays/scalars	
	REAL(kind=8)                                    ::  scalemodes,mu,factor,t=0.0d0
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE      :: kx,ky,kz
	REAL(kind=8), DIMENSION(:), ALLOCATABLE         :: x,y,z,time
	REAL(kind=8), DIMENSION(:), ALLOCATABLE         :: KE,Enstrophy,KEdissipationRate
	REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE     :: u,v,w,&
                                                    ux,uy,uz,&
                                                    vx,vy,vz,&
                                                    wx,wy,wz,&
                                                    omega,nonlin,temp
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	:: uhat,vhat,what,&
                                                    omegahat,temp_c,temp_c1,&
                                                    nonlinhatuh,nonlinhatvh,&
                                                    nonlinhatwh,phat
	!exact solution
	REAL(kind=8), PARAMETER                         :: sl=1.0d0,sk=1.0d0,sm=1.0d0,&
                                                    lamlkm=sqrt(sl**2+sk**2+sm**2)
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	:: uexact,vexact,wexact
	REAL(kind=8)                                    :: temp3
	!MPI and 2DECOMP variables 
	TYPE(DECOMP_INFO)                               :: decomp,sp
	INTEGER(kind=4)                                 :: myid,numprocs,ierr
	INTEGER(kind=MPI_OFFSET_KIND)                   :: filesize, disp
	INTEGER(kind=4)                                 :: p_row=0, p_col=0
	REAL(kind=8), DIMENSION(:), ALLOCATABLE         :: temp1,temp2
	!variables used for saving data and timing
	INTEGER(kind=4)                                 ::  i,j,k,n,nn,kk,allocatestatus
	INTEGER(kind=4)                                 ::  count, iol, ind, plotInt=0
	CHARACTER*100                                   :: name_config
	INTEGER(kind=4)                                 :: start, finish, count_rate
   !initialisation of 2DECOMP&FFT and MPI	
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
	CALL decomp_2d_init(Nx,Ny,Nz,p_row,p_col)
	CALL decomp_info_init(Nx,Ny,Nz,decomp)
   	CALL decomp_info_init(Nx/2+1,Ny,Nz,sp)
	CALL decomp_2d_fft_init
	
	IF (myid.eq.0) THEN
		PRINT *,'Solve 3D incompressible NS equations with Fourier pseudospectral methods'
		PRINT *,'and Carpenter-Kennedy timestepping.'
	   	PRINT *,'Grid:',Nx,'X',Ny,'X',Nz
		PRINT *,'dt:',dt
		PRINT *,'Re:',Re
		PRINT *,'Nt:',Nt
		PRINT *,'plotgap:',plotgap
	END IF		
	ALLOCATE(time(1:numplots+1),kx(1:Nx),ky(1:Ny),kz(1:Nz),x(1:Nx),y(1:Ny),z(1:Nz),temp1(1:9),temp2(1:9),&
	KE(1:numplots+1),Enstrophy(1:numplots+1),KEdissipationRate(1:numplots+1),&
	u(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	v(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	w(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	ux(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	uy(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	uz(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	vx(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	vy(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	vz(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	wx(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	wy(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	wz(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	omega(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	nonlin(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	temp(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	uhat(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	vhat(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	what(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	omegahat(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	temp_c(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	temp_c1(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	nonlinhatuh(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	nonlinhatvh(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	nonlinhatwh(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	phat(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	uexact(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	vexact(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	wexact(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	stat=allocatestatus)
	IF (allocatestatus .ne. 0) STOP 
	IF (myid.eq.0) THEN
	 	PRINT *,'allocated space'
	END IF	
	
	! setup fourier frequencies in x-direction
	DO i=1,Nx/2+1
		kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/Lx  			
	END DO
	kx(1+Nx/2)=0.0d0
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO	
	ind=1
	DO i=-Nx/2,Nx/2-1
		x(ind)=2.0d0*pi*REAL(i,kind(0d0))*Lx/REAL(Nx,kind(0d0))
		ind=ind+1
	END DO
	! setup fourier frequencies in y-direction
	DO j=1,Ny/2+1
		ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly  			
	END DO
	ky(1+Ny/2)=0.0d0
	DO j = 1,Ny/2 -1
		ky(j+1+Ny/2)=-ky(1-j+Ny/2)
	END DO	
	ind=1
	DO j=-Ny/2,Ny/2-1
		y(ind)=2.0d0*pi*REAL(j,kind(0d0))*Ly/REAL(Ny,kind(0d0))
		ind=ind+1
	END DO
	! setup fourier frequencies in z-direction
	DO k=1,Nz/2+1
		kz(k)= cmplx(0.0d0,1.0d0)*REAL(k-1,kind(0d0))/Lz  			
	END DO
	kz(1+Nz/2)=0.0d0
	DO k = 1,Nz/2 -1
		kz(k+1+Nz/2)=-kz(1-k+Nz/2)
	END DO	
	ind=1
	DO k=-Nz/2,Nz/2-1
		z(ind)=2.0d0*pi*REAL(k,kind(0d0))*Lz/REAL(Nz,kind(0d0))
		ind=ind+1
	END DO
	IF (myid.eq.0) THEN
		PRINT *,'Setup grid and fourier frequencies'
	END IF	
	scalemodes=1.0d0/REAL(Nx*Ny*Nz,kind(0d0)); time(1)=0.0d0
	
   IF(readInput.eqv..true.) THEN
		name_config='./data/restartData'
		OPEN(unit=11,file=name_config,status='old')
		REWIND(11)
		READ(11,*) time(1),plotInt
		CLOSE(11)
      name_config='./data/u'
      CALL getName(name_config,plotInt)
      CALL decomp_2d_read_one(1,u,name_config)
      name_config='./data/v'
      CALL getName(name_config,plotInt)
      CALL decomp_2d_read_one(1,v,name_config)
      name_config='./data/w'
      CALL getName(name_config,plotInt)
      CALL decomp_2d_read_one(1,w,name_config)
      IF(myid.eq.0) THEN
         PRINT*,'Read in values'
      END IF
   ELSE If(TGvortex.eqv..true.) THEN
	   !initial conditions for Taylor-Green vortex
	   factor=2.0d0/sqrt(3.0d0)
      DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	      u(i,j,k)=factor*sin(theta+2.0d0*pi/3.0d0)*sin(x(i))*cos(y(j))*cos(z(k))
      END DO; END DO; END DO
      DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	      v(i,j,k)=factor*sin(theta-2.0d0*pi/3.0d0)*cos(x(i))*sin(y(j))*cos(z(k))
      END DO ; END DO ; END DO
      DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	      w(i,j,k)=factor*sin(theta)*cos(x(i))*cos(y(j))*sin(z(k))
      END DO ; END DO ; END DO
      IF(myid.eq.0) THEN
         PRINT*,'Taylor-Green Vortex initial conditions'
      END IF
   ELSE IF(exactsoln.eqv..true.) THEN
		!special exact sol'n usted for testing
		factor=sqrt(3.0d0)
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			u(i,j,k)=-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
							+sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
		END DO; END DO; END DO
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			v(i,j,k)=0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
							-cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
		END DO ; END DO ; END DO
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			w(i,j,k)=cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(1)/Re)
		END DO ; END DO ; END DO
      IF(myid.eq.0) THEN
         PRINT*,'Exact Solution initial conditions'
      END IF
	ELSE
      IF(myid.eq.0) THEN
         PRINT*,'--------------------------------------------------------------------'
		   PRINT*,'Set initial conditions to either readinput, TGvortex, or exactsoln.'
         PRINT*,'--------------------------------------------------------------------'
         goto 111
      END IF
	END IF
	
	CALL decomp_2d_fft_3d(u,uhat)
	CALL decomp_2d_fft_3d(v,vhat)
	CALL decomp_2d_fft_3d(w,what)
	
	! derivative of u with respect to x, y, and z 
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,ux)	
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,uy)	
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,uz)	
    
	! derivative of v with respect to x, y, and z 
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,vx)		
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,vy)	
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,vz)		

	! derivative of w with respect to x, y, and z 
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wx)		
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wy)		
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wz)	
  
   !calculate vorticity
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		temp(i,j,k)=wy(i,j,k)-vz(i,j,k)
	END DO ; END DO ; END DO

	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		omega(i,j,k)=temp(i,j,k)*temp(i,j,k)
	END DO ; END DO ; END DO

	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		temp(i,j,k)=uz(i,j,k)-wx(i,j,k)
	END DO ; END DO ; END DO

	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		omega(i,j,k)=temp(i,j,k)*temp(i,j,k)+omega(i,j,k)
	END DO ; END DO ; END DO

	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		temp(i,j,k)=vx(i,j,k)-uy(i,j,k)
	END DO ; END DO ; END DO
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		omega(i,j,k)=temp(i,j,k)*temp(i,j,k)+omega(i,j,k)
	END DO ; END DO ; END DO

	!Calculate Kinetic Energy
	temp1(1) = sum(u*u)
	temp1(2) = sum(v*v)
	temp1(3) = sum(w*w)
	CALL MPI_ALLREDUCE(temp1(1:3),temp2(1:3),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	temp3=temp2(1)+temp2(2)+temp2(3)
	KE(1)=dx*dy*dz*temp3/(8.0d0*2.0d0*pi*pi*pi)

	!Calculate Enstrophy
   temp1(1)=sum(omega)
	CALL MPI_ALLREDUCE(temp1(1),temp3,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   Enstrophy(1)=dx*dy*dz*temp3/(8.0d0*2.0d0*pi*pi*pi)
	
   temp1=0.0;temp2=0.0
	!Calculate Kinetic Energy Dissipation Rate
	temp1(1)=sum(ux*ux)
	temp1(2)=sum(uy*uy)
	temp1(3)=sum(uz*uz)
	temp1(4)=sum(vx*vx)
	temp1(5)=sum(vy*vy)
	temp1(6)=sum(vz*vz)
	temp1(7)=sum(wx*wx)
	temp1(8)=sum(wy*wy)
	temp1(9)=sum(wz*wz)
	CALL MPI_ALLREDUCE(temp1,temp2,9,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	temp3=temp2(1)+temp2(2)+temp2(3)+temp2(4)+temp2(5)+temp2(6)+temp2(7)+temp2(8)+temp2(9)
	KEdissipationRate(1)=dx*dy*dz*temp3/(8.0d0*Re*pi*pi*pi)

   IF(myid.eq.0) THEN
      PRINT *,'time=',time(1)
	   PRINT*,'KE=',KE(1)
	   PRINT*,'Enstrophy=',Enstrophy(1)
	   PRINT*,'KEdissipationRate=',KEdissipationRate(1)
      PRINT*,'---------------------------------'
		PRINT*,'starting timestepping and timer'
     	PRINT*,'---------------------------------'
   END IF 
	nonlinhatuh(:,:,:)=0.0; nonlinhatvh(:,:,:)=0.0; nonlinhatwh(:,:,:)=0.0;
   phat(:,:,:)=0.0

	CALL system_clock(start,count_rate)   
	DO n=1,numplots
      DO nn=1,plotgap
         t=t+dt
		   DO kk=1,5
			   !nonlinuhat
			   DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   temp(i,j,k)=(u(i,j,k)*ux(i,j,k)+v(i,j,k)*uy(i,j,k)+w(i,j,k)*uz(i,j,k))
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp,temp_c1)
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
			    	temp_c(i,j,k)=kx(i)*temp_c1(i,j,k)
			   END DO ; END DO ; END DO
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   nonlinhatuh(i,j,k)=temp_c1(i,j,k)+REAL(beta(kk))*nonlinhatuh(i,j,k)
			   END DO ; END DO ; END DO
				
			   !nonlinvhat
			   DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   temp(i,j,k)=u(i,j,k)*vx(i,j,k)+v(i,j,k)*vy(i,j,k)+w(i,j,k)*vz(i,j,k)
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp,temp_c1)
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
			    	temp_c(i,j,k)=ky(j)*temp_c1(i,j,k)+temp_c(i,j,k)
			   END DO ; END DO ; END DO
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   nonlinhatvh(i,j,k)=temp_c1(i,j,k)+REAL(beta(kk))*nonlinhatvh(i,j,k)
			   END DO ; END DO ; END DO
			    
			   !nonlinwhat
			   DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   temp(i,j,k)=u(i,j,k)*wx(i,j,k)+v(i,j,k)*wy(i,j,k)+w(i,j,k)*wz(i,j,k)
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp,temp_c1)
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
			    	temp_c(i,j,k)=kz(k)*temp_c1(i,j,k)+temp_c(i,j,k)
			   END DO ; END DO ; END DO
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   nonlinhatwh(i,j,k)=temp_c1(i,j,k)+REAL(beta(kk))*nonlinhatwh(i,j,k)
			   END DO ; END DO ; END DO

			   !phat
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   phat(i,j,k)=-1.0*temp_c(i,j,k)/(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)&
								   + REAL(0.1d0**13,kind(0d0)))+REAL(beta(kk))*phat(i,j,k)
			   END DO ; END DO ; END DO	
			   mu=0.5d0*dt*(alpha(kk+1)-alpha(kk))

			   !uhat,vhat,what
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   factor=(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)+REAL(0.1d0**13,kind(0d0)))
				   uhat(i,j,k)=((uhat(i,j,k)-gamma1(kk)*dt*(nonlinhatuh(i,j,k)+kx(i)*phat(i,j,k))&
					   +mu*ReInv*factor*uhat(i,j,k))/(1.0d0 - mu*ReInv*factor));
				   vhat(i,j,k)=((vhat(i,j,k)-gamma1(kk)*dt*(nonlinhatvh(i,j,k)+ky(j)*phat(i,j,k))&
					   +mu*ReInv*factor*vhat(i,j,k))/(1.0d0 - mu*ReInv*factor));
				   what(i,j,k)=((what(i,j,k)-gamma1(kk)*dt*(nonlinhatwh(i,j,k)+kz(k)*phat(i,j,k))&
					   +mu*ReInv*factor*what(i,j,k))/(1.0d0 - mu*ReInv*factor));
			   END DO ; END DO ; END DO
			
			   !ifft to get u,v,w and derivatives
			
			
			   ! derivative of u with respect to x, y, and z 
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp_c,ux)	
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp_c,uy)	
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp_c,uz)	

			   ! derivative of v with respect to x, y, and z 
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp_c,vx)		
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp_c,vy)	
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp_c,vz)		

			   ! derivative of w with respect to x, y, and z 
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp_c,wx)		
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp_c,wy)		
			   DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
			   END DO ; END DO ; END DO
			   CALL decomp_2d_fft_3d(temp_c,wz)	
				
			   CALL decomp_2d_fft_3d(uhat,u)	
			   DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)	
				   u(i,j,k)=u(i,j,k)*scalemodes
			   END DO; END DO; END DO
			   CALL decomp_2d_fft_3d(vhat,v)
			   DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   v(i,j,k)=v(i,j,k)*scalemodes
			   END DO; END DO; END DO
			   CALL decomp_2d_fft_3d(what,w)
			   DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   w(i,j,k)=w(i,j,k)*scalemodes
			   END DO; END DO; END DO
			END DO

		END DO
		!calculate vorticity
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			temp(i,j,k)=REAL(wy(i,j,k)-vz(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			omega(i,j,k)=REAL(temp(i,j,k)*temp(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			temp(i,j,k)=REAL(uz(i,j,k)-wx(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			omega(i,j,k)=REAL(temp(i,j,k)*temp(i,j,k)+omega(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			temp(i,j,k)=REAL(vx(i,j,k)-uy(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			omega(i,j,k)=REAL(temp(i,j,k)*temp(i,j,k)+omega(i,j,k),KIND=8)
		END DO ; END DO ; END DO

	   !Calculate Kinetic Energy
	   temp1(1) = sum(u*u)
	   temp1(2) = sum(v*v)
	   temp1(3) = sum(w*w)
	   CALL MPI_ALLREDUCE(temp1(1:3),temp2(1:3),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	   temp3=temp2(1)+temp2(2)+temp2(3)
	   KE(n+1)=dx*dy*dz*temp3/(8.0d0*2.0d0*pi*pi*pi)

	   !Calculate Enstrophy
      	   temp1(1)=sum(omega)
	   CALL MPI_ALLREDUCE(temp1(1),temp3,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      	   Enstrophy(n+1)=dx*dy*dz*temp3/(8.0d0*2.0d0*pi*pi*pi)

	   temp1=0.0;temp2=0.0
	   !Calculate Kinetic Energy Dissipation Rate
	   temp1(1)=sum(ux*ux)
	   temp1(2)=sum(uy*uy)
	   temp1(3)=sum(uz*uz)
	   temp1(4)=sum(vx*vx)
	   temp1(5)=sum(vy*vy)
	   temp1(6)=sum(vz*vz)
	   temp1(7)=sum(wx*wx)
	   temp1(8)=sum(wy*wy)
	   temp1(9)=sum(wz*wz)
	   CALL MPI_ALLREDUCE(temp1,temp2,9,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	   temp3=temp2(1)+temp2(2)+temp2(3)+temp2(4)+temp2(5)+temp2(6)+temp2(7)+temp2(8)+temp2(9)
	   KEdissipationRate(n+1)=dx*dy*dz*temp3/(8.0d0*Re*pi*pi*pi)

     	time(n+1)=time(1)+plotgap*dt*n
     	IF(myid.eq.0) THEN
     		PRINT*,'time=',time(n+1)
		PRINT*,'KE=',KE(n+1)
		PRINT*,'Enstrophy=',Enstrophy(n+1)
		PRINT*,'KEdissipationRate=',KEdissipationRate(n+1)
		PRINT*,'---------------------------------'
      	END IF
      	IF(savedata1==1) THEN
            plotInt=plotInt+1
         	name_config='./data/u'
         	CALL savedata(Nx,Ny,Nz,plotInt,name_config,u,decomp)
         	name_config='./data/v'
         	CALL savedata(Nx,Ny,Nz,plotInt,name_config,v,decomp)
         	name_config='./data/w'
         	CALL savedata(Nx,Ny,Nz,plotInt,name_config,w,decomp)
      	END IF
      END DO

	
    CALL system_clock(finish,count_rate)
    IF (myid.eq.0) then
    	PRINT *, 'Program took', REAL(finish-start)/REAL(count_rate), 'for main timestepping loop'
    END IF	
    
	IF(exactsoln.eqv..true.) THEN
		!Calculate error in final numerical solution
      factor=sqrt(3.0d0)
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			uexact(i,j,k)=u(i,j,k) -&
							(-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
							+sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(numplots+1)/Re))
		END DO; END DO; END DO
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			vexact(i,j,k)=v(i,j,k) -&
						(0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
							-cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(numplots+1)/Re))
		END DO ; END DO ; END DO
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			wexact(i,j,k)=w(i,j,k)-&
						(cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(numplots+1)/Re))
		END DO ; END DO ; END DO
		temp1(1) = maxval(abs(uexact)); temp1(2) = maxval(abs(vexact)); temp1(3) = maxval(abs(wexact))
		CALL MPI_ALLREDUCE(temp1(1:3),temp2(1:3),3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
		temp3=temp2(1)+temp2(2)+temp2(3)
		IF (myid.eq.0) THEN
			PRINT*,'The error at time',time(numplots+1),' is=',temp3
		END IF
	END IF

   !save diagnostics
   IF (savedata1==1) THEN
      IF (myid.eq.0) THEN
         name_config='./data/Kinetic Energy' 
         OPEN(UNIT=13, FILE=name_config,position='append')
         DO n=1,(numplots)
	         write(13,*) KE(n)
         END DO
         CLOSE(13)

         name_config='./data/Enstrophy' 
         OPEN(UNIT=13, FILE=name_config,position='append')
         DO n=1,(numplots)
	         write(13,*)  Enstrophy(n)
         END DO
         CLOSE(13)

         name_config='./data/KEdissipationRate' 
         OPEN(UNIT=13, FILE=name_config,position='append')
         DO n=1,(numplots)
	         write(13,*) KEdissipationRate(n)
         END DO
         CLOSE(13)

         name_config='./data/tdata' 
         OPEN(UNIT=13, FILE=name_config,position='append')
         DO n=1,(numplots)
	         write(13,*) time(n)
         END DO
         CLOSE(13)

			!time and plot number needed for restart
         name_config='./data/restartData'
         OPEN(UNIT=13, FILE=name_config, status='unknown')
         write(13,*) time(numplots+1),plotInt
         CLOSE(13)
      END IF
   END IF

   111 continue
	DEALLOCATE(time,kx,ky,kz,x,y,z,temp1,temp2,&
			  KE,Enstrophy,KEdissipationRate,&
			  u,v,w,&
			  ux,uy,uz,&
			  vx,vy,vz,&
			  wx,wy,wz,&
			  omega,nonlin,temp,&
			  uhat,vhat,what,&
			  omegahat,temp_c,temp_c1,&
			  nonlinhatuh,nonlinhatvh,&
			  nonlinhatwh,phat,&
			  uexact,vexact,wexact,&
			  stat=allocatestatus)
	IF (allocatestatus .ne. 0) THEN
                PRINT *,'STOP'
                STOP
    END IF
    CALL decomp_2d_fft_finalize
    CALL decomp_2d_finalize
    IF (myid.eq.0) THEN
		PRINT *,'Program execution complete'
    END IF
	CALL MPI_FINALIZE(ierr)
END PROGRAM main

   SUBROUTINE getName(name_config,plotInt)
      implicit none
      CHARACTER*100,intent(inout)       :: name_config
      INTEGER(kind=4)               :: plotInt
      CHARACTER*100                 :: number_file
      INTEGER(kind=4)               :: ind

      ind = index(name_config,' ') - 1
      WRITE(number_file,'(i0)') 10000000+plotInt
      number_file = name_config(1:ind)//number_file
      ind = index(number_file,' ') - 1
      name_config = number_file(1:ind)//'.datbin'
   END SUBROUTINE getname
