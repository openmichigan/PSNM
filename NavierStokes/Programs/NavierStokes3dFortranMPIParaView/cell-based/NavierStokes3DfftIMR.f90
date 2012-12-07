PROGRAM main    
  !-----------------------------------------------------------------------------------
  !
  !
  ! PURPOSE
  !
  ! This program numerically solves the 3D incompressible Navier-Stokes
  ! on a Cubic Domain [0,2pi]x[0,2pi]x[0,2pi] using pseudo-spectral methods and
  ! Implicit Midpoint rule timestepping. The numerical solution is compared to
  ! an exact solution reported by Shapiro 
  !
  ! Analytical Solution:
  !       u(x,y,z,t)=-0.25*(cos(x)sin(y)sin(z)+sin(x)cos(y)cos(z))exp(-t/Re)
  !       v(x,y,z,t)= 0.25*(sin(x)cos(y)sin(z)-cos(x)sin(y)cos(z))exp(-t/Re)
  !       w(x,y,z,t)= 0.5*cos(x)cos(y)sin(z)exp(-t/Re)
  !
  ! .. Parameters ..
  !  Nx                           = number of modes in x - power of 2 for FFT
  !  Ny                           = number of modes in y - power of 2 for FFT
  !  Nz                           = number of modes in z - power of 2 for FFT
  !  Nt                           = number of timesteps to take
  !  Tmax                         = maximum simulation time
  !  FFTW_IN_PLACE        = value for FFTW input 
  !  FFTW_MEASURE         = value for FFTW input
  !  FFTW_EXHAUSTIVE      = value for FFTW input
  !  FFTW_PATIENT         = value for FFTW input    
  !  FFTW_ESTIMATE        = value for FFTW input
  !  FFTW_FORWARD     = value for FFTW input
  !  FFTW_BACKWARD        = value for FFTW input  
  !  pi = 3.14159265358979323846264338327950288419716939937510d0
  !  Re                           = Reynolds number
  ! .. Scalars ..
  !  i                            = loop counter in x direction
  !  j                            = loop counter in y direction
  !  k                            = loop counter in z direction
  !  n                            = loop counter for timesteps direction  
  !  allocatestatus       = error indicator during allocation
  !  count                        = keep track of information written to disk
  !  iol                          = size of array to write to disk
  !  start                        = variable to record start time of program
  !  finish                       = variable to record end time of program
  !  count_rate           = variable for clock count rate
  !  planfxyz                     = Forward 3d fft plan 
  !  planbxyz                     = Backward 3d fft plan
  !  dt                           = timestep
  ! .. Arrays ..
  !  u                            = velocity in x direction
  !  v                            = velocity in y direction
  !  w                            = velocity in z direction
  !  uold                         = velocity in x direction at previous timestep
  !  vold                         = velocity in y direction at previous timestep
  !  wold                         = velocity in z direction at previous timestep
  !  ux                           = x derivative of velocity in x direction
  !  uy                           = y derivative of velocity in x direction
  !  uz                           = z derivative of velocity in x direction
  !  vx                           = x derivative of velocity in y direction
  !  vy                           = y derivative of velocity in y direction
  !  vz                           = z derivative of velocity in y direction
  !  wx                           = x derivative of velocity in z direction
  !  wy                           = y derivative of velocity in z direction
  !  wz                           = z derivative of velocity in z direction
  !  uxold                        = x derivative of velocity in x direction
  !  uyold                        = y derivative of velocity in x direction
  !  uzold                        = z derivative of velocity in x direction
  !  vxold                        = x derivative of velocity in y direction
  !  vyold                        = y derivative of velocity in y direction
  !  vzold                        = z derivative of velocity in y direction
  !  wxold                        = x derivative of velocity in z direction
  !  wyold                        = y derivative of velocity in z direction
  !  wzold                        = z derivative of velocity in z direction
  !  utemp                        = temporary storage of u to check convergence
  !  vtemp                        = temporary storage of u to check convergence
  !  wtemp                        = temporary storage of u to check convergence
  !  temp_r                       = temporary storage for untransformed variables
  !  uhat                         = Fourier transform of u
  !  vhat                         = Fourier transform of v
  !  what                         = Fourier transform of w
  !  rhsuhatfix           = Fourier transform of righthand side for u for timestepping
  !  rhsvhatfix           = Fourier transform of righthand side for v for timestepping
  !  rhswhatfix           = Fourier transform of righthand side for w for timestepping
  !  nonlinuhat           = Fourier transform of nonlinear term for u
  !  nonlinvhat           = Fourier transform of nonlinear term for u
  !  nonlinwhat           = Fourier transform of nonlinear term for u
  !  phat                         = Fourier transform of nonlinear term for pressure, p
  !  temp_c                       = temporary storage for Fourier transforms
  !  realtemp                     = Real storage
  !
  ! .. Vectors ..
  !  kx                           = fourier frequencies in x direction
  !  ky                           = fourier frequencies in y direction
  !  kz                           = fourier frequencies in z direction
  !  x                            = x locations
  !  y                            = y locations
  !  z                            = y locations
  !  time                         = times at which save data
  !  name_config          = array to store filename for data to be saved
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
  ! but is intended as an example only for which all states can be saved
  !
  !--------------------------------------------------------------------------------
  ! External routines required
  ! 
  ! External libraries required
  ! 2DECOMP&FFT -- Fast Fourier Transform in the West Library
  !                       (http://2decomp.org/)


  !---------------------------------------------------------------------------------
  ! mvm comments: All additions marked with !mvm: comment lines.

  USE decomp_2d
  USE decomp_2d_fft
  USE decomp_2d_io
  USE MPI
  !mvm:
  use NSadaptor_module 
  IMPLICIT NONE   
  ! declare variables
  INTEGER(kind=4), PARAMETER              :: Nx=256
  INTEGER(kind=4), PARAMETER              :: Ny=256
  INTEGER(kind=4), PARAMETER              :: Nz=256
  INTEGER(kind=4), PARAMETER              :: Lx=1
  INTEGER(kind=4), PARAMETER              :: Ly=1
  INTEGER(kind=4), PARAMETER              :: Lz=1
  INTEGER(kind=4), PARAMETER              :: Nt=20
  REAL(kind=8), PARAMETER                 :: dt=0.05d0/Nt
  REAL(kind=8), PARAMETER                 :: Re=1.0d0     
  REAL(kind=8), PARAMETER                 :: tol=0.1d0**10
  REAL(kind=8), PARAMETER                 :: theta=0.0d0

  REAL(kind=8), PARAMETER &
       ::  pi=3.14159265358979323846264338327950288419716939937510d0
  REAL(kind=8), PARAMETER         ::      ReInv=1.0d0/REAL(Re,kind(0d0))
  REAL(kind=8), PARAMETER         ::  dtInv=1.0d0/REAL(dt,kind(0d0)) 
  REAL(kind=8)                            :: scalemodes,chg,factor
  REAL(kind=8), DIMENSION(:), ALLOCATABLE                 :: x, y, z, time,mychg,allchg
  COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE  :: u, v, w,&
       ux, uy, uz,&
       vx, vy, vz,&
       wx, wy, wz,&
       uold, uxold, uyold, uzold,&
       vold, vxold, vyold, vzold,&
       wold, wxold, wyold, wzold,&
       utemp, vtemp, wtemp, temp_r

  COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE      ::      kx, ky, kz                                              
  COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE  ::      uhat, vhat, what,&
       rhsuhatfix, rhsvhatfix,&
       rhswhatfix, nonlinuhat,&
       nonlinvhat, nonlinwhat,&
       phat,temp_c
  !mvm: added x,y,z component arrays 
  REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE     ::  realtemp, realtempx, &
       realtempy, realtempz
  ! MPI and 2DECOMP variables
  TYPE(DECOMP_INFO)                               ::  decomp
  INTEGER(kind=MPI_OFFSET_KIND)                   ::  filesize, disp
  INTEGER(kind=4)                                 ::  p_row=0, p_col=0, numprocs, myid, ierr      

  ! variables used for saving data and timing
  INTEGER(kind=4)                                 :: count, iol 
  INTEGER(kind=4)                                 :: i,j,k,n,t,allocatestatus
  INTEGER(kind=4)                                 :: ind, numberfile
  CHARACTER*100                                   :: name_config
  INTEGER(kind=4)                                 :: start, finish, count_rate 

  ! initialisation of 2DECOMP&FFT
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
  ! do automatic domain decomposition
  CALL decomp_2d_init(Nx,Ny,Nz,p_row,p_col)
  ! get information about domain decomposition choosen
  CALL decomp_info_init(Nx,Ny,Nz,decomp)
  ! initialise FFT library
  CALL decomp_2d_fft_init
  IF (myid.eq.0) THEN
     PRINT *,'Grid:',Nx,'X',Ny,'Y',Nz,'Z'
     PRINT *,'dt:',dt
  END IF

  ALLOCATE(x(1:Nx),y(1:Ny),z(1:Nz),time(1:Nt+1),mychg(1:3),allchg(1:3),&
       u(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),& 
       v(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       w(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       ux(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       uy(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       uz(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vx(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vy(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vz(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       wx(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       wy(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       wz(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       uold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       uxold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       uyold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       uzold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vxold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vyold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vzold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       wold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       wxold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       wyold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       wzold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       utemp(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vtemp(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       wtemp(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       temp_r(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       kx(1:Nx),ky(1:Ny),kz(1:Nz),&
       uhat(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       vhat(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       what(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       rhsuhatfix(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       rhsvhatfix(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       rhswhatfix(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       nonlinuhat(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       nonlinvhat(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       nonlinwhat(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       phat(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       temp_c(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       realtemp(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),& 
       realtempx(decomp%xst(1):decomp%xen(1), &
       decomp%xst(2):decomp%xen(2), &
       decomp%xst(3):decomp%xen(3)),&
       realtempy(decomp%xst(1):decomp%xen(1), &
       decomp%xst(2):decomp%xen(2), &
       decomp%xst(3):decomp%xen(3)),&
       realtempz(decomp%xst(1):decomp%xen(1), &
       decomp%xst(2):decomp%xen(2), &
       decomp%xst(3):decomp%xen(3)),&
       stat=AllocateStatus)      
  IF (AllocateStatus .ne. 0) STOP
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
  scalemodes=1.0d0/REAL(Nx*Ny*Nz,kind(0d0))
  IF (myid.eq.0) THEN
     PRINT *,'Setup grid and fourier frequencies'
  END IF

  !initial conditions for Taylor-Green vortex
  !       factor=2.0d0/sqrt(3.0d0)
  !       DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
  !               u(i,j,k)=factor*sin(theta+2.0d0*pi/3.0d0)*sin(x(i))*cos(y(j))*cos(z(k))
  !       END DO; END DO; END DO
  !       DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
  !               v(i,j,k)=factor*sin(theta-2.0d0*pi/3.0d0)*cos(x(i))*sin(y(j))*cos(z(k))
  !       END DO ; END DO ; END DO
  !       DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
  !               w(i,j,k)=factor*sin(theta)*cos(x(i))*cos(y(j))*sin(z(k))
  !       END DO ; END DO ; END DO

  time(1)=0.0d0
  factor=sqrt(3.0d0)
  DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     u(i,j,k)=-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
          +sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
  END DO; END DO; END DO
  DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     v(i,j,k)=0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
          -cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
  END DO; END DO ; END DO
  DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     w(i,j,k)=cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(1)/Re)
  END DO; END DO ; END DO

  CALL decomp_2d_fft_3d(u,uhat,DECOMP_2D_FFT_FORWARD)
  CALL decomp_2d_fft_3d(v,vhat,DECOMP_2D_FFT_FORWARD)
  CALL decomp_2d_fft_3d(w,what,DECOMP_2D_FFT_FORWARD)


  ! derivative of u with respect to x, y, and z 
  DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
     temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
  END DO; END DO ; END DO
  CALL decomp_2d_fft_3d(temp_c,ux,DECOMP_2D_FFT_BACKWARD) 
  DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
     temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
  END DO; END DO ; END DO
  CALL decomp_2d_fft_3d(temp_c,uy,DECOMP_2D_FFT_BACKWARD) 
  DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
     temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
  END DO; END DO ; END DO
  CALL decomp_2d_fft_3d(temp_c,uz,DECOMP_2D_FFT_BACKWARD) 

  ! derivative of v with respect to x, y, and z 
  DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
     temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
  END DO; END DO ; END DO
  CALL decomp_2d_fft_3d(temp_c,vx,DECOMP_2D_FFT_BACKWARD)         
  DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
     temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
  END DO; END DO ; END DO
  CALL decomp_2d_fft_3d(temp_c,vy,DECOMP_2D_FFT_BACKWARD) 
  DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
     temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
  END DO; END DO ; END DO
  CALL decomp_2d_fft_3d(temp_c,vz,DECOMP_2D_FFT_BACKWARD)         

  ! derivative of w with respect to x, y, and z 
  DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
     temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
  END DO; END DO ; END DO
  CALL decomp_2d_fft_3d(temp_c,wx,DECOMP_2D_FFT_BACKWARD)         
  DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
     temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
  END DO; END DO ; END DO
  CALL decomp_2d_fft_3d(temp_c,wy,DECOMP_2D_FFT_BACKWARD)         
  DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
     temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
  END DO; END DO ; END DO
  CALL decomp_2d_fft_3d(temp_c,wz,DECOMP_2D_FFT_BACKWARD)         

  ! save initial data
  n = 0
  !mvm: removed savedata calls from coprocessing version
  !mvm: could have this after the coprocessor initialize and
  !mvm: followed by a call to the coprocessor to write out the initial data
  !mvm: image. Otherwise, these are just wasted calls at this spot.
  call coprocessorinitialize("pipeline.py", 11)

  DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     realtempx(i,j,k)=REAL(wy(i,j,k)-vz(i,j,k),KIND=8)
  END DO; END DO ; END DO

  DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     realtempy(i,j,k)=REAL(uz(i,j,k)-wx(i,j,k),KIND=8)
  END DO; END DO ; END DO

  DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     realtempz(i,j,k)=REAL(vx(i,j,k)-uy(i,j,k),KIND=8)
  END DO; END DO ; END DO

  call NSadaptor(Nx,Ny,Nz,decomp%xst(1),decomp%xen(1), &
          decomp%xst(2),decomp%xen(2),decomp%xst(3),decomp%xen(3), n, n*dt, &
          realtempx, realtempy, realtempz)

  !start timer
  CALL system_clock(start,count_rate)
  ! mvm: Simulation loop starts here
  DO n=1,Nt
     !fixed point
     DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
        uold(i,j,k)=u(i,j,k)
        uxold(i,j,k)=ux(i,j,k)
        uyold(i,j,k)=uy(i,j,k)
        uzold(i,j,k)=uz(i,j,k)
     END DO; END DO ; END DO
     DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
        vold(i,j,k)=v(i,j,k)
        vxold(i,j,k)=vx(i,j,k)
        vyold(i,j,k)=vy(i,j,k)
        vzold(i,j,k)=vz(i,j,k)
     END DO; END DO ; END DO
     DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
        wold(i,j,k)=w(i,j,k)
        wxold(i,j,k)=wx(i,j,k)
        wyold(i,j,k)=wy(i,j,k)
        wzold(i,j,k)=wz(i,j,k)
     END DO; END DO ; END DO
     DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
        rhsuhatfix(i,j,k) = (dtInv+(0.5*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*uhat(i,j,k) 
     END DO; END DO ; END DO
     DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
        rhsvhatfix(i,j,k) = (dtInv+(0.5*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*vhat(i,j,k) 
     END DO; END DO ; END DO
     DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
        rhswhatfix(i,j,k) = (dtInv+(0.5*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*what(i,j,k) 
     END DO; END DO ; END DO

     chg=1
     DO WHILE (chg .gt. tol)
        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
           temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(ux(i,j,k)+uxold(i,j,k))&
                +(v(i,j,k)+vold(i,j,k))*(uy(i,j,k)+uyold(i,j,k))&
                +(w(i,j,k)+wold(i,j,k))*(uz(i,j,k)+uzold(i,j,k)))
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_r,nonlinuhat,DECOMP_2D_FFT_FORWARD)
        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
           temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(vx(i,j,k)+vxold(i,j,k))&
                +(v(i,j,k)+vold(i,j,k))*(vy(i,j,k)+vyold(i,j,k))&
                +(w(i,j,k)+wold(i,j,k))*(vz(i,j,k)+vzold(i,j,k)))
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_r,nonlinvhat,DECOMP_2D_FFT_FORWARD)
        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
           temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(wx(i,j,k)+wxold(i,j,k))&
                +(v(i,j,k)+vold(i,j,k))*(wy(i,j,k)+wyold(i,j,k))&
                +(w(i,j,k)+wold(i,j,k))*(wz(i,j,k)+wzold(i,j,k)))
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_r,nonlinwhat,DECOMP_2D_FFT_FORWARD)
        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           phat(i,j,k)=-1.0d0*( kx(i)*nonlinuhat(i,j,k)&
                +ky(j)*nonlinvhat(i,j,k)&
                +kz(k)*nonlinwhat(i,j,k))&
                /(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)+0.1d0**13)
        END DO; END DO ; END DO

        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           uhat(i,j,k)=(rhsuhatfix(i,j,k)-nonlinuhat(i,j,k)-kx(i)*phat(i,j,k))/&
                (dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
        END DO; END DO ; END DO
        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           vhat(i,j,k)=(rhsvhatfix(i,j,k)-nonlinvhat(i,j,k)-ky(j)*phat(i,j,k))/&
                (dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
        END DO; END DO ; END DO
        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           what(i,j,k)=(rhswhatfix(i,j,k)-nonlinwhat(i,j,k)-kz(k)*phat(i,j,k))/&
                (dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
        END DO; END DO ; END DO

        ! derivative of u with respect to x, y, and z 
        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_c,ux,DECOMP_2D_FFT_BACKWARD) 
        DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_c,uy,DECOMP_2D_FFT_BACKWARD) 
        DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_c,uz,DECOMP_2D_FFT_BACKWARD) 

        ! derivative of v with respect to x, y, and z 
        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_c,vx,DECOMP_2D_FFT_BACKWARD) 
        DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_c,vy,DECOMP_2D_FFT_BACKWARD) 
        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_c,vz,DECOMP_2D_FFT_BACKWARD) 

        ! derivative of w with respect to x, y, and z 
        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_c,wx,DECOMP_2D_FFT_BACKWARD) 
        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_c,wy,DECOMP_2D_FFT_BACKWARD) 
        DO k=decomp%zst(3),decomp%zen(3) ; DO j=decomp%zst(2),decomp%zen(2) ; DO i=decomp%zst(1),decomp%zen(1)
           temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
        END DO; END DO ; END DO
        CALL decomp_2d_fft_3d(temp_c,wz,DECOMP_2D_FFT_BACKWARD) 

        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
           utemp(i,j,k)=u(i,j,k)
        END DO; END DO ; END DO
        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
           vtemp(i,j,k)=v(i,j,k)
        END DO; END DO ; END DO
        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
           wtemp(i,j,k)=w(i,j,k)
        END DO; END DO ; END DO

        CALL decomp_2d_fft_3d(uhat,u,DECOMP_2D_FFT_BACKWARD)    
        CALL decomp_2d_fft_3d(vhat,v,DECOMP_2D_FFT_BACKWARD)    
        CALL decomp_2d_fft_3d(what,w,DECOMP_2D_FFT_BACKWARD)    

        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
           u(i,j,k)=u(i,j,k)*scalemodes
        END DO; END DO ; END DO
        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
           v(i,j,k)=v(i,j,k)*scalemodes
        END DO; END DO ; END DO
        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
           w(i,j,k)=w(i,j,k)*scalemodes
        END DO; END DO ; END DO

        mychg(1) =maxval(abs(utemp-u))
        mychg(2) =maxval(abs(vtemp-v))
        mychg(3) =maxval(abs(wtemp-w))
        CALL MPI_ALLREDUCE(mychg,allchg,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
        chg=allchg(1)+allchg(2)+allchg(3)
        IF (myid.eq.0) THEN
           PRINT *,'chg:',chg
        END IF
     END DO
     time(n+1)=n*dt

     !goto 5100
     IF (myid.eq.0) THEN     
        PRINT *,'time',n*dt
     END IF

     !mvm: populating the arrays sent to the coprocessor
     DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
        realtempx(i,j,k)=REAL(wy(i,j,k)-vz(i,j,k),KIND=8)
     END DO; END DO ; END DO

     DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
        realtempy(i,j,k)=REAL(uz(i,j,k)-wx(i,j,k),KIND=8)
     END DO; END DO ; END DO

     DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
        realtempz(i,j,k)=REAL(vx(i,j,k)-uy(i,j,k),KIND=8)
     END DO; END DO ; END DO

     call NSadaptor(Nx,Ny,Nz,decomp%xst(1),decomp%xen(1), &
             decomp%xst(2),decomp%xen(2),decomp%xst(3),decomp%xen(3), n, n*dt, &
             realtempx, realtempy, realtempz) 
  END DO
  !mvm:
  call coprocessorfinalize()

  CALL system_clock(finish,count_rate)

  IF (myid.eq.0) then
     PRINT *, 'Program took', REAL(finish-start)/REAL(count_rate), 'for main timestepping loop'
  END IF

  IF (myid.eq.0) THEN
     name_config = './data/tdata.dat' 
     OPEN(unit=11,FILE=name_config,status="UNKNOWN")         
     REWIND(11)
     DO n=1,1+Nt
        WRITE(11,*) time(n)
     END DO 
     CLOSE(11)

     name_config = './data/xcoord.dat' 
     OPEN(unit=11,FILE=name_config,status="UNKNOWN")         
     REWIND(11)
     DO i=1,Nx
        WRITE(11,*) x(i)
     END DO
     CLOSE(11)       

     name_config = './data/ycoord.dat' 
     OPEN(unit=11,FILE=name_config,status="UNKNOWN")         
     REWIND(11)
     DO j=1,Ny
        WRITE(11,*) y(j)
     END DO
     CLOSE(11)

     name_config = './data/zcoord.dat' 
     OPEN(unit=11,FILE=name_config,status="UNKNOWN")         
     REWIND(11)
     DO k=1,Nz
        WRITE(11,*) z(k)
     END DO
     CLOSE(11)
     PRINT *,'Saved data'
  END IF

  ! Calculate error in final numerical solution
  DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     utemp(i,j,k)=u(i,j,k) -&
          (-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
          +sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
  END DO; END DO; END DO
  DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     vtemp(i,j,k)=v(i,j,k) -&
          (0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
          -cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
  END DO; END DO ; END DO
  DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
     wtemp(i,j,k)=w(i,j,k)-&
          (cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(Nt+1)/Re))
  END DO; END DO ; END DO
  mychg(1) = maxval(abs(utemp))
  mychg(2) = maxval(abs(vtemp))
  mychg(3) = maxval(abs(wtemp))
  CALL MPI_ALLREDUCE(mychg,allchg,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  chg=allchg(1)+allchg(2)+allchg(3)
  IF (myid.eq.0) THEN
     PRINT*,'The error at the final timestep is',chg
  END IF

  ! clean up 
  CALL decomp_2d_fft_finalize
  CALL decomp_2d_finalize

  DEALLOCATE(x,y,z,time,mychg,allchg,u,v,w,ux,uy,uz,vx,vy,vz,wx,wy,wz,uold,uxold,uyold,uzold,&
       vold,vxold,vyold,vzold,wold,wxold,wyold,wzold,utemp,vtemp,wtemp,&
       temp_r,kx,ky,kz,uhat,vhat,what,rhsuhatfix,rhsvhatfix,&
       rhswhatfix,phat,nonlinuhat,nonlinvhat,nonlinwhat,temp_c,&
       realtemp,stat=AllocateStatus)         
  !mvm: deallocate Coprocessing specific arrays
  deallocate(realtempx, realtempy, realtempz)
  IF (AllocateStatus .ne. 0) STOP
  IF (myid.eq.0) THEN
     PRINT *,'Program execution complete'
  END IF
  CALL MPI_FINALIZE(ierr)         

END PROGRAM main

