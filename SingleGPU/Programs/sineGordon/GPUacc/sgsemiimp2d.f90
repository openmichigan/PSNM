!--------------------------------------------------------------------
!
!
! PURPOSE
!
! This program solves nonlinear sine-Gordon equation in 2 dimensions
! u_{tt}-u_{xx}-u_{yy}=-sin(u)
! using a second order implicit-explicit time stepping scheme.
!
! The boundary conditions are u(x=0,y)=u(2*Lx*\pi,y),
!       u(x,y=0)=u(x,y=2*Ly*\pi)
! The initial condition is set in initialdata.f90
!
! AUTHORS
!
! B. Cloutier, B.K. Muite, P. Rigge
! 4 June 2012
!
! .. Parameters ..
!  Nx                           = number of modes in x - power of 2 for FFT
!  Ny                           = number of modes in y - power of 2 for FFT
!  Nt                           = number of timesteps to take
!  plotgap                      = number of timesteps between plots
!  FFTW_IN_PLACE                = value for FFTW input
!  FFTW_MEASURE                 = value for FFTW input
!  FFTW_EXHAUSTIVE              = value for FFTW input
!  FFTW_PATIENT                 = value for FFTW input
!  FFTW_ESTIMATE                = value for FFTW input
!  FFTW_FORWARD                 = value for FFTW input
!  FFTW_BACKWARD                = value for FFTW input
!  pi                           = 3.1415926535...
!  Lx                           = width of box in x direction
!  Ly                           = width of box in y direction
! .. Scalars ..
!  i                            = loop counter in x direction
!  j                            = loop counter in y direction
!  n                            = loop counter for timesteps direction
!  allocatestatus               = error indicator during allocation
!  start                        = variable to record start time of program
!  finish                       = variable to record end time of program
!  count_rate                   = variable for clock count rate
!  planfxy                      = Forward 2d fft plan  (FFTW)
!  planbxy                      = Backward 2d fft plan (FFTW)
!  planf                        = Forward 2d fft plan  (CUFFT)
!  planb                        = Backward 2d fft plan (CUFFT)
!  dt                           = timestep
!  ierr                         = error code
!  plotnum                      = number of plot
! .. Arrays ..
!  u                            = approximate solution
!  uold                         = approximate solution
!  v                            = Fourier transform of approximate solution
!  vold                         = Fourier transform of approximate solution
!  nonlinhat                    = Fourier transform of nonlinear term, sin(u)
!  temp1                        = extra space for energy computation
!  temp2                        = extra space for energy computation
!  savearray                    = temp array to save out to disk
! .. Vectors ..
!  kx                           = fourier frequencies in x direction
!  ky                           = fourier frequencies in y direction
!  x                            = x locations
!  y                            = y locations
!  time                         = times at which save data
!  en                           = total energy
!  enstr                        = strain energy
!  enpot                        = potential energy
!  enkin                        = kinetic energy
!  name_config                  = array to store filename for data to be saved
!
! REFERENCES
!
! ACKNOWLEDGEMENTS
!
! This program is based on example code to demonstrate usage of Fortran and 
! CUDA FFT routines taken from 
! http://cudamusing.blogspot.com/2010/05/CALLing-cufft-from-cuda-fortran.html
! 
! and
!
! http://cudamusing.blogspot.com/search?q=cublas
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
!       getgrid.f90     -- Get initial grid of points
!       initialdata.f90 -- Get initial data
!       enercalc.f90    -- Subroutine to calculate the energy
!       savedata.f90    -- Save initial data
! External libraries required
!       Cuda FFT
!       OpenACC
!       FFTW3           -- Fastest Fourier Transform in the West
!                       (http://www.fftw.org/)
!       OpenMP
module precision
  ! Precision control
  integer, parameter, public :: Single = kind(0.0) ! Single precision
  integer, parameter, public :: Double = kind(0.0d0) ! Double precision
  !
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
  ! cufftPlan2d(cufftHandle *plan, int nx,int ny, cufftType type,int batch)
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
  ! cufftExecD2Z(cufftHandle plan,
  ! cufftDoubleReal    *idata,
  ! cufftDoubleComplex *odata)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface cufftExecD2Z
     subroutine cufftExecD2Z(plan, idata, odata) &
          & bind(C,name='cufftExecD2Z')
       use iso_c_binding
       use precision
       integer(c_int),  value  :: plan
       real(fp_kind),   device :: idata(1:nx,1:ny)
       complex(fp_kind),device :: odata(1:nx,1:ny)
     end subroutine cufftExecD2Z
  end interface cufftExecD2Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! cufftExecD2Z(cufftHandle plan,
  ! cufftDoubleComplex *idata,
  ! cufftDoubleReal    *odata)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  interface cufftExecZ2D
     subroutine cufftExecZ2D(plan, idata, odata) &
          & bind(C,name='cufftExecZ2D')
       use iso_c_binding
       use precision
       integer(c_int),value   :: plan
       complex(fp_kind),device:: idata(1:nx,1:ny)
       real(fp_kind),device   :: odata(1:nx,1:ny)
     end subroutine cufftExecZ2D
  end interface cufftExecZ2D
end module cufft

PROGRAM sg2d
  USE precision
  USE cufft
  USE openacc
  ! Declare variables
  IMPLICIT NONE
  INTEGER(kind=4), PARAMETER                           :: Nx=1024
  INTEGER(kind=4), PARAMETER                           :: Ny=Nx
  INTEGER(kind=4), PARAMETER                           :: Nt=500
  INTEGER(kind=4), PARAMETER                           :: plotgap=Nt+1
  REAL(kind=8), PARAMETER                              :: &
       pi=3.14159265358979323846264338327950288419716939937510d0
  REAL(kind=8), PARAMETER                              :: Lx=5.0d0
  REAL(kind=8), PARAMETER                              :: Ly=5.0d0
  REAL(kind=8)                                         :: dt=0.001d0
  COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE           :: kx,ky
  REAL(kind=8),          DIMENSION(:), ALLOCATABLE     :: x,y
  REAL   (kind=8), DIMENSION(:,:), ALLOCATABLE         :: u,uold
  COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE         :: temp1,temp2,v,vold,nonlinhat
  REAL(kind=8), DIMENSION(:,:), ALLOCATABLE            :: savearray
  REAL(kind=8), DIMENSION(:), ALLOCATABLE              :: time,enkin,enstr,enpot,en
  INTEGER(kind=4)                                      :: ierr,i,j,n,allocatestatus
  INTEGER(kind=4)                                      :: start, finish, count_rate, plotnum
  INTEGER(kind=4), PARAMETER                           :: FFTW_IN_PLACE = 8, FFTW_MEASURE = 0, &
       FFTW_EXHAUSTIVE = 8, FFTW_PATIENT = 32, FFTW_ESTIMATE = 64
  INTEGER(kind=4),PARAMETER                            :: FFTW_FORWARD = -1, FFTW_BACKWARD=1
  INTEGER(kind=8)                                      :: planfxy,planbxy
  CHARACTER*100                                        :: name_config
  INTEGER(kind=4)                                      :: planf,planb
  ! print run information
  PRINT *,"Nx=", Nx
  PRINT *,"Ny=", Ny
  PRINT *,"Nt=", Nt
  PRINT *,"Lx=", Lx
  PRINT *,"Ly=", Ly
  PRINT *,"dt=", dt
  ALLOCATE(kx(1:Nx),ky(1:Ny),x(1:Nx),y(1:Ny),u(1:Nx,1:Ny),uold(1:Nx,1:Ny),&
       v(1:Nx/2+1,1:Ny),vold(1:Nx/2+1,1:Ny),nonlinhat(1:Nx/2+1,1:Ny),&
       savearray(1:Nx,1:Ny),time(1:1+Nt/plotgap),enkin(1:1+Nt/plotgap+1),&
       enstr(1:1+Nt/plotgap+1),enpot(1:1+Nt/plotgap+1),en(1:1+Nt/plotgap),&
       temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny),&
       stat=allocatestatus)
  IF (allocatestatus .ne. 0) stop
  PRINT *,'allocated arrays'
  ! set up cuda ffts
  call cufftPlan2D(planf,nx,ny,CUFFT_D2Z)
  call cufftPlan2D(planb,nx,ny,CUFFT_Z2D)
  ! set up fftw ffts
  CALL dfftw_plan_dft_2d_(planfxy,Nx,Ny,u,temp2,FFTW_FORWARD,FFTW_ESTIMATE)
  CALL dfftw_plan_dft_2d_(planbxy,Nx,Ny,temp2,u,FFTW_BACKWARD,FFTW_ESTIMATE)
  PRINT *,'Setup FFTs'
  ! setup grid, wave numbers
  !$acc data copy(x, y, kx, ky, vold, v, nonlinhat, uold, u)
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
     x(i)=(-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*Lx
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
     y(j)=(-1.0d0 + 2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)))*pi*Ly
  END DO
  !$acc end kernels
  PRINT *,'Got grid and fourier frequencies'
  !$acc kernels loop
  DO j=1,Ny
    DO i=1,Nx
        u(i,j)=0.5d0*exp(-1.0d0*(x(i)**2 +y(j)**2))
    END DO
  END DO
  !$acc end kernels
  !$acc kernels loop
  DO j=1,Ny
    DO i=1,Nx
      uold(i,j)=0.5d0*exp(-1.0d0*(x(i)**2 +y(j)**2))
    END DO
  END DO
  !$acc end kernels
  savearray=REAL(u)
  plotnum=1
  name_config = 'data/u'
  ! CALL savedata(Nx,Ny,plotnum,name_config,savearray) ! disabled for benchmarking
  PRINT *,'data saved'
  !$acc end data
  CALL enercalc(Nx,Ny,planfxy,planbxy,dt,enkin(plotnum),enstr(plotnum),&
       enpot(plotnum),en(plotnum),kx(1:Nx),ky(1:Ny),temp1,temp2,&
       u(1:Nx,1:Ny),uold(1:Nx,1:Ny))
  !$acc data copy(x, y, kx, ky, vold, v, nonlinhat, uold, u)
  call cufftExecD2Z(planf,u,v)
  call cufftExecD2Z(planf,uold,vold)
  PRINT *,'Got initial data, starting timestepping'
  time(plotnum)=0.0d0
  CALL system_clock(start,count_rate)
  DO n=1,Nt
     !$acc kernels loop
     DO j=1,Ny
        DO i=1,Nx
           uold(i,j)=u(i,j)
           u(i,j)=sin(u(i,j))
        END DO
     END DO
     !$acc end kernels
     call cufftExecD2Z(planf,u,nonlinhat)
     !$acc kernels loop
     DO j=1,Ny
        DO i=1,Nx/2+1
           nonlinhat(i,j)=( 0.25*(kx(i)*kx(i) + ky(j)*ky(j))&
                *(2.0d0*v(i,j)+vold(i,j))+(2.0d0*v(i,j)-vold(i,j))/(dt*dt)&
                -nonlinhat(i,j) )/(1/(dt*dt)-0.25*(kx(i)*kx(i) + ky(j)*ky(j)))
           vold(i,j)=v(i,j)
           v(i,j)=nonlinhat(i,j)
           ! prescale nonlinhat
           nonlinhat(i,j)=nonlinhat(i,j)/REAL(Nx*Ny,kind(0d0))
        END DO
     END DO
     !$acc end kernels
     call cufftExecZ2D(planb,nonlinhat,u)
  END DO
  CALL system_clock(finish,count_rate)
  !$acc end data
  PRINT *,'Finished time stepping'
  ! compute energy at the end
  ! savearray=REAL(u(1:Nx,1:Ny),kind(0d0)) ! disabled for benchmarking
  ! CALL savedata(Nx,Ny,plotnum+1,name_config,savearray)
  CALL enercalc(Nx,Ny,planfxy,planbxy,dt,enkin(plotnum+1),enstr(plotnum+1),&
       enpot(plotnum+1),en(plotnum+1),kx,ky,temp1,temp2,u(1:Nx,1:Ny),uold(1:Nx,1:Ny))
  PRINT*,'Program took ',&
       REAL(finish-start,kind(0d0))/REAL(count_rate,kind(0d0)),&
       'for Time stepping'
  CALL saveresults(Nt,plotgap,time(1:1+n/plotgap),en(1:1+n/plotgap+1),&
       enstr(1:1+n/plotgap+1),enkin(1:1+n/plotgap+1),enpot(1:1+n/plotgap+1))

  ! Save times at which output was made in text format
  PRINT *,'Saved data'

  call cufftDestroy(planf)
  call cufftDestroy(planb)
  PRINT *,'Destroy CUFFT Plans'
  call dfftw_destroy_plan_(planbxy)
  call dfftw_destroy_plan_(planfxy)
  PRINT *,'Destroy FFTW Plans'
  DEALLOCATE(kx,ky,x,y,u,uold,time,enkin,enstr,enpot,en,savearray,temp1,temp2,&
       stat=allocatestatus)
  IF (allocatestatus .ne. 0) STOP
  PRINT *,'Deallocated host arrays'
  PRINT *,'Program execution complete'
END PROGRAM sg2d
