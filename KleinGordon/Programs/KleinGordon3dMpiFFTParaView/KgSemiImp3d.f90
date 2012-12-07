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
!       u(x,y=-Ly*pi,z)=u(x,y=Ly*pi,z),u(x,y,z=-Ly*pi)=u(x,y,z=Ly*pi),
! The initial condition is u=0.5*exp(-x^2-y^2-z^2)*sin(10*x+12*y)
!
! .. Parameters ..
!  Nx                           = number of modes in x - power of 2 for FFT
!  Ny                           = number of modes in y - power of 2 for FFT
!  Nz                           = number of modes in z - power of 2 for FFT
!  Nt                           = number of timesteps to take
!  Tmax                         = maximum simulation time
!  plotgap                      = number of timesteps between plots
!  pi = 3.14159265358979323846264338327950288419716939937510d0
!  Lx                           = width of box in x direction
!  Ly                           = width of box in y direction
!  Lz                           = width of box in z direction
!  ES                           = +1 for focusing and -1 for defocusing
! .. Scalars ..
!  i                            = loop counter in x direction
!  j                            = loop counter in y direction
!  k                            = loop counter in z direction
!  n                            = loop counter for timesteps direction  
!  allocatestatus       = error indicator during allocation
!  start                        = variable to record start time of program
!  finish                       = variable to record end time of program
!  count_rate           = variable for clock count rate
!  dt                           = timestep
!  modescalereal        = Number to scale after backward FFT 
!  ierr                         = error code
!  plotnum                      = number of plot
!  myid                         = Process id
!  p_row                        = number of rows for domain decomposition
!  p_col                        = number of columns for domain decomposition
!  filesize                     = total filesize
!  disp                         = displacement to start writing data from
! .. Arrays ..
!  unew                         = approximate solution
!  vnew                         = Fourier transform of approximate solution
!  u                            = approximate solution
!  v                            = Fourier transform of approximate solution
!  uold                         = approximate solution
!  vold                         = Fourier transform of approximate solution
!  nonlin                       = nonlinear term, u^3
!  nonlinhat            = Fourier transform of nonlinear term, u^3
! .. Vectors ..
!  kx                           = fourier frequencies in x direction
!  ky                           = fourier frequencies in y direction
!  kz                           = fourier frequencies in z direction
!  x                            = x locations
!  y                            = y locations
!  z                            = z locations
!  time                         = times at which save data
!  en                           = total energy  
!  enstr                        = strain energy
!  enpot                        = potential energy
!  enkin                        = kinetic energy
!  name_config          = array to store filename for data to be saved                  
!  fftfxyz                      = array to setup 2D Fourier transform
!  fftbxyz                      = array to setup 2D Fourier transform
! .. Special Structures ..
!  decomp                       = contains information on domain decomposition
!                                       see http://www.2decomp.org/ for more information
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
!       getgrid.f90     -- Get initial grid of points
!       initialdata.f90 -- Get initial data
!       enercalc.f90 -- Subroutine to calculate the energy
!       savedata.f90 -- Save initial data
!       storeold.f90 -- Store old data
! External libraries required
!       2DECOMP&FFT      -- Domain decomposition and Fast Fourier Library
!                       (http://www.2decomp.org/index.html)
! MPI library

PROGRAM Kg
  use decomp_2d
  use decomp_2d_fft
  use decomp_2d_io
  !coprocessing:
  use KGadaptor_module

  implicit none
  INCLUDE 'mpif.h'
  ! Declare variables
  INTEGER(kind=4)                         :: Nx, Ny, Nz, Nt, plotgap      
  REAL(kind=8), PARAMETER         :: &
       pi=3.14159265358979323846264338327950288419716939937510d0
  REAL(kind=8)                            ::  Lx,Ly,Lz,Es,dt,starttime,modescalereal      
  COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE      ::  kx,ky,kz 
  REAL(kind=8),    DIMENSION(:), ALLOCATABLE      ::  x,y,z       
  COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  u,nonlin
  COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  v,nonlinhat
  COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  uold
  COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  vold
  COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  unew
  COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE::  vnew
  REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: savearray
  REAL(kind=8), DIMENSION(:), ALLOCATABLE ::  time,enkin,enstr,enpot,en
  INTEGER(kind=4)         ::  ierr,i,j,k,n,allocatestatus,myid,numprocs
  INTEGER(kind=4)         ::  start, finish, count_rate, plotnum
  TYPE(DECOMP_INFO)       ::  decomp
  INTEGER(kind=MPI_OFFSET_KIND) :: filesize, disp
  INTEGER(kind=4) ::  p_row=0, p_col=0    
  CHARACTER*100   ::  name_config
  ! initialisation of 2DECOMP&FFT
  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 

  CALL readinputfile(Nx,Ny,Nz,Nt,plotgap,Lx,Ly,Lz, &
       Es,DT,starttime,myid,ierr)
  ! do automatic domain decomposition
  CALL decomp_2d_init(Nx,Ny,Nz,p_row,p_col)
  ! get information about domain decomposition choosen
  CALL decomp_info_init(Nx,Ny,Nz,decomp)
  ! initialise FFT library
  CALL decomp_2d_fft_init
  ALLOCATE(kx(decomp%zst(1):decomp%zen(1)),&
       ky(decomp%zst(2):decomp%zen(2)),&
       kz(decomp%zst(3):decomp%zen(3)),&
       x(decomp%xst(1):decomp%xen(1)),&
       y(decomp%xst(2):decomp%xen(2)),&
       z(decomp%xst(3):decomp%xen(3)),&
       u(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       v(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       nonlin(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       nonlinhat(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       uold(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vold(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       unew(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       vnew(decomp%zst(1):decomp%zen(1),&
       decomp%zst(2):decomp%zen(2),&
       decomp%zst(3):decomp%zen(3)),&
       savearray(decomp%xst(1):decomp%xen(1),&
       decomp%xst(2):decomp%xen(2),&
       decomp%xst(3):decomp%xen(3)),&
       time(1:1+Nt/plotgap),enkin(1:1+Nt/plotgap),&
       enstr(1:1+Nt/plotgap),enpot(1:1+Nt/plotgap),&
       en(1:1+Nt/plotgap),stat=allocatestatus) 
  IF (allocatestatus .ne. 0) stop 
  IF (myid.eq.0) THEN     
     PRINT *,'allocated arrays'
  END IF
  ! setup fourier frequencies
  CALL getgrid(myid,Nx,Ny,Nz,Lx,Ly,Lz,pi,name_config,x,y,z,kx,ky,kz,decomp)
  IF (myid.eq.0) THEN     
     PRINT *,'Setup grid and fourier frequencies'
  END IF
  CALL initialdata(Nx,Ny,Nz,x,y,z,u,uold,decomp)
  plotnum=1       
  name_config = 'data/u' 
  savearray=REAL(u)
! CALL savedata(Nx,Ny,Nz,plotnum,name_config,savearray,decomp)

  CALL decomp_2d_fft_3d(u,v,DECOMP_2D_FFT_FORWARD)
  CALL decomp_2d_fft_3d(uold,vold,DECOMP_2D_FFT_FORWARD)

  modescalereal=1.0d0/REAL(Nx,KIND(0d0))
  modescalereal=modescalereal/REAL(Ny,KIND(0d0))
  modescalereal=modescalereal/REAL(Nz,KIND(0d0))

  CALL enercalc(myid,Nx,Ny,Nz,dt,Es,modescalereal,&
       enkin(plotnum),enstr(plotnum),&
       enpot(plotnum),en(plotnum),&
       kx,ky,kz,nonlin,nonlinhat,&
       v,vold,u,uold,decomp)

  IF (myid.eq.0) THEN                             
     PRINT *,'Got initial data, starting timestepping'
  END IF
  time(plotnum)=0.0d0+starttime
  CALL system_clock(start,count_rate)     

  !coprocessing:
  call coprocessorinitialize("pipeline.py", 11)

  DO n=1,Nt                                       
     DO k=decomp%xst(3),decomp%xen(3)
        DO j=decomp%xst(2),decomp%xen(2)
           DO i=decomp%xst(1),decomp%xen(1)
              nonlin(i,j,k)=(abs(u(i,j,k))*2)*u(i,j,k)
           END DO
        END DO
     END DO
     CALL decomp_2d_fft_3d(nonlin,nonlinhat,DECOMP_2D_FFT_FORWARD)
     DO k=decomp%zst(3),decomp%zen(3)
        DO j=decomp%zst(2),decomp%zen(2)
           DO i=decomp%zst(1),decomp%zen(1)
              vnew(i,j,k)=&
                   ( 0.25*(kx(i)*kx(i) + ky(j)*ky(j)+ kz(k)*kz(k)-1.0d0)&
                   *(2.0d0*v(i,j,k)+vold(i,j,k))&
                   +(2.0d0*v(i,j,k)-vold(i,j,k))/(dt*dt)&
                   +Es*nonlinhat(i,j,k) )&
                   /(1/(dt*dt)-0.25*(kx(i)*kx(i)+ ky(j)*ky(j)+ kz(k)*kz(k)-1.0d0))
           END DO
        END DO
     END DO
     CALL decomp_2d_fft_3d(vnew,unew,DECOMP_2D_FFT_BACKWARD)
     ! normalize result
     DO k=decomp%xst(3),decomp%xen(3)
        DO j=decomp%xst(2),decomp%xen(2)
           DO i=decomp%xst(1),decomp%xen(1)
              unew(i,j,k)=unew(i,j,k)*modescalereal
           END DO
        END DO
     END DO
     IF (mod(n,plotgap)==0) THEN
        plotnum=plotnum+1
        time(plotnum)=n*dt+starttime
        IF (myid.eq.0) THEN
           PRINT *,'time',n*dt+starttime
        END IF
        CALL enercalc(myid,Nx,Ny,Nz,dt,Es,modescalereal,&
             enkin(plotnum),enstr(plotnum),&
             enpot(plotnum),en(plotnum),&
             kx,ky,kz,nonlin,nonlinhat,&
             vnew,v,unew,u,decomp)
        savearray=REAL(unew,kind(0d0))
!       CALL savedata(Nx,Ny,Nz,plotnum,name_config,savearray,decomp)
     END IF
     !coprocessing:
!    print *, "n*dt+starttime: ", n*dt+starttime
     call KGadaptor(Nx, Ny, Nz, decomp%xst(1), decomp%xen(1), &
                  decomp%xst(2), decomp%xen(2), decomp%xst(3), decomp%xen(3), &
                  n, n*dt+starttime, savearray)      

     ! .. Update old values ..
     CALL storeold(Nx,Ny,Nz,unew,u,uold,vnew,v,vold,decomp)
  END DO
  !coprocessing:
  call coprocessorfinalize()

  CALL system_clock(finish,count_rate)
  IF (myid.eq.0) THEN
     PRINT *,'Finished time stepping'
     PRINT*,'Program took ',&
          REAL(finish-start,kind(0d0))/REAL(count_rate,kind(0d0)),&
          'for Time stepping'
     CALL saveresults(Nt,plotgap,time,en,enstr,enkin,enpot)          
     ! Save times at which output was made in text format
     PRINT *,'Saved data'
  END IF
  CALL decomp_2d_fft_finalize
  CALL decomp_2d_finalize

  DEALLOCATE(kx,ky,kz,x,y,z,u,v,nonlin,nonlinhat,savearray,&
       uold,vold,unew,vnew,time,enkin,enstr,enpot,en,&
       stat=allocatestatus)    
  IF (allocatestatus .ne. 0) STOP 
  IF (myid.eq.0) THEN
     PRINT *,'Deallocated arrays'
     PRINT *,'Program execution complete'
  END IF
  CALL MPI_FINALIZE(ierr)

END PROGRAM Kg
