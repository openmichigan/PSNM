!--------------------------------------------------------------------
!
!
! PURPOSE
!
! This subroutine program calculates the energy for the nonlinear
! sine-Gordon equation in 2 dimensions
! u_{tt}-u_{xx}+u_{yy}=-sin(u)
!
! The energy density is given by
! 0.5u_t^2+0.5u_x^2+0.5u_y^2+(1-cos(u))
!
! AUTHORS
!
! B. Cloutier, B.K. Muite, P. Rigge
! 4 June 2012
!
! INPUT
!
! .. Scalars ..
!  Nx                           = number of modes in x - power of 2 for FFT
!  Ny                           = number of modes in y - power of 2 for FFT
!  planfxy                      = Forward 2d fft plan
!  planbxy                      = Backward 2d fft plan
!  dt                           = timestep
! .. Arrays ..
!  u                            = approximate solution
!  v                            = Fourier transform of approximate solution
!  uold                         = approximate solution
!  vold                         = Fourier transform of approximate solution
!  temp1                        = array to hold temporary values
!  temp2                        = array to hold temporary values
! .. Vectors ..
!  kx                           = fourier frequencies in x direction
!  ky                           = fourier frequencies in y direction
!
! OUTPUT
!
! .. Scalars ..
!  enkin                        = Kinetic energy
!  enstr                        = Strain energy
!  enpot                        = Potential energy
!  en                           = Total energy
!
! LOCAL VARIABLES
!
! .. Scalars ..
!  j                            = loop counter in y direction
!  i                            = loop counter in x direction
!
! FURTHER COMMENTS
! Check that the initial iterate is consistent with the
! boundary conditions for the domain specified
!--------------------------------------------------------------------
! External routines required
!
! External libraries required
!       FFTW3    -- Fast Fourier Transform in the West Library
!                       (http://www.fftw.org/)
!       OpenMP library
SUBROUTINE enercalc(Nx,Ny,planfxy,planbxy,dt,enkin,enstr,&
     enpot,en,kx_p,ky_p,temp1_p,temp2_p,u_p,uold_p) bind(C,name="enercalc")
  USE omp_lib
  USE ISO_C_BINDING
  IMPLICIT NONE
  INTEGER(KIND=C_INT), INTENT(IN)                                        :: Nx,Ny
  REAL(KIND=C_DOUBLE), INTENT(IN)                                        :: dt
  INTEGER(KIND=C_LONG), INTENT(IN)                                       :: planfxy
  INTEGER(KIND=C_LONG), INTENT(IN)                                       :: planbxy
  TYPE(c_ptr), INTENT(IN)                                                :: kx_p,ky_p,u_p,uold_p
  TYPE(c_ptr), INTENT(INOUT)                                             :: temp1_p,temp2_p
  COMPLEX(KIND=C_DOUBLE_COMPLEX), POINTER                               :: kx(:),ky(:)
  REAL   (KIND=C_DOUBLE), POINTER                                        :: u(:,:),uold(:,:)
  COMPLEX(KIND=C_DOUBLE_COMPLEX), POINTER                                :: temp1(:,:),temp2(:,:)
  REAL(KIND=C_DOUBLE), INTENT(OUT)                                       :: enkin,enstr
  REAL(KIND=C_DOUBLE), INTENT(OUT)                                       :: enpot,en
  INTEGER(KIND=4)                                                        :: i,j
  CALL c_f_pointer(kx_p,kx,[Nx])
  CALL c_f_pointer(ky_p,ky,[Ny])
  CALL c_f_pointer(temp1_p,temp1,[Nx,Ny])
  CALL c_f_pointer(temp2_p,temp2,[Nx,Ny])
  CALL c_f_pointer(u_p,u,[Nx,Ny])
  CALL c_f_pointer(uold_p,uold,[Nx,Ny])
  !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
  DO j=1,Ny
  	DO i=1,Nx
	     temp1(i,j)=0.5d0*(uold(i,j)+u(i,j))
	END DO
  END DO
  CALL dfftw_execute_dft_(planfxy,temp1,temp2)
  !.. Strain energy ..
  !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
  DO j=1,Ny
  	DO i=1,Nx
     temp1(i,j)=kx(i)*temp2(i,j)
     END DO
  END DO
  !$OMP END PARALLEL DO
  CALL dfftw_execute_dft_(planbxy,temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny))
  !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
  DO j=1,Ny
  	DO i=1,Nx
     temp2(i,j)=abs(temp2(i,j)/REAL(Nx*Ny,kind(0d0)))**2
     END DO
  END DO
  !$OMP END PARALLEL DO
  CALL dfftw_execute_dft_(planfxy,temp2(1:Nx,1:Ny),temp1(1:Nx,1:Ny))
  enstr=0.5d0*REAL(abs(temp1(1,1)),kind(0d0))/REAL(Nx*Ny,kind(0d0))
  !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
  DO j=1,Ny
  	DO i=1,Nx
     temp1(i,j)=0.5d0*(uold(i,j)+u(i,j))
	END DO
  END DO
  CALL dfftw_execute_dft_(planfxy,temp1,temp2)
  !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
  DO j=1,Ny
  	DO i=1,Nx
     temp1(i,j)=ky(j)*temp2(i,j)
    END DO
  END DO
  !$OMP END PARALLEL DO
  CALL dfftw_execute_dft_(planbxy,temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny))
  !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
  DO j=1,Ny
  	DO i=1,Nx
     temp2(i,j)=abs(temp2(i,j)/REAL(Nx*Ny,kind(0d0)))**2
    END DO
  END DO
  !$OMP END PARALLEL DO
  CALL dfftw_execute_dft_(planfxy,temp2(1:Nx,1:Ny),temp1(1:Nx,1:Ny))
  enstr=enstr+0.5d0*REAL(abs(temp1(1,1)),kind(0d0))/REAL(Nx*Ny,kind(0d0))
  ! .. Kinetic Energy ..
  !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
  DO j=1,Ny
  	DO i=1,Nx
     temp2(i,j)=( abs(u(i,j)-uold(i,j))/dt )**2
	END DO
  END DO
  !$OMP END PARALLEL DO
  CALL dfftw_execute_dft_(planfxy,temp2(1:Nx,1:Ny),temp1(1:Nx,1:Ny))
  enkin=0.5d0*REAL(abs(temp1(1,1)),kind(0d0))/REAL(Nx*Ny,kind(0d0))
  ! .. Potential Energy ..
  !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
  DO j=1,Ny
  	DO i=1,Nx
     temp2(i,j)=1-cos(0.5*u(i,j)+0.5*uold(i,j))
    END DO
  END DO
  !$OMP END PARALLEL DO
  CALL dfftw_execute_dft_(planfxy,temp2(1:Nx,1:Ny),temp1(1:Nx,1:Ny))
  enpot=REAL(abs(temp1(1,1)),kind(0d0))/REAL(Nx*Ny,kind(0d0))
  en=enpot+enkin+enstr
END SUBROUTINE enercalc
