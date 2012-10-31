                                                                                                                                                                                                                                                     
 	!--------------------------------------------------------------------
	!
	! PURPOSE
	!
	! This program numerically solves the 2D incompressible Navier-Stokes
	! on a Square Domain [0,1]x[0,1] using pseudo-spectral methods and
	! Crank-Nicolson timestepping. The numerical solution is compared to
	! the exact Taylor-Green Vortex Solution. 
	!
	! AUTHORS
	!
	! B. Cloutier, B.K. Muite, P. Rigge
	! 4 June 2012
	!
	! Periodic free-slip boundary conditions and Initial conditions:
	!	u(x,y,0)=sin(2*pi*x)cos(2*pi*y)
	!	v(x,y,0)=-cos(2*pi*x)sin(2*pi*y)
	! Analytical Solution (subscript denote derivatives):
	!	u(x,y,t)=sin(2*pi*x)cos(2*pi*y)exp(-8*pi^2*t/Re)
	!	v(x,y,t)=-cos(2*pi*x)sin(2*pi*y)exp(-8*pi^2*t/Re)
	!   u_y(x,y,t)=-2*pi*sin(2*pi*x)sin(2*pi*y)exp(-8*pi^2*t/Re)
	!	v_x(x,y,t)=2*pi*sin(2*pi*x)sin(2*pi*y)exp(-8*pi^2*t/Re)
	!	omega=v_x-u_y
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  nplots			= number of plots produced
	!  plotgap			= number of timesteps inbetween plots
	!  Re 				= dimensionless Renold's number
	!  ReInv			= 1/Re for optimization
	!  dt				= timestep size 
	!  dtInv			= 1/dt for optimization
	!  tol				= determines when convergences is reached
	!  scalemodes		= 1/(Nx*Ny), scaling after preforming FFTs
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  n				= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  time				= times at which data is saved
	!  chg				= error at each iteration	
	! .. Arrays ..
	!  omeg				= vorticity	in real space
	!  omeghat			= 2D Fourier transform of vorticity
	!						at next iterate
	!  omegoldhat		= 2D Fourier transform of vorticity at previous
	!						iterate
	!  nl				= nonlinear term 
	!  nlhat			= nonlinear term in Fourier space
	!  nloldhat			= nonlinear term in Fourier space
	!						at previous iterate
	!  omegexact		= taylor-green vorticity at
	!						at final step
	!  psihat			= 2D Fourier transform of streamfunction
	!						at next iteration
	!  temp1/temp2/temp3= reusable complex/real space used for
	!						calculations. This reduces number of
	!						arrays stored.
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  x				= x locations
	!  y				= y locations
	!  name_config		= array to store filename for data to be saved    		
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

	interface cufftExecD2Z
    subroutine cufftExecD2Z(plan, idata, odata) &
          & bind(C,name='cufftExecD2Z')
       use iso_c_binding
       use precision
       integer(c_int),  value  :: plan
       real(fp_kind),   device :: idata(1:nx,1:ny)
       complex(fp_kind),device :: odata(1:nx/2+1,1:ny)
    end subroutine cufftExecD2Z
	end interface cufftExecD2Z

	interface cufftExecZ2D
    subroutine cufftExecZ2D(plan, idata, odata) &
          & bind(C,name='cufftExecZ2D')
       use iso_c_binding
       use precision
       integer(c_int),value:: plan
       complex(fp_kind),device:: idata(1:nx/2+1,1:ny)
       real(fp_kind),device :: odata(1:nx,1:ny)
    end subroutine cufftExecZ2D
	end interface cufftExecZ2D	
	end module cufft
	
	
	PROGRAM main	
	USE precision
	USE cufft	
	USE openacc
	
	IMPLICIT NONE		
    INTEGER(kind=4), PARAMETER 			::  Nx=512
	INTEGER(kind=4), PARAMETER 			::  Ny=512
	REAL(kind=8), PARAMETER			::  dt=0.000125d0
	REAL(kind=8), PARAMETER			::  dtInv=1.0d0/REAL(dt,kind(0d0))   
	REAL(kind=8), PARAMETER	&
		::  pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER			::  Re=1.0d0 		
	REAL(kind=8), PARAMETER			::	ReInv=1.0d0/REAL(Re,kind(0d0))
	REAL(kind=8), PARAMETER			::	tol=0.1d0**10
	REAL(kind=8)					::  scalemodes
	REAL(kind=8)					::  chg
	INTEGER(kind=4), PARAMETER		::	nplots=1, plotgap=20
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE		::  kx,ky						
	REAL(kind=8), DIMENSION(:), ALLOCATABLE			::  x,y,time
	REAL(kind=8), DIMENSION(:,:), ALLOCATABLE		:: omeg,nl, temp2, temp3, omegexact
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE	:: omegoldhat, nloldhat,&
													   omeghat,nlhat, psihat,temp1
	INTEGER(kind=4)								::  i,j,n,t, allocatestatus
	INTEGER(kind=4)								::  pland2z,planz2d
	INTEGER(kind=4)								::  count, iol 
	CHARACTER*100			 					::  name_config
	INTEGER(kind=4)								::  start, finish, count_rate 
	INTEGER(kind=4)								::  ierr,vecsize,gangsize
	INTEGER(kind=8)								::  planfxy,planbxy
	
	vecsize=32
	gangsize=16
  	PRINT *,'Grid:',Nx,'X',Ny
	PRINT *,'dt:',dt	
	ALLOCATE(time(1:nplots+1),kx(1:Nx),ky(1:Ny),x(1:Nx),y(1:Ny),&
			omeg(1:Nx,1:Ny),omegoldhat(1:Nx/2+1,1:Ny),&
			nloldhat(1:Nx/2+1,1:Ny),temp3(1:Nx,1:Ny),omeghat(1:Nx/2+1,1:Ny),&
			nl(1:Nx,1:Ny),nlhat(1:Nx/2+1,1:Ny), psihat(1:Nx/2+1,1:Ny),&
			temp1(1:Nx/2+1,1:Ny),omegexact(1:Nx,1:Ny),temp2(1:Nx,1:Ny),&
			stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 
	PRINT *,'allocated space'
	
	CALL cufftPlan2D(pland2z,nx,ny,CUFFT_D2Z)
	CALL cufftPlan2D(planz2d,nx,ny,CUFFT_Z2D)
 	
	PRINT *,'Setup 2D FFTs'

	! setup fourier frequencies in x-direction
	!$acc data copy(kx,ky,x,y,time,temp3,omeg,nl,temp1,temp2,omegoldhat,nloldhat,omeghat,nlhat,psihat)
	PRINT *, 'Copied arrays over to device'
	!$acc kernels loop
	DO i=1,Nx/2+1
		kx(i)= 2.0d0*pi*cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))  			
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
		x(i)=REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)) 
	END DO
	!$acc end kernels
	! setup fourier frequencies in y-direction
	!$acc kernels loop
	DO j=1,Ny/2+1
		ky(j)= 2.0d0*pi*cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))  			
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
		y(j)=REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)) 
	END DO
	!$acc end kernels
	scalemodes=1.0d0/REAL(Nx*Ny,kind(0d0))
	PRINT *,'Setup grid and fourier frequencies'

	!initial data
	!$acc kernels loop
	DO j=1,Ny
		DO i=1,NX
			omeg(i,j)=4.0d0*pi*sin(2.0d0*pi*x(i))*sin(2.0d0*pi*y(j))!+0.01d0*cos(2.0d0*pi*y(j))
		END DO
	END DO
	!$acc end kernels
	!\hat{\omega^{n,k}}
	CALL cufftExecD2Z(pland2z,omeg,omeghat)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!get initial nonlinear term using omeghat, psihat, u, and v!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!\hat{\psi^{n+1,k+1}}
	!$acc kernels loop gang(gangsize), vector(vecsize)
	DO j=1,Ny
		DO i=1,Nx/2+1
			psihat(i,j)=-omeghat(i,j)/(kx(i)*kx(i)+ky(j)*ky(j) + 0.1d0**14)
		END DO
	END DO
	!$acc end kernels
	!\omega^{n+1,k+1}
	CALL cufftExecZ2D(planz2d,omeghat,omeg)
	
	!get \hat{\psi_y^{n+1,k+1}} used to get u, NOTE: u=\psi_y	
	!$acc kernels loop gang(gangsize), vector(vecsize)
	DO j=1,Ny
		DO i=1,Nx/2+1
			temp1(i,j)=psihat(i,j)*ky(j)*scalemodes
		END DO
	END DO	
	!$acc end kernels
	CALL cufftExecZ2D(planz2d,temp1,temp3) !u
	
	! \hat{\omega_x^{n,k}}
	!$acc kernels loop
	DO j=1,Ny
		DO i=1,Nx/2+1
			temp1(i,j)=omeghat(i,j)*kx(i)
		END DO
	END DO
	!$acc end kernels
	! \omega_x^{n,k}
	CALL cufftExecZ2D(planz2d,temp1,temp2)
	
	! first part nonlinear term in real space
	!$acc kernels loop
	DO j=1,Ny
		DO i=1,Nx
			nl(i,j)=temp3(i,j)*temp2(i,j)
		END DO
	END DO
	!$acc end kernels
	
	!get \hat{\psi_x^{n+1,k+1}} used to get v, NOTE: v=-\psi_x
	!$acc kernels loop gang(gangsize), vector(vecsize)
	DO j=1,Ny
		DO i=1,Nx/2+1
			temp1(i,j)=-psihat(i,j)*kx(i)*scalemodes
		END DO
	END DO
	!$acc end kernels
	CALL cufftExecZ2D(planz2d,temp1,temp3) !v
	
	! \hat{\omega_y^{n,k}}
	!$acc kernels loop
	DO j=1,Ny
		DO i=1,Nx/2+1
			temp1(i,j)=omeghat(i,j)*ky(j)
		END DO
	END DO
	!$acc end kernels
	! \omega_y^{n,k}
	CALL cufftExecZ2D(planz2d,temp1,temp2)

	! get the rest of nonlinear term in real space
	!$acc kernels loop
	DO j=1,Ny
		DO i=1,Nx
			nl(i,j)=(nl(i,j)+temp3(i,j)*temp2(i,j))*scalemodes
		END DO
	END DO
	!$acc end kernels
	! transform nonlinear term into fourier space
	CALL cufftExecD2Z(pland2z,nl,nlhat)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!$acc kernels loop
	DO j=1,Ny
		DO i=1,Nx				
			temp2(i,j)=omeg(i,j)	
		END DO
	END DO
	!$acc end kernels 

	PRINT *,'Got initial data, starting timestepping'	
	time(1)=0.0d0
	CALL system_clock(start,count_rate)
	DO t=1,nplots
		DO n=1,plotgap
			chg=1.0d0
			! save old values(__^{n,k} terms in equation)
			!$acc kernels loop gang(gangsize), vector(vecsize)
			DO j=1,Ny
				DO i=1,Nx/2+1
					nloldhat(i,j)=nlhat(i,j)
				END DO
			END DO
			!$acc end kernels
			!$acc kernels loop gang(gangsize), vector(vecsize)
			DO j=1,Ny
				DO i=1,Nx/2+1
					omegoldhat(i,j)=omeghat(i,j)
				END DO
			END DO
			!$acc end kernels
			DO WHILE (chg>tol)				
				!Crank-Nicolson timestepping to get \hat{\omega^{n+1,k+1}}
				!$acc kernels loop gang(gangsize), vector(vecsize)
				DO j=1,Ny
					DO i=1,Nx/2+1
						omeghat(i,j)=( (dtInv+0.5d0*ReInv*(kx(i)*kx(i)+ky(j)*ky(j)))&
								*omegoldhat(i,j) - 0.5d0*(nloldhat(i,j)+nlhat(i,j)))/&
								(dtInv-0.5d0*ReInv*(kx(i)*kx(i)+ky(j)*ky(j)))  
					END DO
				END DO
				!$acc end kernels
				CALL cufftExecZ2D(planz2d,omeghat,omeg)
				
				! check for convergence
				chg=0.0d0
				!$acc kernels loop gang(gangsize), vector(vecsize)
				DO j=1,Ny
					DO i=1,Nx
						chg=chg+(omeg(i,j)-temp2(i,j))*(omeg(i,j)-temp2(i,j))&
						*scalemodes*scalemodes
					END DO
				END DO
				!$acc end kernels
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!get nonlinear term using omeghat, psihat, u, and v!
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!\hat{\psi^{n+1,k+1}}
				!$acc kernels loop gang(gangsize), vector(vecsize)
				DO j=1,Ny
					DO i=1,Nx/2+1
						psihat(i,j)=-omeghat(i,j)/(kx(i)*kx(i)+ky(j)*ky(j) + 0.1d0**14)
					END DO
				END DO
				!$acc end kernels
				!\omega^{n+1,k+1}
				CALL cufftExecZ2D(planz2d,omeghat,omeg)
				
				!get \hat{\psi_y^{n+1,k+1}} used to get u, NOTE: u=\psi_y	
				!$acc kernels loop gang(gangsize), vector(vecsize)
				DO j=1,Ny
					DO i=1,Nx/2+1
						temp1(i,j)=psihat(i,j)*ky(j)*scalemodes
					END DO
				END DO	
				!$acc end kernels
				CALL cufftExecZ2D(planz2d,temp1,temp3) !u
				
				! \hat{\omega_x^{n,k}}
				!$acc kernels loop
				DO j=1,Ny
					DO i=1,Nx/2+1
						temp1(i,j)=omeghat(i,j)*kx(i)
					END DO
				END DO
				!$acc end kernels
				! \omega_x^{n,k}
				CALL cufftExecZ2D(planz2d,temp1,temp2)
				
				! first part nonlinear term in real space
				!$acc kernels loop
				DO j=1,Ny
					DO i=1,Nx
						nl(i,j)=temp3(i,j)*temp2(i,j)
					END DO
				END DO
				!$acc end kernels
				
				!get \hat{\psi_x^{n+1,k+1}} used to get v, NOTE: v=-\psi_x
				!$acc kernels loop gang(gangsize), vector(vecsize)
				DO j=1,Ny
					DO i=1,Nx/2+1
						temp1(i,j)=-psihat(i,j)*kx(i)*scalemodes
					END DO
				END DO
				!$acc end kernels
				CALL cufftExecZ2D(planz2d,temp1,temp3)
				
				! \hat{\omega_y^{n,k}}
				!$acc kernels loop
				DO j=1,Ny
					DO i=1,Nx/2+1
						temp1(i,j)=omeghat(i,j)*ky(j)
					END DO
				END DO
				!$acc end kernels
				! \omega_y^{n,k}
				CALL cufftExecZ2D(planz2d,temp1,temp2)

				! get the rest of nonlinear term in real space
				!$acc kernels loop
				DO j=1,Ny
					DO i=1,Nx
						nl(i,j)=(nl(i,j)+temp3(i,j)*temp2(i,j))*scalemodes
					END DO
				END DO
				!$acc end kernels
				! transform nonlinear term into fourier space
				CALL cufftExecD2Z(pland2z,nl,nlhat)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				!\omega^{n+1,k+1} is saved for next iteration
				!$acc kernels loop gang(gangsize), vector(vecsize)
				DO j=1,Ny
					DO i=1,Nx
						temp2(i,j)=omeg(i,j)
					END DO
				END DO
				!$acc end kernels
			END DO
		END DO	
		time(t+1)=time(t)+dt*plotgap
		!PRINT *, time(t+1)
	END DO
	CALL system_clock(finish,count_rate)
	PRINT*,'Program took ',REAL(finish-start)/REAL(count_rate),&
			'for Time stepping'
		
	!get exact omega
	!$acc kernels loop gang(gangsize), vector(vecsize)
	DO j=1,Ny
		DO i=1,Nx
			omegexact(i,j)=4.0d0*pi*sin(2.0d0*pi*x(i))*&
				sin(2.0d0*pi*y(j))*exp(-8.0d0*ReInv*pi**2*nplots*plotgap*dt)
		END DO
	END DO
	!$acc end kernels
	!$acc end data
	
	!compute max error
	PRINT *,'Max Error:',maxval(abs(omeg*scalemodes-omegexact))		
	
	!!!!!!!!!!!!!!!!!!!!!!!!
	!copy over data to disk!
	!!!!!!!!!!!!!!!!!!!!!!!!
	write(name_config,'(a,i0,a)') 'omega',1,'.datbin'
	INQUIRE(iolength=iol) omeg(1,1)
	OPEN(unit=11,FILE=name_config,form="unformatted", access="direct",recl=iol) 
	count = 1	
	DO j=1,Ny
		DO i=1,Nx
			WRITE(11,rec=count) omeg(i,j)*scalemodes
			count=count+1
		END DO
	END DO
	CLOSE(11)
	
	name_config = 'time.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO i=1,Nplots+1
		WRITE(11,*) time(i)
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
	!!!!!!!!!!!!!!!!!!!!!!!!

	CALL cufftDestroy(pland2z)
	CALL cufftDestroy(planz2d)	
	
	DEALLOCATE(time,temp1,temp2,temp3,kx,ky,x,y,&
			omeg,omegoldhat,omegexact, nloldhat,&
			omeghat,nl, nlhat, psihat,&
			stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP
	PRINT *,'Program execution complete'	
	
	END PROGRAM main