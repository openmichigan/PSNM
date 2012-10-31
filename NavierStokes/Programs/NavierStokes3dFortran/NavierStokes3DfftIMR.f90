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
	!	u(x,y,z,t)=-0.25*(cos(x)sin(y)sin(z)+sin(x)cos(y)cos(z))exp(-t/Re)
	!	v(x,y,z,t)= 0.25*(sin(x)cos(y)sin(z)-cos(x)sin(y)cos(z))exp(-t/Re)
	!	w(x,y,z,t)= 0.5*cos(x)cos(y)sin(z)exp(-t/Re)
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  Nz				= number of modes in z - power of 2 for FFT
	!  Nt				= number of timesteps to take
	!  Tmax				= maximum simulation time
	!  FFTW_IN_PLACE 	= value for FFTW input 
	!  FFTW_MEASURE 	= value for FFTW input
	!  FFTW_EXHAUSTIVE 	= value for FFTW input
	!  FFTW_PATIENT 	= value for FFTW input    
	!  FFTW_ESTIMATE 	= value for FFTW input
	!  FFTW_FORWARD     = value for FFTW input
	!  FFTW_BACKWARD	= value for FFTW input	
	!  pi = 3.14159265358979323846264338327950288419716939937510d0
	!  Re				= Reynolds number
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
	!  n				= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  count			= keep track of information written to disk
	!  iol				= size of array to write to disk
	!  start			= variable to record start time of program
	!  finish			= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  planfxyz			= Forward 3d fft plan 
	!  planbxyz			= Backward 3d fft plan
	!  dt				= timestep
	! .. Arrays ..
	!  u				= velocity in x direction
	!  v				= velocity in y direction
	!  w				= velocity in z direction
	!  uold				= velocity in x direction at previous timestep
	!  vold				= velocity in y direction at previous timestep
	!  wold				= velocity in z direction at previous timestep
	!  ux				= x derivative of velocity in x direction
	!  uy				= y derivative of velocity in x direction
	!  uz				= z derivative of velocity in x direction
	!  vx				= x derivative of velocity in y direction
	!  vy				= y derivative of velocity in y direction
	!  vz				= z derivative of velocity in y direction
	!  wx				= x derivative of velocity in z direction
	!  wy				= y derivative of velocity in z direction
	!  wz				= z derivative of velocity in z direction
	!  uxold			= x derivative of velocity in x direction
	!  uyold			= y derivative of velocity in x direction
	!  uzold			= z derivative of velocity in x direction
	!  vxold			= x derivative of velocity in y direction
	!  vyold			= y derivative of velocity in y direction
	!  vzold			= z derivative of velocity in y direction
	!  wxold			= x derivative of velocity in z direction
	!  wyold			= y derivative of velocity in z direction
	!  wzold			= z derivative of velocity in z direction
	!  omeg				= vorticity	in real space
	!  omegold			= vorticity in real space at previous
	!						iterate
	!  omegcheck		= store of vorticity at previous iterate
	!  omegoldhat		= 2D Fourier transform of vorticity at previous
	!						iterate
	!  omegoldhat_x		= x-derivative of vorticity in Fourier space
	!						at previous iterate
	!  omegold_x		= x-derivative of vorticity in real space 
	!						at previous iterate
	!  omegoldhat_y 	= y-derivative of vorticity in Fourier space 
	!						at previous iterate
	!  omegold_y		= y-derivative of vorticity in real space
	!						 at previous iterate
	!  nlold			= nonlinear term in real space
	!						at previous iterate
	!  nloldhat			= nonlinear term in Fourier space
	!						at previous iterate
	!  omeghat			= 2D Fourier transform of vorticity
	!						at next iterate
	!  omeghat_x		= x-derivative of vorticity in Fourier space
	!						at next timestep
	!  omeghat_y		= y-derivative of vorticity in Fourier space
	!						at next timestep
	!  omeg_x			= x-derivative of vorticity in real space
	!						at next timestep
	!  omeg_y			= y-derivative of vorticity in real space
	!						at next timestep
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  kz				= fourier frequencies in z direction
	!  x				= x locations
	!  y				= y locations
	!  z				= y locations
	!  time				= times at which save data
	!  name_config		= array to store filename for data to be saved  
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
	! FFTW3	 -- Fast Fourier Transform in the West Library
	!			(http://www.fftw.org/)
	IMPLICIT NONE		
   	!declare variables
   	INTEGER(kind=4), PARAMETER 		:: Nx=64
	INTEGER(kind=4), PARAMETER 		:: Ny=64
	INTEGER(kind=4), PARAMETER 		:: Nz=64
   	INTEGER(kind=4), PARAMETER 		:: Lx=1
	INTEGER(kind=4), PARAMETER 		:: Ly=1
	INTEGER(kind=4), PARAMETER 		:: Lz=1
	INTEGER(kind=4), PARAMETER		:: Nt=20
	REAL(kind=8), PARAMETER			:: dt=0.2d0/Nt
	REAL(kind=8), PARAMETER			:: Re=1.0d0	
	REAL(kind=8), PARAMETER			:: tol=0.1d0**10
	REAL(kind=8), PARAMETER			:: theta=0.0d0

	REAL(kind=8), PARAMETER	&
		::  pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER		::	ReInv=1.0d0/REAL(Re,kind(0d0))
	REAL(kind=8), PARAMETER		:: dtInv=1.0d0/REAL(dt,kind(0d0)) 
	REAL(kind=8)									:: scalemodes,chg,factor
	REAL(kind=8), DIMENSION(:), ALLOCATABLE			:: x, y, z, time
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	:: u, v, w,&
														ux, uy, uz,&
														vx, vy, vz,&
														wx, wy, wz,&
														uold, uxold, uyold, uzold,&
														vold, vxold, vyold, vzold,&
														wold, wxold, wyold, wzold,&
 														utemp, vtemp, wtemp, temp_r
																	
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE		:: kx, ky, kz						
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	:: uhat, vhat, what,&
													rhsuhatfix, rhsvhatfix,&
													rhswhatfix, nonlinuhat,&
													nonlinvhat, nonlinwhat,&
													phat,temp_c
	REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE 	:: realtemp
	!FFTW variables
	INTEGER(kind=4)							:: ierr
	INTEGER(kind=4), PARAMETER      		:: FFTW_IN_PLACE = 8,&
												FFTW_MEASURE = 0,&
										    	FFTW_EXHAUSTIVE = 8,&
												FFTW_PATIENT = 32,&
                    	 				    	FFTW_ESTIMATE = 64
   	INTEGER(kind=4),PARAMETER       		:: FFTW_FORWARD = -1,&
												FFTW_BACKWARD=1	
	INTEGER(kind=8)							:: planfxyz,planbxyz
	
	!variables used for saving data and timing
	INTEGER(kind=4)							:: count, iol 
	INTEGER(kind=4)							:: i,j,k,n,t,allocatestatus
	INTEGER(kind=4)							:: ind, numberfile
	CHARACTER*100			 				:: name_config
	INTEGER(kind=4)							::	start, finish, count_rate 
	
    PRINT *,'Grid:',Nx,'X',Ny,'Y',Nz,'Z'
	PRINT *,'dt:',dt	
	ALLOCATE(x(1:Nx),y(1:Ny),z(1:Nz),time(1:Nt+1),u(1:Nx,1:Ny,1:Nz),& 
 				v(1:Nx,1:Ny,1:Nz), w(1:Nx,1:Ny,1:Nz), ux(1:Nx,1:Ny,1:Nz),&
				uy(1:Nx,1:Ny,1:Nz), uz(1:Nx,1:Ny,1:Nz), vx(1:Nx,1:Ny,1:Nz),&
				vy(1:Nx,1:Ny,1:Nz), vz(1:Nx,1:Ny,1:Nz), wx(1:Nx,1:Ny,1:Nz),&
				wy(1:Nx,1:Ny,1:Nz), wz(1:Nx,1:Ny,1:Nz), uold(1:Nx,1:Ny,1:Nz),&
				uxold(1:Nx,1:Ny,1:Nz), uyold(1:Nx,1:Ny,1:Nz), uzold(1:Nx,1:Ny,1:Nz),&
				vold(1:Nx,1:Ny,1:Nz), vxold(1:Nx,1:Ny,1:Nz), vyold(1:Nx,1:Ny,1:Nz),&
				vzold(1:Nx,1:Ny,1:Nz), wold(1:Nx,1:Ny,1:Nz), wxold(1:Nx,1:Ny,1:Nz),&
				wyold(1:Nx,1:Ny,1:Nz), wzold(1:Nx,1:Ny,1:Nz), utemp(1:Nx,1:Ny,1:Nz),&
 				vtemp(1:Nx,1:Ny,1:Nz), wtemp(1:Nx,1:Ny,1:Nz), temp_r(1:Nx,1:Ny,1:Nz),&
				kx(1:Nx),ky(1:Ny),kz(1:Nz),uhat(1:Nx,1:Ny,1:Nz), vhat(1:Nx,1:Ny,1:Nz),&
	 			what(1:Nx,1:Ny,1:Nz), rhsuhatfix(1:Nx,1:Ny,1:Nz),&
 				rhsvhatfix(1:Nx,1:Ny,1:Nz), rhswhatfix(1:Nx,1:Ny,1:Nz),&
 				nonlinuhat(1:Nx,1:Ny,1:Nz), nonlinvhat(1:Nx,1:Ny,1:Nz),&
 				nonlinwhat(1:Nx,1:Ny,1:Nz), phat(1:Nx,1:Ny,1:Nz),temp_c(1:Nx,1:Ny,1:Nz),&
				realtemp(1:Nx,1:Ny,1:Nz), stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 
	PRINT *,'allocated space'

	CALL dfftw_plan_dft_3d_(planfxyz,Nx,Ny,Nz,temp_r(1:Nx,1:Ny,1:Nz),&
			temp_c(1:Nx,1:Ny,1:Nz),FFTW_FORWARD,FFTW_ESTIMATE)
	CALL dfftw_plan_dft_3d_(planbxyz,Nx,Ny,Nz,temp_c(1:Nx,1:Ny,1:Nz),&
			temp_r(1:Nx,1:Ny,1:Nz),FFTW_BACKWARD,FFTW_ESTIMATE)
	PRINT *,'Setup 3D FFTs'

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
	PRINT *,'Setup grid and fourier frequencies'

	!initial conditions for Taylor-Green vortex
!	factor=2.0d0/sqrt(3.0d0)
!	DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
!		u(i,j,k)=factor*sin(theta+2.0d0*pi/3.0d0)*sin(x(i))*cos(y(j))*cos(z(k))
!	END DO; END DO; END DO
!	DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
!		v(i,j,k)=factor*sin(theta-2.0d0*pi/3.0d0)*cos(x(i))*sin(y(j))*cos(z(k))
!	END DO ; END DO ; END DO
!	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
!		w(i,j,k)=factor*sin(theta)*cos(x(i))*cos(y(j))*sin(z(k))
!	END DO ; END DO ; END DO

	! Initial conditions for exact solution
	time(1)=0.0d0
	factor=sqrt(3.0d0)
	DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
		u(i,j,k)=-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
						+sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
	END DO; END DO; END DO
	DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
		v(i,j,k)=0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
						-cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
	END DO ; END DO ; END DO
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		w(i,j,k)=cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(1)/Re)
	END DO ; END DO ; END DO

	CALL dfftw_execute_dft_(planfxyz,u(1:Nx,1:Ny,1:Nz),uhat(1:Nx,1:Ny,1:Nz))
	CALL dfftw_execute_dft_(planfxyz,v(1:Nx,1:Ny,1:Nz),vhat(1:Nx,1:Ny,1:Nz))
	CALL dfftw_execute_dft_(planfxyz,w(1:Nx,1:Ny,1:Nz),what(1:Nx,1:Ny,1:Nz))
	
	! derivative of u with respect to x, y, and z 
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),ux(1:Nx,1:Ny,1:Nz))
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),uy(1:Nx,1:Ny,1:Nz))
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),uz(1:Nx,1:Ny,1:Nz))

	! derivative of v with respect to x, y, and z 
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),vx(1:Nx,1:Ny,1:Nz))
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),vy(1:Nx,1:Ny,1:Nz))
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),vz(1:Nx,1:Ny,1:Nz))

	! derivative of w with respect to x, y, and z 
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),wx(1:Nx,1:Ny,1:Nz))
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),wy(1:Nx,1:Ny,1:Nz))
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),wz(1:Nx,1:Ny,1:Nz))
	! save initial data
	time(1)=0.0
	n=0
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		realtemp(i,j,k)=REAL(wy(i,j,k)-vz(i,j,k),KIND=8)
	END DO ; END DO ; END DO
	name_config='./data/omegax'
	CALL savedata(Nx,Ny,Nz,n,name_config,realtemp)
	!omegay
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		realtemp(i,j,k)=REAL(uz(i,j,k)-wx(i,j,k),KIND=8)
	END DO ; END DO ; END DO
	name_config='./data/omegay'
	CALL savedata(Nx,Ny,Nz,n,name_config,realtemp)
	!omegaz
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		realtemp(i,j,k)=REAL(vx(i,j,k)-uy(i,j,k),KIND=8)
	END DO ; END DO ; END DO
	name_config='./data/omegaz'
	CALL savedata(Nx,Ny,Nz,n,name_config,realtemp)

	DO n=1,Nt
		!fixed point
		DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
			uold(i,j,k)=u(i,j,k)
			uxold(i,j,k)=ux(i,j,k)
			uyold(i,j,k)=uy(i,j,k)
			uzold(i,j,k)=uz(i,j,k)
		END DO ; END DO ; END DO
		DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
			vold(i,j,k)=v(i,j,k)
			vxold(i,j,k)=vx(i,j,k)
			vyold(i,j,k)=vy(i,j,k)
			vzold(i,j,k)=vz(i,j,k)
		END DO ; END DO ; END DO
		DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
			wold(i,j,k)=w(i,j,k)
			wxold(i,j,k)=wx(i,j,k)
			wyold(i,j,k)=wy(i,j,k)
			wzold(i,j,k)=wz(i,j,k)
		END DO ; END DO ; END DO
		DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
			rhsuhatfix(i,j,k) = (dtInv+(0.5d0*ReInv)*&
			(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*uhat(i,j,k) 
		END DO ; END DO ; END DO
		DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
			rhsvhatfix(i,j,k) = (dtInv+(0.5d0*ReInv)*&
			(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*vhat(i,j,k) 
		END DO ; END DO ; END DO
		DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
			rhswhatfix(i,j,k) = (dtInv+(0.5d0*ReInv)*&
			(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))*what(i,j,k) 
		END DO ; END DO ; END DO
		
		chg=1
		DO WHILE (chg .gt. tol)
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(ux(i,j,k)+uxold(i,j,k))&
										+(v(i,j,k)+vold(i,j,k))*(uy(i,j,k)+uyold(i,j,k))&
										+(w(i,j,k)+wold(i,j,k))*(uz(i,j,k)+uzold(i,j,k)))
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planfxyz,temp_r(1:Nx,1:Ny,1:Nz),nonlinuhat(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(vx(i,j,k)+vxold(i,j,k))&
										+(v(i,j,k)+vold(i,j,k))*(vy(i,j,k)+vyold(i,j,k))&
										+(w(i,j,k)+wold(i,j,k))*(vz(i,j,k)+vzold(i,j,k)))
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planfxyz,temp_r(1:Nx,1:Ny,1:Nz),nonlinvhat(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_r(i,j,k)=0.25d0*((u(i,j,k)+uold(i,j,k))*(wx(i,j,k)+wxold(i,j,k))&
										+(v(i,j,k)+vold(i,j,k))*(wy(i,j,k)+wyold(i,j,k))&
										+(w(i,j,k)+wold(i,j,k))*(wz(i,j,k)+wzold(i,j,k)))
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planfxyz,temp_r(1:Nx,1:Ny,1:Nz),nonlinwhat(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				phat(i,j,k)=-1.0d0*( kx(i)*nonlinuhat(i,j,k)&
									+ky(j)*nonlinvhat(i,j,k)&
									+kz(k)*nonlinwhat(i,j,k))&
									/(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)+0.1d0**13)
			END DO ; END DO ; END DO

			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				uhat(i,j,k)=(rhsuhatfix(i,j,k)-nonlinuhat(i,j,k)-kx(i)*phat(i,j,k))/&
							(dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
			END DO ; END DO ; END DO
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				vhat(i,j,k)=(rhsvhatfix(i,j,k)-nonlinvhat(i,j,k)-ky(j)*phat(i,j,k))/&
							(dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
			END DO ; END DO ; END DO
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				what(i,j,k)=(rhswhatfix(i,j,k)-nonlinwhat(i,j,k)-kz(k)*phat(i,j,k))/&
							(dtInv-(0.5d0*ReInv)*(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))) !*scalemodes
			END DO ; END DO ; END DO

			! derivative of u with respect to x, y, and z 
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),ux(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),uy(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),uz(1:Nx,1:Ny,1:Nz))

			! derivative of v with respect to x, y, and z 
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),vx(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),vy(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),vz(1:Nx,1:Ny,1:Nz))

			! derivative of w with respect to x, y, and z 
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),wx(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),wy(1:Nx,1:Ny,1:Nz))
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
			END DO ; END DO ; END DO
			CALL dfftw_execute_dft_(planbxyz,temp_c(1:Nx,1:Ny,1:Nz),wz(1:Nx,1:Ny,1:Nz))

			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				utemp(i,j,k)=u(i,j,k)
			END DO ; END DO ; END DO
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				vtemp(i,j,k)=v(i,j,k)
			END DO ; END DO ; END DO
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				wtemp(i,j,k)=w(i,j,k)
			END DO ; END DO ; END DO

			CALL dfftw_execute_dft_(planbxyz,uhat(1:Nx,1:Ny,1:Nz),u(1:Nx,1:Ny,1:Nz))
			CALL dfftw_execute_dft_(planbxyz,vhat(1:Nx,1:Ny,1:Nz),v(1:Nx,1:Ny,1:Nz))
			CALL dfftw_execute_dft_(planbxyz,what(1:Nx,1:Ny,1:Nz),w(1:Nx,1:Ny,1:Nz))

			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				u(i,j,k)=u(i,j,k)*scalemodes
			END DO ; END DO ; END DO
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				v(i,j,k)=v(i,j,k)*scalemodes
			END DO ; END DO ; END DO
			DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
				w(i,j,k)=w(i,j,k)*scalemodes
			END DO ; END DO ; END DO
						
			chg =maxval(abs(utemp-u))+maxval(abs(vtemp-v))+maxval(abs(wtemp-w))
			PRINT *,'chg:',chg
		END DO
		time(n+1)=n*dt
		PRINT *,'time',n*dt
		!NOTE: utemp, vtemp, and wtemp are just temporary space that can be used
		!		 instead of creating new arrays.	
		!omegax
		DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
			realtemp(i,j,k)=REAL(wy(i,j,k)-vz(i,j,k),KIND=8)
		END DO ; END DO ; END DO
		name_config='./data/omegax'
		CALL savedata(Nx,Ny,Nz,n,name_config,realtemp)
		!omegay
		DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
			realtemp(i,j,k)=REAL(uz(i,j,k)-wx(i,j,k),KIND=8)
		END DO ; END DO ; END DO
		name_config='./data/omegay'
		CALL savedata(Nx,Ny,Nz,n,name_config,realtemp)
		!omegaz
		DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
			realtemp(i,j,k)=REAL(vx(i,j,k)-uy(i,j,k),KIND=8)
		END DO ; END DO ; END DO
		name_config='./data/omegaz'
		CALL savedata(Nx,Ny,Nz,n,name_config,realtemp)
	END DO

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

	! Calculate error in final numerical solution
	DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
		utemp(i,j,k)=u(i,j,k) -&
						(-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
						+sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
	END DO; END DO; END DO
	DO k=1,Nz; DO j=1,Ny; DO i=1,Nx
		vtemp(i,j,k)=v(i,j,k) -&
					(0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
						-cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
	END DO ; END DO ; END DO
	DO k=1,Nz ; DO j=1,Ny ; DO i=1,Nx
		wtemp(i,j,k)=w(i,j,k)-&
					(cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(Nt+1)/Re))
	END DO ; END DO ; END DO
	chg=maxval(abs(utemp))+maxval(abs(vtemp))+maxval(abs(wtemp))
	PRINT*,'The error at the final timestep is',chg

	CALL dfftw_destroy_plan_(planfxyz)
	CALL dfftw_destroy_plan_(planbxyz)
	DEALLOCATE(x,y,z,time,u,v,w,ux,uy,uz,vx,vy,vz,wx,wy,wz,uold,uxold,uyold,uzold,&
				  vold,vxold,vyold,vzold,wold,wxold,wyold,wzold,utemp,vtemp,wtemp,&
				  temp_r,kx,ky,kz,uhat,vhat,what,rhsuhatfix,rhsvhatfix,&
 				  rhswhatfix,phat,nonlinuhat,nonlinvhat,nonlinwhat,temp_c,&
				  realtemp,stat=AllocateStatus)		
	IF (AllocateStatus .ne. 0) STOP
	PRINT *,'Program execution complete'
	END PROGRAM main

