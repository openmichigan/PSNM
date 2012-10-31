	PROGRAM main
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program numerically solves the 2D incompressible Navier-Stokes
	! on a Square Domain [0,1]x[0,1] using pseudo-spectral methods and
	! Crank-Nicolson timestepping. The numerical solution is compared to
	! the exact Taylor-Green Vortex Solution. 
	!
	! Periodic free-slip boundary conditions and Initial conditions:
	!	u(x,y,0)=sin(2*pi*x)cos(2*pi*y)
	!	v(x,y,0)=-cos(2*pi*x)sin(2*pi*y)
	! Analytical Solution:
	!	u(x,y,t)=sin(2*pi*x)cos(2*pi*y)exp(-8*pi^2*nu*t)
	!	v(x,y,t)=-cos(2*pi*x)sin(2*pi*y)exp(-8*pi^2*nu*t)
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
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
	!  mu				= viscosity
	!  rho				= density
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  n				= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  count			= keep track of information written to disk
	!  iol				= size of array to write to disk
	!  start			= variable to record start time of program
	!  finish			= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  planfx			= Forward 1d fft plan in x
	!  planbx			= Backward 1d fft plan in x
	!  planfy			= Forward 1d fft plan in y
	!  planby			= Backward 1d fft plan in y
	!  dt				= timestep
	! .. Arrays ..
	!  u				= velocity in x direction
	!  uold				= velocity in x direction at previous timestep
	!  v				= velocity in y direction
	!  vold				= velocity in y direction at previous timestep
	!  u_y				= y derivative of velocity in x direction
	!  v_x				= x derivative of velocity in y direction
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
	!  kxx				= square of fourier frequencies in x direction
	!  kyy				= square of fourier frequencies in y direction
	!  x				= x locations
	!  y				= y locations
	!  time				= times at which save data
	!  name_config		= array to store filename for data to be saved    		
	!  fftfx			= array to setup x Fourier transform
	!  fftbx			= array to setup y Fourier transform
	! REFERENCES
	!
	! ACKNOWLEDGEMENTS
	!
	! ACCURACY
	!		
	! ERROR INDICATORS AND WARNINGS
	!
	! FURTHER COMMENTS
	! This program has not been optimized to use the least amount of memory
	! but is intended as an example only for which all states can be saved
	!--------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	! FFTW3	 -- Fast Fourier Transform in the West Library
	!			(http://www.fftw.org/)				 	   
	! declare variables
	
	IMPLICIT NONE		
	INTEGER(kind=4), PARAMETER 	::  Nx=256			
	INTEGER(kind=4), PARAMETER 	::  Ny=256			
	REAL(kind=8), PARAMETER		::  dt=0.00125	    
	REAL(kind=8), PARAMETER	&
		::  pi=3.14159265358979323846264338327950288419716939937510
	REAL(kind=8), PARAMETER		::  rho=1.0d0 		
	REAL(kind=8), PARAMETER		::  mu=1.0d0		
	REAL(kind=8), PARAMETER		::	tol=0.1d0**10		
	REAL(kind=8)				::	chg		
	INTEGER(kind=4), PARAMETER	::	nplots=50
	REAL(kind=8), DIMENSION(:), ALLOCATABLE			::	time
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE		::  kx,kxx				
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE		::  ky,kyy				
	REAL(kind=8), DIMENSION(:), ALLOCATABLE			::  x					
	REAL(kind=8), DIMENSION(:), ALLOCATABLE			::  y						
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE	:: & 
		u,uold,v,vold,u_y,v_x,omegold, omegcheck, omeg,&
		omegoldhat, omegoldhat_x, omegold_x,&	
	  	omegoldhat_y, omegold_y, nlold, nloldhat,&
		omeghat, omeghat_x, omeghat_y, omeg_x, omeg_y,&
		nl, nlhat, psihat, psihat_x, psi_x, psihat_y, psi_y
	REAL(kind=8),DIMENSION(:,:), ALLOCATABLE	::	uexact_y,vexact_x,omegexact
	INTEGER(kind=4)								::  i,j,k,n, allocatestatus, count, iol
	INTEGER(kind=4)			 					::	start, finish, count_rate
	INTEGER(kind=4), PARAMETER      			::  FFTW_IN_PLACE = 8, FFTW_MEASURE = 0, &
										    		FFTW_EXHAUSTIVE = 8, FFTW_PATIENT = 32,    &
                    	 				    		FFTW_ESTIMATE = 64
    INTEGER(kind=4),PARAMETER       			::  FFTW_FORWARD = -1, FFTW_BACKWARD=1	
	COMPLEX(kind=8), DIMENSION(:,:), ALLOCATABLE	::  fftfx,fftbx
	INTEGER(kind=8)								::  planfxy,planbxy
	CHARACTER*100			 					::  name_config
	    		
	CALL system_clock(start,count_rate)		
	ALLOCATE(time(1:nplots),kx(1:Nx),kxx(1:Nx),ky(1:Ny),kyy(1:Ny),x(1:Nx),y(1:Ny),&
			u(1:Nx,1:Ny),uold(1:Nx,1:Ny),v(1:Nx,1:Ny),vold(1:Nx,1:Ny),u_y(1:Nx,1:Ny),&
			v_x(1:Nx,1:Ny),omegold(1:Nx,1:Ny),omegcheck(1:Nx,1:Ny), omeg(1:Nx,1:Ny),&
			omegoldhat(1:Nx,1:Ny),omegoldhat_x(1:Nx,1:Ny), omegold_x(1:Nx,1:Ny), &
			omegoldhat_y(1:Nx,1:Ny),omegold_y(1:Nx,1:Ny), nlold(1:Nx,1:Ny), nloldhat(1:Nx,1:Ny),&
			omeghat(1:Nx,1:Ny), omeghat_x(1:Nx,1:Ny), omeghat_y(1:Nx,1:Ny), omeg_x(1:Nx,1:Ny),&
			omeg_y(1:Nx,1:Ny), nl(1:Nx,1:Ny), nlhat(1:Nx,1:Ny), psihat(1:Nx,1:Ny), &
			psihat_x(1:Nx,1:Ny), psi_x(1:Nx,1:Ny), psihat_y(1:Nx,1:Ny), psi_y(1:Nx,1:Ny),&
			uexact_y(1:Nx,1:Ny), vexact_x(1:Nx,1:Ny), omegexact(1:Nx,1:Ny),fftfx(1:Nx,1:Ny),&
			fftbx(1:Nx,1:Ny),stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 
	PRINT *,'allocated space'
		
	! set up ffts
	CALL dfftw_plan_dft_2d_(planfxy,Nx,Ny,fftfx(1:Nx,1:Ny),fftbx(1:Nx,1:Ny),&
    		FFTW_FORWARD,FFTW_EXHAUSTIVE)
	CALL dfftw_plan_dft_2d_(planbxy,Nx,Ny,fftbx(1:Nx,1:Ny),fftfx(1:Nx,1:Ny),&
			FFTW_BACKWARD,FFTW_EXHAUSTIVE)
		
	! setup fourier frequencies in x-direction
	DO i=1,1+Nx/2
		kx(i)= 2.0d0*pi*cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))  			
	END DO
	kx(1+Nx/2)=0.0d0
	DO i = 1,Nx/2 -1
		kx(i+1+Nx/2)=-kx(1-i+Nx/2)
	END DO
	DO i=1,Nx
		kxx(i)=kx(i)*kx(i)
	END DO
	DO i=1,Nx
		x(i)=REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)) 
	END DO
		
	! setup fourier frequencies in y-direction
	DO j=1,1+Ny/2
		ky(j)= 2.0d0*pi*cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))  			
	END DO
	ky(1+Ny/2)=0.0d0
	DO j = 1,Ny/2 -1
		ky(j+1+Ny/2)=-ky(1-j+Ny/2)
	END DO
	DO j=1,Ny
		kyy(j)=ky(j)*ky(j)
	END DO
	DO j=1,Ny
		y(j)=REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)) 
	END DO
	PRINT *,'Setup grid and fourier frequencies'
	
		
	DO j=1,Ny
		DO i=1,Nx
			u(i,j)=sin(2.0d0*pi*x(i))*cos(2.0d0*pi*y(j))
			v(i,j)=-cos(2.0d0*pi*x(i))*sin(2.0d0*pi*y(j))
			u_y(i,j)=-2.0d0*pi*sin(2.0d0*pi*x(i))*sin(2.0d0*pi*y(j))
			v_x(i,j)=2.0d0*pi*sin(2.0d0*pi*x(i))*sin(2.0d0*pi*y(j))
			omeg(i,j)=v_x(i,j)-u_y(i,j)
		END DO
	END DO
	
	! Vorticity to Fourier Space
	CALL dfftw_execute_dft_(planfxy,omeg(1:Nx,1:Ny),omeghat(1:Nx,1:Ny))	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!Initial nonlinear term !!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	! obtain \hat{\omega}_x^{n,k}
	DO j=1,Ny
		omeghat_x(1:Nx,j)=omeghat(1:Nx,j)*kx(1:Nx)
	END DO
	! obtain \hat{\omega}_y^{n,k}
	DO i=1,Nx
		omeghat_y(i,1:Ny)=omeghat(i,1:Ny)*ky(1:Ny)
	END DO
	! convert to real space 
	CALL dfftw_execute_dft_(planbxy,omeghat_x(1:Nx,1:Ny),omeg_x(1:Nx,1:Ny))
	CALL dfftw_execute_dft_(planbxy,omeghat_y(1:Nx,1:Ny),omeg_y(1:Nx,1:Ny))
	! compute nonlinear term in real space
	DO j=1,Ny
		nl(1:Nx,j)=u(1:Nx,j)*omeg_x(1:Nx,j)/REAL(Nx*Ny,kind(0d0))+&
								v(1:Nx,j)*omeg_y(1:Nx,j)/REAL(Nx*Ny,kind(0d0))	
	END DO
	CALL dfftw_execute_dft_(planfxy,nl(1:Nx,1:Ny),nlhat(1:Nx,1:Ny)) 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	time(1)=0.0d0
	PRINT *,'Got initial data, starting timestepping'
	DO n=1,nplots
		chg=1
		! save old values
		uold(1:Nx,1:Ny)=u(1:Nx,1:Ny)
		vold(1:Nx,1:Ny)=v(1:Nx,1:Ny)
		omegold(1:Nx,1:Ny)=omeg(1:Nx,1:Ny)
		omegcheck(1:Nx,1:Ny)=omeg(1:Nx,1:Ny)		
		omegoldhat(1:Nx,1:Ny)=omeghat(1:Nx,1:Ny)
		nloldhat(1:Nx,1:Ny)=nlhat(1:Nx,1:Ny)
		DO WHILE (chg>tol) 
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!nonlinear fixed (n,k+1)!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			! obtain \hat{\omega}_x^{n+1,k}
			DO j=1,Ny
				omeghat_x(1:Nx,j)=omeghat(1:Nx,j)*kx(1:Nx) 
			END DO
			! obtain \hat{\omega}_y^{n+1,k}
			DO i=1,Nx
				omeghat_y(i,1:Ny)=omeghat(i,1:Ny)*ky(1:Ny)
			END DO
			! convert back to real space 
			CALL dfftw_execute_dft_(planbxy,omeghat_x(1:Nx,1:Ny),omeg_x(1:Nx,1:Ny)) 
			CALL dfftw_execute_dft_(planbxy,omeghat_y(1:Nx,1:Ny),omeg_y(1:Nx,1:Ny))
			! calculate nonlinear term in real space
			DO j=1,Ny
				nl(1:Nx,j)=u(1:Nx,j)*omeg_x(1:Nx,j)/REAL(Nx*Ny,kind(0d0))+&
					v(1:Nx,j)*omeg_y(1:Nx,j)/REAL(Nx*Ny,kind(0d0))
			END DO
			! convert back to fourier
			CALL dfftw_execute_dft_(planfxy,nl(1:Nx,1:Ny),nlhat(1:Nx,1:Ny)) 
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			! obtain \hat{\omega}^{n+1,k+1} with Crank Nicolson timestepping
			DO j=1,Ny
				omeghat(1:Nx,j)=( (1.0d0/dt+0.5d0*(mu/rho)*(kxx(1:Nx)+kyy(j)))&
						*omegoldhat(1:Nx,j) - 0.5d0*(nloldhat(1:Nx,j)+nlhat(1:Nx,j)))/&
				     	(1.0d0/dt-0.5d0*(mu/rho)*(kxx(1:Nx)+kyy(j)))   
			END DO
			
			! calculate \hat{\psi}^{n+1,k+1}
			DO j=1,Ny
				psihat(1:Nx,j)=-omeghat(1:Nx,j)/(kxx(1:Nx)+kyy(j))
			END DO
			psihat(1,1)=0.0d0
		    psihat(Nx/2+1,Ny/2+1)=0.0d0
			psihat(Nx/2+1,1)=0.0d0
			psihat(1,Ny/2+1)=0.0d0
			
			! obtain \psi_x^{n+1,k+1} and \psi_y^{n+1,k+1}
			DO j=1,Ny
				psihat_x(1:Nx,j)=psihat(1:Nx,j)*kx(1:Nx)
			END DO
			CALL dfftw_execute_dft_(planbxy,psihat_x(1:Nx,1:Ny),psi_x(1:Nx,1:Ny)) 
			DO i=1,Nx
				psihat_y(i,1:Ny)=psihat(i,1:Ny)*ky(1:Ny)
			END DO
			CALL dfftw_execute_dft_(planbxy,psihat_y(1:Ny,1:Ny),psi_y(1:Ny,1:Ny)) 
			DO j=1,Ny
				psi_x(1:Nx,j)=psi_x(1:Nx,j)/REAL(Nx*Ny,kind(0d0))
				psi_y(1:Nx,j)=psi_y(1:Nx,j)/REAL(Nx*Ny,kind(0d0))
			END DO	
			
			! obtain \omega^{n+1,k+1}
			CALL dfftw_execute_dft_(planbxy,omeghat(1:Nx,1:Ny),omeg(1:Nx,1:Ny)) 
			DO j=1,Ny
				omeg(1:Nx,j)=omeg(1:Nx,j)/REAL(Nx*Ny,kind(0d0))
			END DO
			
			! obtain u^{n+1,k+1} and v^{n+1,k+1} using stream function (\psi) in real space
			DO j=1,Ny
				u(1:Nx,j)=psi_y(1:Nx,j)
				v(1:Nx,j)=-psi_x(1:Nx,j)
			END DO	
			
			! check for convergence	
			chg=maxval(abs(omeg-omegcheck)) 
			! saves {n+1,k+1} to {n,k} for next iteration
			omegcheck=omeg 	
		END DO
		time(n+1)=time(n)+dt
		PRINT *,'TIME ',time(n+1)
	END DO
		
	DO j=1,Ny
		DO i=1,Nx
			uexact_y(i,j)=-2.0d0*pi*sin(2.0d0*pi*x(i))*sin(2.0d0*pi*y(j))*&
								exp(-8.0d0*mu*(pi**2)*nplots*dt)
			vexact_x(i,j)=2.0d0*pi*sin(2.0d0*pi*x(i))*sin(2.0d0*pi*y(j))*&
								exp(-8.0d0*mu*(pi**2)*nplots*dt)
			omegexact(i,j)=vexact_x(i,j)-uexact_y(i,j)
		END DO
	END DO
		
	name_config = 'omegafinal.datbin' 
	INQUIRE(iolength=iol) omegexact(1,1)
	OPEN(unit=11,FILE=name_config,form="unformatted", access="direct",recl=iol) 
	count = 1	
	DO j=1,Ny
		DO i=1,Nx
			WRITE(11,rec=count) REAL(omeg(i,j),KIND(0d0))
			count=count+1
		END DO
	END DO
	CLOSE(11)
		
	name_config = 'omegaexactfinal.datbin' 
	OPEN(unit=11,FILE=name_config,form="unformatted", access="direct",recl=iol) 
	count = 1	
	DO j=1,Ny
		DO i=1,Nx
			WRITE(11,rec=count) omegexact(i,j)
			count=count+1
		END DO
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
	
	CALL dfftw_destroy_plan_(planfxy)
	CALL dfftw_destroy_plan_(planbxy)
	CALL dfftw_cleanup_()
		
	DEALLOCATE(time,kx,kxx,ky,kyy,x,y,&
			u,uold,v,vold,u_y,v_x,omegold, omegcheck, omeg, &
			omegoldhat, omegoldhat_x, omegold_x,&	
			omegoldhat_y, omegold_y, nlold, nloldhat,&
			omeghat, omeghat_x, omeghat_y, omeg_x, omeg_y,&
			nl, nlhat, psihat, psihat_x, psi_x, psihat_y, psi_y,&
			uexact_y,vexact_x,omegexact, &
			fftfx,fftbx,stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 
	PRINT *,'Program execution complete'
	END PROGRAM main

