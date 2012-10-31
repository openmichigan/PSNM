	SUBROUTINE enercalc(Nx,Ny,planfxy,planbxy,dt,Es,enkin,enstr,&
				enpot,en,kx,ky,temp1,temp2,v,vold,u,uold)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine program calculates the energy for the nonlinear 
	! Klein-Gordon equation in 2 dimensions
	! u_{tt}-u_{xx}+u_{yy}+u=Es*|u|^2u
	!
	! The energy density is given by 
	! 0.5u_t^2+0.5u_x^2+0.5u_y^2+0.5u^2+Es*0.25u^4
	!
	! INPUT 
	!
	! .. Scalars ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  planfxy			= Forward 2d fft plan
	!  planbxy			= Backward 2d fft plan
	!  dt				= timestep
	!  Es				= +1 for focusing, -1 for defocusing
	! .. Arrays ..
	!  u 				= approximate solution
	!  v 				= Fourier transform of approximate solution
	!  uold 			= approximate solution
	!  vold 			= Fourier transform of approximate solution
	!  temp1			= array to hold temporary values
	!  temp2			= array to hold temporary values
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!
	! OUTPUT
	!
	! .. Scalars ..
	!  enkin			= Kinetic energy
	!  enstr			= Strain energy
	!  enpot			= Potential energy
	!  en				= Total energy
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  j				= loop counter in y direction
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
	! 	FFTW3	 -- Fast Fourier Transform in the West Library
	!			(http://www.fftw.org/)
	! 	OpenMP library
	USE omp_lib		 	   
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)							:: Nx,Ny
	REAL(KIND=8), INTENT(IN)							:: dt,Es
	INTEGER(KIND=8), INTENT(IN)							:: planfxy	
	INTEGER(KIND=8), INTENT(IN)							:: planbxy	
	COMPLEX(KIND=8), DIMENSION(1:Nx),INTENT(IN)			:: kx
	COMPLEX(KIND=8), DIMENSION(1:Ny),INTENT(IN)			:: ky
	COMPLEX(KIND=8), DIMENSION(1:Nx,1:Ny),INTENT(IN)	:: u,v,uold,vold
	COMPLEX(KIND=8), DIMENSION(1:Nx,1:Ny),INTENT(INOUT)	:: temp1,temp2
	REAL(KIND=8), INTENT(OUT)							:: enkin,enstr
	REAL(KIND=8), INTENT(OUT)							:: enpot,en	
	INTEGER(KIND=4)										:: j
	
	!.. Strain energy ..
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		temp1(1:Nx,j)=0.5d0*kx(1:Nx)*(vold(1:Nx,j)+v(1:Nx,j))
	END DO
	!$OMP END PARALLEL DO	
	CALL dfftw_execute_dft_(planbxy,temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny))
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		temp1(1:Nx,j)=abs(temp2(1:Nx,j)/REAL(Nx*Ny,kind(0d0)))**2
	END DO
	!$OMP END PARALLEL DO
	CALL dfftw_execute_dft_(planfxy,temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny))
	enstr=0.5d0*REAL(abs(temp2(1,1)),kind(0d0))/REAL(Nx*Ny,kind(0d0))
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		temp1(1:Nx,j)=0.5d0*ky(j)*(vold(1:Nx,j)+v(1:Nx,j))
	END DO
	!$OMP END PARALLEL DO	
	CALL dfftw_execute_dft_(planbxy,temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny))
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		temp1(1:Nx,j)=abs(temp2(1:Nx,j)/REAL(Nx*Ny,kind(0d0)))**2
	END DO
	!$OMP END PARALLEL DO
	CALL dfftw_execute_dft_(planfxy,temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny))
	enstr=enstr+0.5d0*REAL(abs(temp2(1,1)),kind(0d0))/REAL(Nx*Ny,kind(0d0))

	! .. Kinetic Energy ..
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		temp1(1:Nx,j)=( abs(u(1:Nx,j)-uold(1:Nx,j))/dt )**2
	END DO
	!$OMP END PARALLEL DO
	CALL dfftw_execute_dft_(planfxy,temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny))
	enkin=0.5d0*REAL(abs(temp2(1,1)),kind(0d0))/REAL(Nx*Ny,kind(0d0))
	
	! .. Potential Energy ..
	!$OMP PARALLEL DO PRIVATE(j) SCHEDULE(static)
	DO j=1,Ny
		temp1(1:Nx,j)=0.5d0*(abs((u(1:Nx,j)+uold(1:Nx,j))*0.50d0))**2&
					-0.125d0*Es*(abs(u(1:Nx,j))**4+abs(uold(1:Nx,j))**4)
	END DO
	!$OMP END PARALLEL DO
	CALL dfftw_execute_dft_(planfxy,temp1(1:Nx,1:Ny),temp2(1:Nx,1:Ny))
	enpot=REAL(abs(temp2(1,1)),kind(0d0))/REAL(Nx*Ny,kind(0d0))
	
	en=enpot+enkin+enstr
	
	END SUBROUTINE enercalc