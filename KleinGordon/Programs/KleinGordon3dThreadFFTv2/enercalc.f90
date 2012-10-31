	SUBROUTINE enercalc(Nx,Ny,Nz,planbxy,planbz,planfxy,planfz,dt,Es,&
				modescalereal,enkin,enstr,enpot,en,mythreadid,numthreads,&
				lowxv,lowzv,upxv,upzv,kx,ky,kz,&
				temp1,temp2,vtemp1,vtemp2,v,vold,u,uold)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine program calculates the energy for the nonlinear 
	! Klein-Gordon equation in 3 dimensions
	! u_{tt}-(u_{xx}+u_{yy}+u_{zz})+u=Es*|u|^2u
	!
	! The energy density is given by 
	! 0.5u_t^2+0.5u_x^2+0.5u_y^2+0.5u_z^2+0.5u^2+Es*0.25u^4
	!
	! INPUT 
	!
	! .. Scalars ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  Nz				= number of modes in z - power of 2 for FFT
	!  planfxyz			= Forward 2d fft plan
	!  planbxyz			= Backward 2d fft plan
	!  dt				= timestep
	!  Es				= +1 for focusing, -1 for defocusing
	!  modescalereal	= Number to scale after backward FFT 
	!  numthreads		= total number of OpenMP threads
	!  mythreadid		= thread number
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
	!  kz				= fourier frequencies in z direction
	!  lowxv			= array with entries of array to be used by each thread
	!  upxv				= array with entries of array to be used by each thread
	!  lowzv			= array with entries of array to be used by each thread
	!  upzv				= array with entries of array to be used by each thread
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
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
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
	INTEGER(KIND=4), INTENT(IN)		:: Nx,Ny,Nz,numthreads,mythreadid
	INTEGER(KIND=4), DIMENSION(1:numthreads), INTENT(IN) :: &
												lowxv,upxv,lowzv,upzv											
	REAL(KIND=8), INTENT(IN)						:: dt,Es,modescalereal
	INTEGER(KIND=8), INTENT(IN)						:: planbxy,planbz,planfxy,planfz	
	COMPLEX(KIND=8), DIMENSION(1:Nx),INTENT(IN)		:: kx
	COMPLEX(KIND=8), DIMENSION(1:Ny),INTENT(IN)		:: ky
	COMPLEX(KIND=8), DIMENSION(1:Nz),INTENT(IN)		:: kz
	COMPLEX(KIND=8), DIMENSION(1:Nx,1:Ny,1:Nz),INTENT(IN)	:: u,uold
	COMPLEX(KIND=8), DIMENSION(1:Nz,1:Ny,1:Nx),INTENT(IN)	:: v,vold
	COMPLEX(KIND=8), DIMENSION(1:Nx,1:Ny,1:Nz),INTENT(INOUT)	:: temp1,vtemp1
	COMPLEX(KIND=8), DIMENSION(1:Nz,1:Ny,1:Nx),INTENT(INOUT)	:: temp2,vtemp2
	REAL(KIND=8), INTENT(OUT)							:: enkin,enstr
	REAL(KIND=8), INTENT(OUT)							:: enpot,en	
	INTEGER(KIND=4)										:: i,j,k
	
	!.. Strain energy ..
	DO i=lowxv(mythreadid),upxv(mythreadid)
		DO j=1,Ny
			DO k=1,Nz
				temp2(k,j,i)=0.5d0*kx(i)*(vold(k,j,i)+v(k,j,i))
			END DO
		END DO
	END DO
	!$OMP FLUSH
	CALL backwardfft(Nx,Ny,Nz,planbz,planbxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,temp2,temp1,vtemp1,vtemp2)
	
	DO k=lowzv(mythreadid),upzv(mythreadid)
		DO j=1,Ny
			DO i=1,Nx
				temp1(i,j,k)=abs(temp1(i,j,k)*modescalereal)**2
			END DO
		END DO
	END DO
	!$OMP FLUSH
	CALL forwardfft(Nx,Ny,Nz,planfz,planfxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,temp1,temp2,vtemp1,vtemp2)

	enstr=0.5d0*REAL(abs(temp2(1,1,1)),kind(0d0))
	DO i=lowxv(mythreadid),upxv(mythreadid)
		DO j=1,Ny
			DO k=1,Nz
				temp2(k,j,i)=0.5d0*ky(j)*(vold(k,j,i)+v(k,j,i))
			END DO
		END DO
	END DO
	!$OMP FLUSH
	CALL backwardfft(Nx,Ny,Nz,planbz,planbxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,temp2,temp1,vtemp1,vtemp2)
	
	DO k=lowzv(mythreadid),upzv(mythreadid)
		DO j=1,Ny
			DO i=1,Nx
				temp1(i,j,k)=abs(temp1(i,j,k)*modescalereal)**2
			END DO
		END DO
	END DO
	!$OMP FLUSH
	
	CALL forwardfft(Nx,Ny,Nz,planfz,planfxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,temp1,temp2,vtemp1,vtemp2)
	
	enstr=enstr+0.5d0*REAL(abs(temp2(1,1,1)),kind(0d0))
	DO i=lowxv(mythreadid),upxv(mythreadid)
		DO j=1,Ny
			DO k=1,Nz
				temp2(k,j,i)=0.5d0*kz(k)*(vold(k,j,i)+v(k,j,i))
			END DO
		END DO
	END DO
	!$OMP FLUSH	
	CALL backwardfft(Nx,Ny,Nz,planbz,planbxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,temp2,temp1,vtemp1,vtemp2)
	DO k=lowzv(mythreadid),upzv(mythreadid)
		DO j=1,Ny
			DO i=1,Nx
				temp1(i,j,k)=abs(temp1(i,j,k)*modescalereal)**2
			END DO
		END DO
	END DO
	!$OMP FLUSH
	CALL forwardfft(Nx,Ny,Nz,planfz,planfxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,temp1,temp2,vtemp1,vtemp2)
	
	enstr=enstr+0.5d0*REAL(abs(temp2(1,1,1)),kind(0d0))

	! .. Kinetic Energy ..
	DO k=lowzv(mythreadid),upzv(mythreadid)
		DO j=1,Ny
			DO i=1,Nx
				temp1(i,j,k)=( abs(u(i,j,k)-uold(i,j,k))/dt )**2
			END DO
		END DO
	END DO
	!$OMP FLUSH
	CALL forwardfft(Nx,Ny,Nz,planfz,planfxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,temp1,temp2,vtemp1,vtemp2)
	
	enkin=0.5d0*REAL(abs(temp2(1,1,1)),kind(0d0))
	
	! .. Potential Energy ..
	DO k=lowzv(mythreadid),upzv(mythreadid)
		DO j=1,Ny
			DO i=1,Nx
				temp1(i,j,k)=0.5d0*(abs((u(i,j,k)+uold(i,j,k))*0.50d0))**2&
					-0.125d0*Es*(abs(u(i,j,k))**4+abs(uold(i,j,k))**4)
			END DO
		END DO
	END DO
	!$OMP FLUSH
	CALL forwardfft(Nx,Ny,Nz,planfz,planfxy,mythreadid,numthreads,&
					lowxv,lowzv,upxv,upzv,temp1,temp2,vtemp1,vtemp2)
	
	enpot=REAL(abs(temp2(1,1,1)),kind(0d0))
	
	en=enpot+enkin+enstr
	
	END SUBROUTINE enercalc