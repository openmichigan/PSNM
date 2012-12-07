


	SUBROUTINE enercalc(myid,Nx,Ny,Nz,dt,Es,modescalereal,enkin,enstr,&
				enpot,en,kx,ky,kz,tempu,tempv,v,vold,u,uold,decomp)
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
	!  dt				= timestep
	!  Es				= +1 for focusing, -1 for defocusing
	!  modescalereal	= parameter to scale after doing backward FFT
	!  myid				= Process id
	! .. Arrays ..
	!  u 				= approximate solution
	!  v 				= Fourier transform of approximate solution
	!  uold 			= approximate solution
	!  vold 			= Fourier transform of approximate solution
	!  tempu			= array to hold temporary values - real space
	!  tempv			= array to hold temporary values - fourier space
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction
	!  ky				= fourier frequencies in y direction
	!  kz				= fourier frequencies in z direction
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!					see http://www.2decomp.org/ for more information
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
	! 2DECOMP&FFT	 -- Domain decomposition and Fast Fourier Library
	!			(http://www.2decomp.org/index.html)
	! MPI library

	USE decomp_2d
	USE decomp_2d_fft
	USE decomp_2d_io
	IMPLICIT NONE					 
	INCLUDE 'mpif.h'
	
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)								:: Nx,Ny,Nz,myid
	REAL(KIND=8), INTENT(IN)								:: dt,Es,modescalereal
	TYPE(DECOMP_INFO), INTENT(IN)							:: decomp
	COMPLEX(KIND=8), DIMENSION(decomp%zst(1):decomp%zen(1)),INTENT(IN)	:: kx
	COMPLEX(KIND=8), DIMENSION(decomp%zst(2):decomp%zen(2)),INTENT(IN)	:: ky
	COMPLEX(KIND=8), DIMENSION(decomp%zst(3):decomp%zen(3)),INTENT(IN)	:: kz
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   								decomp%xst(2):decomp%xen(2),&
   								decomp%xst(3):decomp%xen(3)),&
   													INTENT(IN)	:: u,uold
	COMPLEX(KIND=8), DIMENSION(decomp%zst(1):decomp%zen(1),&
            					decomp%zst(2):decomp%zen(2),&
            					decomp%zst(3):decomp%zen(3)),&
            										INTENT(IN)	:: v,vold
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   								decomp%xst(2):decomp%xen(2),&
   								decomp%xst(3):decomp%xen(3)),&
   													INTENT(INOUT) :: tempu
	COMPLEX(KIND=8), DIMENSION(decomp%zst(1):decomp%zen(1),&
            					decomp%zst(2):decomp%zen(2),&
            					decomp%zst(3):decomp%zen(3)),&
            										INTENT(INOUT):: tempv
	REAL(KIND=8), INTENT(OUT)							:: enkin,enstr
	REAL(KIND=8), INTENT(OUT)							:: enpot,en	
	INTEGER(KIND=4)										:: i,j,k
	
	!.. Strain energy ..
	DO k=decomp%zst(3),decomp%zen(3)
		DO j=decomp%zst(2),decomp%zen(2)
			DO i=decomp%zst(1),decomp%zen(1)
				tempv(i,j,k)=0.5d0*kx(i)*(vold(i,j,k)+v(i,j,k))
			END DO
		END DO
	END DO
	CALL decomp_2d_fft_3d(tempv,tempu,DECOMP_2D_FFT_BACKWARD)
								
	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				tempu(i,j,k)=abs(tempu(i,j,k)*modescalereal)**2
			END DO
		END DO
	END DO
	CALL decomp_2d_fft_3d(tempu,tempv,DECOMP_2D_FFT_FORWARD)
	IF(myid.eq.0) THEN
		enstr=0.5d0*REAL(abs(tempv(1,1,1)),kind(0d0))
	END IF
	DO k=decomp%zst(3),decomp%zen(3)
		DO j=decomp%zst(2),decomp%zen(2)
			DO i=decomp%zst(1),decomp%zen(1)
				tempv(i,j,k)=0.5d0*ky(j)*(vold(i,j,k)+v(i,j,k))
			END DO
		END DO
	END DO
	CALL decomp_2d_fft_3d(tempv,tempu,DECOMP_2D_FFT_BACKWARD)
	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				tempu(i,j,k)=abs(tempu(i,j,k)*modescalereal)**2
			END DO
		END DO
	END DO
	CALL decomp_2d_fft_3d(tempu,tempv,DECOMP_2D_FFT_FORWARD)
	IF(myid.eq.0) THEN
		enstr=enstr+0.5d0*REAL(abs(tempv(1,1,1)),kind(0d0))
	END IF
	DO k=decomp%zst(3),decomp%zen(3)
		DO j=decomp%zst(2),decomp%zen(2)
			DO i=decomp%zst(1),decomp%zen(1)
				tempv(i,j,k)=0.5d0*kz(k)*(vold(i,j,k)+v(i,j,k))
			END DO
		END DO
	END DO
	CALL decomp_2d_fft_3d(tempv,tempu,DECOMP_2D_FFT_BACKWARD)
	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				tempu(i,j,k)=abs(tempu(i,j,k)*modescalereal)**2
			END DO
		END DO
	END DO
	CALL decomp_2d_fft_3d(tempu,tempv,DECOMP_2D_FFT_FORWARD)
	IF(myid.eq.0) THEN
		enstr=enstr+0.5d0*REAL(abs(tempv(1,1,1)),kind(0d0))
	END IF
	! .. Kinetic Energy ..
	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				tempu(i,j,k)=( abs(u(i,j,k)-uold(i,j,k))/dt )**2
			END DO
		END DO
	END DO
	CALL decomp_2d_fft_3d(tempu,tempv,DECOMP_2D_FFT_FORWARD)
	IF(myid.eq.0) THEN
		enkin=0.5d0*REAL(abs(tempv(1,1,1)),kind(0d0))
	END IF
	! .. Potential Energy ..
	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				tempu(i,j,k)=0.5d0*(abs((u(i,j,k)+uold(i,j,k))*0.50d0))**2&
					-0.125d0*Es*(abs(u(i,j,k))**4+abs(uold(i,j,k))**4)
			END DO
		END DO
	END DO
	CALL decomp_2d_fft_3d(tempu,tempv,DECOMP_2D_FFT_FORWARD)
	IF(myid.eq.0) THEN
		enpot=REAL(abs(tempv(1,1,1)),kind(0d0))	
		en=enpot+enkin+enstr
	END IF
	END SUBROUTINE enercalc
