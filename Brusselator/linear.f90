SUBROUTINE linear(a,b,Du,Dv,dt,modescalereal,coef,kx,ky,kz,u,v,uhat,vhat,decomp) 
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine solves the linear equation in a splitting 
	! method for the Brusselator equation in 3 dimensions
	! u_t = a + u^2v - (b+1)u + Du (u_xx + u_yy + u_zz)
    ! v_t = bu - u^2v + Dv (v_xx + v_yy + v_zz)
	!
	!
	! INPUT
	!
	! .. Parameters ..
	!  a                = equation parameter
	!  b                = equation parameter
	!  Du               = equation diffusion coefficient
	!  Dv               = equation diffusion coefficient
    !  modescalereal    = scaling of terms after forward and backward FFT
	!  dt				= timestep
	!  coef				= (complex) timestep coefficient
	! .. Arrays ..
	!  kx               = Fourier wave numbers in x direction
	!  ky               = Fourier wave numbers in y direction
	!  kz               = Fourier wave numbers in z direction	
	!  u				= u field
	!  v				= v field
	!  uhat				= storage for u field in Fourier space 
	!  vhat				= storage for v field in Fourier space
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!					see http://www.2decomp.org/ for more information
	! OUTPUT
	!
	! .. Arrays ..
	!  u 				= updated u field
	!  v 				= updated v field
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
	!  temp             = temporary storage variable
	!  uhattemp         = temporary storage variable
	!  vhattemp         = temporary storage variable
	!  uhatstore        = temporary storage variable
	!  vhatstore        = temporary storage variable
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
	IMPLICIT NONE					 
	INCLUDE 'mpif.h'
	! Declare variables
	TYPE(DECOMP_INFO), INTENT(IN)						    :: decomp
	REAL(KIND=8), INTENT(IN)                                :: dt,a,b,Du,Dv,modescalereal
	COMPLEX(KIND=8), INTENT(IN)                             :: coef
	COMPLEX(kind=8), DIMENSION(decomp%zst(1):decomp%zen(1)), INTENT(IN) :: kx
	COMPLEX(kind=8), DIMENSION(decomp%zst(2):decomp%zen(2)), INTENT(IN) :: ky
	COMPLEX(kind=8), DIMENSION(decomp%zst(3):decomp%zen(3)), INTENT(IN) :: kz	
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   							   decomp%xst(2):decomp%xen(2),&
   							   decomp%xst(3):decomp%xen(3)),&
   											   INTENT(INOUT) :: u,v
	COMPLEX(KIND=8), DIMENSION(decomp%zst(1):decomp%zen(1),&
                               decomp%zst(2):decomp%zen(2),&
                               decomp%zst(3):decomp%zen(3)),&
   											INTENT(OUT)      :: uhat,vhat
	INTEGER(kind=4)										     :: i,j,k
	COMPLEX(kind=8)                   :: uhattemp,vhattemp,uhatstore,vhatstore
	CALL decomp_2d_fft_3d(u,uhat,DECOMP_2D_FFT_FORWARD)
	CALL decomp_2d_fft_3d(v,vhat,DECOMP_2D_FFT_FORWARD)
	IF ((decomp%zst(3).eq.1).and.(decomp%zst(2).eq.1).and.(decomp%zst(1).eq.1)) THEN
	 uhattemp=uhat(1,1,1)
	 vhattemp=vhat(1,1,1)
	END IF
 	DO k=decomp%zst(3),decomp%zen(3); DO j=decomp%zst(2),decomp%zen(2); DO i=decomp%zst(1),decomp%zen(1)
 	 uhatstore=uhat(i,j,k)
 	 vhatstore=vhat(i,j,k)
 	 vhat(i,j,k)=vhatstore*exp(coef*dt*Dv*(kx(i)**2+ky(j)**2+kz(k)**2)) &
 	 	            +b*uhatstore*(exp(coef*dt*Dv*(kx(i)**2+ky(j)**2+kz(k)**2)) &
 	 		                     -exp(coef*dt*(Du*(kx(i)**2+ky(j)**2+kz(k)**2)-b-1.0d0) )) &
 	 	                         /(b+1.0d0-(Du-Dv)*(kx(i)**2+ky(j)**2+kz(k)**2))
 	 uhat(i,j,k)=uhatstore*exp( coef*dt*(-b-1.0d0+Du*( kx(i)**2 + ky(j)**2 + kz(k)**2 )) )                
 	END DO; END DO; END DO
	IF ((decomp%zst(3).eq.1).and.(decomp%zst(2).eq.1).and.(decomp%zst(1).eq.1)) THEN
	 uhat(1,1,1)=uhattemp*exp(-coef*dt*(b+1.0d0))+&
	               (a/modescalereal)*(1.0d0-exp(-coef*dt*(b+1.0d0)))/(b+1.0d0)
	 vhat(1,1,1)=vhattemp+b*( (exp(-(b+1)*coef*dt)-1.0d0)* &
			                (-uhattemp+(a/modescalereal)/(b+1.0d0))/(b+1.0d0) &
			                    +(a/modescalereal)*coef*dt/(b+1.0d0))
	END IF
	! don't forget to scale
	CALL decomp_2d_fft_3d(uhat,u,DECOMP_2D_FFT_BACKWARD)
	CALL decomp_2d_fft_3d(vhat,v,DECOMP_2D_FFT_BACKWARD)
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	 u(i,j,k)=u(i,j,k)*modescalereal
	 v(i,j,k)=v(i,j,k)*modescalereal
	END DO; END DO; END DO
	END SUBROUTINE linear
