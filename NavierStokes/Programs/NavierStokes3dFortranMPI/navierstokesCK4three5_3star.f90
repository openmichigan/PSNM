PROGRAM main
	USE decomp_2d
	USE decomp_2d_fft
	USE decomp_2d_io
		
	IMPLICIT NONE
	INCLUDE "mpif.h"	
	!Grid & Stepping
   	INTEGER(kind=4), PARAMETER 		:: Nx=512, Ny=512, Nz=512		
   	INTEGER(kind=4), PARAMETER 		:: Lx=1, Ly=1, Lz=1
	INTEGER(kind=4), PARAMETER		:: Nt=20,plotgap=10
	INTEGER(kind=4)				:: numplots=Nt/plotgap
	!choose initial conditions
	logical					::readinput=.false.,&
						  TGvortex=.true.,&
						  exactsoln=.false.
	!runime parameters
	INTEGER(kind=4), PARAMETER 	  	:: savedata1=0
	REAL(kind=8), PARAMETER			:: dt=0.005d0
	REAL(kind=8), PARAMETER			:: Re=1600.0d0 		
	REAL(kind=8)				:: theta=0.0d0
	!additional parameters & constants
	REAL(kind=8), PARAMETER	&
		::  pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8), PARAMETER		:: dx=2.0d0*pi*Lx/REAL(Nx,kind(0d0))
	REAL(kind=8), PARAMETER		:: dy=2.0d0*pi*Ly/REAL(Ny,kind(0d0))
	REAL(kind=8), PARAMETER		:: dz=2.0d0*pi*Lz/REAL(Nz,kind(0d0))
	REAL(kind=8), PARAMETER		::	ReInv=1.0d0/REAL(Re,kind(0d0))
	REAL(kind=8), PARAMETER		:: dtInv=1.0d0/REAL(dt,kind(0d0)) 
	!Carpenter-Kennedy coefficients, Notation is as in Algorithm 3 of Kectheson(2010)
   	!Butcher table coefficents
	REAL(kind=8), DIMENSION(6)	:: beta,gamma1,gamma2,gamma3
        REAL(kind=8), DIMENSION(7)      :: delta
	!computational arrays/scalars	
	REAL(kind=8)						::  scalemodes,mu,factor,t=0.0d0
	COMPLEX(kind=8), DIMENSION(:), ALLOCATABLE		:: kx,ky,kz 				
	REAL(kind=8), DIMENSION(:), ALLOCATABLE			:: x,y,z,time
	REAL(kind=8), DIMENSION(:), ALLOCATABLE			:: KE,Enstrophy,&
							            KEdissipationRate
	REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE		:: u,v,w,&
							            ux,uy,uz,&
							            vx,vy,vz,&
							            wx,wy,wz,&
							            omega,nonlin,temp
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE		:: uhat,vhat,what,&
								   omegahat,temp_c,temp_c1,&
								   nonlinhatuh,nonlinhatvh,&
								   nonlinhatwh,phat
	!exact solution
	REAL(kind=8), PARAMETER					:: sl=1.0d0,sk=1.0d0,sm=1.0d0,&
								   lamlkm=sqrt(sl**2+sk**2+sm**2)
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE		:: utemp,vtemp,wtemp,utemp2,vtemp2,wtemp2
	REAL(kind=8)						:: temp3
	!MPI and 2DECOMP variables 
	TYPE(DECOMP_INFO)					:: decomp,sp
	INTEGER(kind=4)						:: myid,numprocs,ierr
	INTEGER(kind=MPI_OFFSET_KIND)				:: filesize, disp
	INTEGER(kind=4)						:: p_row=0, p_col=0
	REAL(kind=8), DIMENSION(:), ALLOCATABLE			:: temp1,temp2
	!variables used for saving data and timing
	INTEGER(kind=4)						::  i,j,k,n,nn,kk,allocatestatus
	INTEGER(kind=4)						::  count, iol, ind, plotInt=0
	CHARACTER*100			 			:: name_config
	INTEGER(kind=4)						::	start, finish, count_rate
   	!initialisation of 2DECOMP&FFT and MPI	
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
	CALL decomp_2d_init(Nx,Ny,Nz,p_row,p_col)
	CALL decomp_info_init(Nx,Ny,Nz,decomp)
   	CALL decomp_info_init(Nx/2+1,Ny,Nz,sp)
	CALL decomp_2d_fft_init
	
	IF (myid.eq.0) THEN
		PRINT *,'Solve 3D incompressible NS equations with Fourier pseudospectral methods'
		PRINT *,'and fourth order five stage three storage location embedded Runge Kutta method.'
	   	PRINT *,'Grid:',Nx,'X',Ny,'X',Nz
		PRINT *,'dt:',dt
		PRINT *,'Re:',Re
		PRINT *,'Nt:',Nt
      		PRINT *,'plotgap:',plotgap
	END IF		
	ALLOCATE(time(1:Nt+1),kx(1:Nx),ky(1:Ny),kz(1:Nz),x(1:Nx),y(1:Ny),z(1:Nz),temp1(1:9),temp2(1:9),&
	KE(1:Nt+1),Enstrophy(1:Nt+1),KEdissipationRate(1:Nt+1),&
	u(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	v(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	w(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	ux(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	uy(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	uz(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	vx(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	vy(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	vz(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	wx(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	wy(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	wz(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	omega(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	nonlin(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	temp(decomp%xst(1):decomp%xen(1),decomp%xst(2):decomp%xen(2),decomp%xst(3):decomp%xen(3)),&
	uhat(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	vhat(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	what(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	omegahat(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	temp_c(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	temp_c1(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	nonlinhatuh(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	nonlinhatvh(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	nonlinhatwh(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	phat(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	utemp(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	vtemp(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	wtemp(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	utemp2(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	vtemp2(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	wtemp2(sp%zst(1):sp%zen(1),sp%zst(2):sp%zen(2),sp%zst(3):sp%zen(3)),&
	stat=allocatestatus)
	IF (allocatestatus .ne. 0) STOP 
	IF (myid.eq.0) THEN
	 	PRINT *,'allocated space'
	END IF	
	beta=(/0.0, 0.075152045700771, 0.211361016946069, 1.100713347634329, 0.728537814675568, 0.393172889823198/)
   	gamma1=(/0.0000000000000000, 0.0, -0.497531095840104, 1.010070514199942, -3.196559004608766, 1.717835630267259/)	
	gamma2=(/0.0, 1.0, 1.384996869124138, 3.878155713328178, -2.324512951813145, -0.514633322274467/)
	gamma3=(/0.0, 0.0, 0.0, 0.0, 1.642598936063715, 0.188295940828347/)
	delta=(/1.0, 0.081252332929194, -1.083849060586449, -1.096110881845602, 2.859440022030827, -0.655568367959557, -0.194421504490852/)
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
	IF (myid.eq.0) THEN
		PRINT *,'Setup grid and fourier frequencies'
	END IF	
	scalemodes=1.0d0/REAL(Nx*Ny*Nz,kind(0d0)); time(1)=0.0d0
	
   	IF(readInput.eqv..true.) THEN
		name_config='./data/restartData.dat'
		OPEN(unit=11,file=name_config,status='old')
		REWIND(11)
		READ(11,*) time(1),plotInt
		CLOSE(11)
      		name_config='./data/u'
      		CALL getName(name_config,plotInt)
      		CALL decomp_2d_read_one(1,u,name_config)
      		name_config='./data/v'
      		CALL getName(name_config,plotInt)
      		CALL decomp_2d_read_one(1,v,name_config)
      		name_config='./data/w'
      		CALL getName(name_config,plotInt)
      		CALL decomp_2d_read_one(1,w,name_config)
      		IF(myid.eq.0) THEN
         		PRINT*,'Read in values'
      		END IF
   		ELSE If(TGvortex.eqv..true.) THEN
	   		!initial conditions for Taylor-Green vortex
	   		factor=2.0/sqrt(3.0)
      			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	      			u(i,j,k)=factor*sin(theta+2.0*pi/3.0)*sin(x(i))*cos(y(j))*cos(z(k))
      			END DO; END DO; END DO
      			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	      			v(i,j,k)=factor*sin(theta-2.0*pi/3.0)*cos(x(i))*sin(y(j))*cos(z(k))
      			END DO ; END DO ; END DO
      			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	      			w(i,j,k)=factor*sin(theta)*cos(x(i))*cos(y(j))*sin(z(k))
      			END DO ; END DO ; END DO
      			IF(myid.eq.0) THEN
         			PRINT*,'Taylor-Green Vortex initial conditions'
      			END IF
   		ELSE IF(exactsoln.eqv..true.) THEN
			!special exact sol'n usted for testing
			factor=sqrt(3.0d0)
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				u(i,j,k)=-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
							+sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
			END DO; END DO; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				v(i,j,k)=0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
							-cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(1)/Re)
			END DO ; END DO ; END DO
			DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				w(i,j,k)=cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(1)/Re)
			END DO ; END DO ; END DO
      			IF(myid.eq.0) THEN
         			PRINT*,'Exact Solution initial conditions'
      			END IF
		ELSE
      			IF(myid.eq.0) THEN
         		PRINT*,'--------------------------------------------------------------------'
		   	PRINT*,'Set initial conditions to either readinput, TGvortex, or exactsoln.'
         		PRINT*,'--------------------------------------------------------------------'
         		goto 111
      		END IF
	END IF
	
	CALL decomp_2d_fft_3d(u,uhat)
	CALL decomp_2d_fft_3d(v,vhat)
	CALL decomp_2d_fft_3d(w,what)
	
	! derivative of u with respect to x, y, and z 
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,ux)	
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,uy)	
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,uz)	
    
	! derivative of v with respect to x, y, and z 
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,vx)		
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,vy)	
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,vz)		

	! derivative of w with respect to x, y, and z 
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wx)		
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wy)		
	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
		temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
	END DO ; END DO ; END DO
	CALL decomp_2d_fft_3d(temp_c,wz)	
  
   	!calculate vorticity
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		temp(i,j,k)=wy(i,j,k)-vz(i,j,k)
	END DO ; END DO ; END DO

	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		omega(i,j,k)=temp(i,j,k)*temp(i,j,k)
	END DO ; END DO ; END DO

	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		temp(i,j,k)=uz(i,j,k)-wx(i,j,k)
	END DO ; END DO ; END DO

	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		omega(i,j,k)=temp(i,j,k)*temp(i,j,k)+omega(i,j,k)
	END DO ; END DO ; END DO

	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		temp(i,j,k)=vx(i,j,k)-uy(i,j,k)
	END DO ; END DO ; END DO
	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
		omega(i,j,k)=temp(i,j,k)*temp(i,j,k)+omega(i,j,k)
	END DO ; END DO ; END DO

	!Calculate Kinetic Energy
	temp1(1) = sum(u*u)
	temp1(2) = sum(v*v)
	temp1(3) = sum(w*w)
	CALL MPI_ALLREDUCE(temp1(1:3),temp2(1:3),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	temp3=temp2(1)+temp2(2)+temp2(3)
	KE(1)=0.5d0*dx*dy*dz*temp3

	!Calculate Enstrophy
   	temp1(1)=sum(omega)
	CALL MPI_ALLREDUCE(temp1(1),temp3,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
   	Enstrophy(1)=0.5d0*dx*dy*dz*temp3
	
   	temp1=0.0;temp2=0.0
	!Calculate Kinetic Energy Dissipation Rate
	temp1(1)=sum(ux*ux)
	temp1(2)=sum(uy*uy)
	temp1(3)=sum(uz*uz)
	temp1(4)=sum(vx*vx)
	temp1(5)=sum(vy*vy)
	temp1(6)=sum(vz*vz)
	temp1(7)=sum(wx*wx)
	temp1(8)=sum(wy*wy)
	temp1(9)=sum(wz*wz)
	CALL MPI_ALLREDUCE(temp1,temp2,9,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	temp3=temp2(1)+temp2(2)+temp2(3)+temp2(4)+temp2(5)+temp2(6)+temp2(7)+temp2(8)+temp2(9)
	KEdissipationRate(1)=2.0d0*dx*dy*dz*temp3

   	IF(myid.eq.0) THEN
      		PRINT *,'time=',time(1)
	   	PRINT*,'KE=',KE(1)
	   	PRINT*,'Enstrophy=',Enstrophy(1)
	   	PRINT*,'KEdissipationRate=',KEdissipationRate(1)
      		PRINT*,'---------------------------------'
		PRINT*,'starting timestepping and timer'
     		PRINT*,'---------------------------------'
   	END IF 
	nonlinhatuh(:,:,:)=0.0d0; nonlinhatvh(:,:,:)=0.0d0; nonlinhatwh(:,:,:)=0.0d0;
   	phat(:,:,:)=0.0d0

	CALL system_clock(start,count_rate)   
	DO n=1,numplots
      		DO nn=1,plotgap
         		t=t+dt
			DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
         			utemp(i,j,k)=0.0d0
				vtemp(i,j,k)=0.0d0
				wtemp(i,j,j)=0.0d0
				utemp2(i,j,k)=uhat(i,j,k)
				vtemp2(i,j,k)=vhat(i,j,k)
				wtemp2(i,j,k)=what(i,j,k)
			END DO ; END DO ; END DO
		   	DO kk=2,6
            			DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
               				utemp2(i,j,k)=utemp2(i,j,k)+delta(kk-1)*uhat(i,j,k)
               				vtemp2(i,j,k)=vtemp2(i,j,k)+delta(kk-1)*vhat(i,j,k)
               				wtemp2(i,j,k)=wtemp2(i,j,k)+delta(kk-1)*what(i,j,k)
            			END DO ; END DO ; END DO

			   	!nonlinuhat
			   	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   	temp(i,j,k)=(u(i,j,k)*ux(i,j,k)+v(i,j,k)*uy(i,j,k)+w(i,j,k)*uz(i,j,k))
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp,nonlinhatuh)
				
			   	!nonlinvhat
			   	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   	temp(i,j,k)=u(i,j,k)*vx(i,j,k)+v(i,j,k)*vy(i,j,k)+w(i,j,k)*vz(i,j,k)
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp,nonlinhatvh)
			    
			   	!nonlinwhat
			   	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   	temp(i,j,k)=u(i,j,k)*wx(i,j,k)+v(i,j,k)*wy(i,j,k)+w(i,j,k)*wz(i,j,k)
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp,nonlinhatwh)

			   	!phat
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	phat(i,j,k)=-1.0*(kx(i)*nonlinhatuh(i,j,k)+ky(j)*nonlinhatvh(i,j,k)+kz(k)*nonlinhatwh(i,j,k))&
                           		/(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)+REAL(0.1d0**13,kind(0d0)))
			   	END DO ; END DO ; END DO	

			   	!uhat,vhat,what
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	factor=(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)+REAL(0.1d0**13,kind(0d0)))
				   	uhat(i,j,k)=gamma1(kk)*uhat(i,j,k)+gamma2(kk)*utemp(i,j,k)+gamma3(kk)*utemp2(i,j,k)+beta(kk)*dt&
                           		*(-nonlinhatuh(i,j,k)-kx(i)*phat(i,j,k)+ReInv*factor*uhat(i,j,k));
				   	vhat(i,j,k)=gamma1(kk)*vhat(i,j,k)+gamma2(kk)*vtemp(i,j,k)+gamma3(kk)*vtemp2(i,j,k)+beta(kk)*dt&
                           		*(-nonlinhatvh(i,j,k)-ky(j)*phat(i,j,k)+ReInv*factor*vhat(i,j,k));
				   	what(i,j,k)=gamma1(kk)*what(i,j,k)+gamma2(kk)*wtemp(i,j,k)+gamma3(kk)*wtemp2(i,j,k)+beta(kk)*dt&
                           		*(-nonlinhatwh(i,j,k)-kz(k)*phat(i,j,k)+ReInv*factor*what(i,j,k));
			   	END DO ; END DO ; END DO
			
			   	!ifft to get u,v,w and derivatives
			
			
			   	! derivative of u with respect to x, y, and z 
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	temp_c(i,j,k)=uhat(i,j,k)*kx(i)*scalemodes
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp_c,ux)	
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	temp_c(i,j,k)=uhat(i,j,k)*ky(j)*scalemodes
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp_c,uy)	
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	temp_c(i,j,k)=uhat(i,j,k)*kz(k)*scalemodes
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp_c,uz)	

			   	! derivative of v with respect to x, y, and z 
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	temp_c(i,j,k)=vhat(i,j,k)*kx(i)*scalemodes
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp_c,vx)		
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	temp_c(i,j,k)=vhat(i,j,k)*ky(j)*scalemodes
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp_c,vy)	
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	temp_c(i,j,k)=vhat(i,j,k)*kz(k)*scalemodes
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp_c,vz)		

			   	! derivative of w with respect to x, y, and z 
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	temp_c(i,j,k)=what(i,j,k)*kx(i)*scalemodes
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp_c,wx)		
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	temp_c(i,j,k)=what(i,j,k)*ky(j)*scalemodes
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp_c,wy)		
			   	DO k=sp%zst(3),sp%zen(3) ; DO j=sp%zst(2),sp%zen(2) ; DO i=sp%zst(1),sp%zen(1)
				   	temp_c(i,j,k)=what(i,j,k)*kz(k)*scalemodes
			   	END DO ; END DO ; END DO
			   	CALL decomp_2d_fft_3d(temp_c,wz)	
				
			   	CALL decomp_2d_fft_3d(uhat,u)	
			   	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)	
				   	u(i,j,k)=u(i,j,k)*scalemodes
			   	END DO; END DO; END DO
			   	CALL decomp_2d_fft_3d(vhat,v)
			   	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   	v(i,j,k)=v(i,j,k)*scalemodes
			   	END DO; END DO; END DO
			   	CALL decomp_2d_fft_3d(what,w)
			   	DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
				   	w(i,j,k)=w(i,j,k)*scalemodes
			   	END DO; END DO; END DO
			END DO

		END DO
		!calculate vorticity
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			temp(i,j,k)=REAL(wy(i,j,k)-vz(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			omega(i,j,k)=REAL(temp(i,j,k)*temp(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			temp(i,j,k)=REAL(uz(i,j,k)-wx(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			omega(i,j,k)=REAL(temp(i,j,k)*temp(i,j,k)+omega(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			temp(i,j,k)=REAL(vx(i,j,k)-uy(i,j,k),KIND=8)
		END DO ; END DO ; END DO

		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			omega(i,j,k)=REAL(temp(i,j,k)*temp(i,j,k)+omega(i,j,k),KIND=8)
		END DO ; END DO ; END DO

	   	!Calculate Kinetic Energy
	   	temp1(1) = sum(u*u)
	   	temp1(2) = sum(v*v)
	   	temp1(3) = sum(w*w)
	   	CALL MPI_ALLREDUCE(temp1(1:3),temp2(1:3),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	   	temp3=temp2(1)+temp2(2)+temp2(3)
	   	KE(n+1)=0.5d0*dx*dy*dz*temp3

	   	!Calculate Enstrophy
      		temp1(1)=sum(omega)
	   	CALL MPI_ALLREDUCE(temp1(1),temp3,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      		Enstrophy(n+1)=0.5d0*dx*dy*dz*temp3

	   	temp1=0.0;temp2=0.0
	   	!Calculate Kinetic Energy Dissipation Rate
	   	temp1(1)=sum(ux*ux)
	   	temp1(2)=sum(uy*uy)
	   	temp1(3)=sum(uz*uz)
	   	temp1(4)=sum(vx*vx)
	   	temp1(5)=sum(vy*vy)
	   	temp1(6)=sum(vz*vz)
	   	temp1(7)=sum(wx*wx)
	   	temp1(8)=sum(wy*wy)
	   	temp1(9)=sum(wz*wz)
	   	CALL MPI_ALLREDUCE(temp1,temp2,9,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
	   	temp3=temp2(1)+temp2(2)+temp2(3)+temp2(4)+temp2(5)+temp2(6)+temp2(7)+temp2(8)+temp2(9)
	   	KEdissipationRate(n+1)=2.0d0*dx*dy*dz*temp3

      		time(n+1)=time(1)+plotgap*dt*n
     		IF(myid.eq.0) THEN
         		PRINT*,'time=',time(n+1)
	      		PRINT*,'KE=',KE(n+1)
	      		PRINT*,'Enstrophy=',Enstrophy(n+1)
	      		PRINT*,'KEdissipationRate=',KEdissipationRate(n+1)
         		PRINT*,'---------------------------------'
      		END IF
      		IF(savedata1==1) THEN
		   	plotInt=plotInt+1
         		name_config='./data/u'
         		CALL savedata(Nx,Ny,Nz,plotInt,name_config,u,decomp)
         		name_config='./data/v'
         		CALL savedata(Nx,Ny,Nz,plotInt,name_config,v,decomp)
         		name_config='./data/w'
         		CALL savedata(Nx,Ny,Nz,plotInt,name_config,w,decomp)
			name_config='./data/omega'
			CALL savedata(Nx,Ny,Nz,plotInt,name_config,omega,decomp)
      		END IF
		   	IF (savedata1==1) THEN
      		IF (myid.eq.0) THEN
         		name_config='./data/KineticEnergy.dat' 
	      		OPEN(UNIT=13, FILE=name_config,position='append') 
	         	write(13,*) KE(n+1)
         		CLOSE(13)

         		name_config='./data/Enstrophy.dat' 
	      		OPEN(UNIT=13, FILE=name_config,position='append') 
	         	write(13,*)  Enstrophy(n+1)
         		CLOSE(13)

         		name_config='./data/KEdissipationRate.dat' 
	      		OPEN(UNIT=13, FILE=name_config,position='append') 
	         	write(13,*) KEdissipationRate(n+1)
         		CLOSE(13)

         		name_config='./data/tdata.dat' 
	      		OPEN(UNIT=13, FILE=name_config,position='append') 
	         	write(13,*) time(n+1)
         		CLOSE(13)

			!time and plot number needed for restart
         		name_config='./data/restartData.dat'
	      		OPEN(UNIT=13, FILE=name_config, status='unknown') 
			REWIND(13)
	      		write(13,*) time(Nt+1),plotInt
         		CLOSE(13)
      		END IF
   	END IF


	END DO

	
    	CALL system_clock(finish,count_rate)
    	IF (myid.eq.0) then
    		PRINT *, 'Program took', REAL(finish-start)/REAL(count_rate), 'for main timestepping loop'
    	END IF	
    
	IF(exactsoln.eqv..true.) THEN
		!Calculate error in final numerical solution
      		factor=sqrt(3.0d0)
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			temp(i,j,k)=u(i,j,k) -&
							(-0.5*( factor*cos(x(i))*sin(y(j))*sin(z(k))&
							+sin(x(i))*cos(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
		END DO; END DO; END DO
      		temp1(1) = maxval(abs(temp));
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			temp(i,j,k)=v(i,j,k) -&
						(0.5*(  factor*sin(x(i))*cos(y(j))*sin(z(k))&
							-cos(x(i))*sin(y(j))*cos(z(k)) )*exp(-(factor**2)*time(Nt+1)/Re))
		END DO ; END DO ; END DO
      		temp1(2) = maxval(abs(temp));
		DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
			temp(i,j,k)=w(i,j,k)-&
						(cos(x(i))*cos(y(j))*sin(z(k))*exp(-(factor**2)*time(Nt+1)/Re))
		END DO ; END DO ; END DO
		temp1(3) = maxval(abs(temp));
		CALL MPI_ALLREDUCE(temp1(1:3),temp2(1:3),3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
		temp3=temp2(1)+temp2(2)+temp2(3)
		IF (myid.eq.0) THEN
			PRINT*,'The error at time',time(Nt+1),'is=',temp3
		END IF
	END IF

   111 continue
	DEALLOCATE(time,kx,ky,kz,x,y,z,temp1,temp2,&
			  KE,Enstrophy,KEdissipationRate,&
			  u,v,w,&
			  ux,uy,uz,&
			  vx,vy,vz,&
			  wx,wy,wz,&
			  omega,nonlin,temp,&
			  uhat,vhat,what,&
			  omegahat,temp_c,temp_c1,&
			  nonlinhatuh,nonlinhatvh,&
			  nonlinhatwh,phat,&
			  utemp,vtemp,wtemp,&
			  utemp2,vtemp2,wtemp2,&
			  stat=allocatestatus)
	IF (allocatestatus .ne. 0) THEN
                PRINT *,'STOP'
                STOP
   	END IF
   	CALL decomp_2d_fft_finalize
   	CALL decomp_2d_finalize 
	IF (myid.eq.0) THEN
		PRINT *,'Program execution complete'
	END IF
	CALL MPI_FINALIZE(ierr)
END PROGRAM main

   SUBROUTINE getName(name_config,plotInt)
      implicit none
      CHARACTER*100,intent(inout)   :: name_config
      INTEGER(kind=4)               :: plotInt
      CHARACTER*100                 :: number_file
      INTEGER(kind=4)               :: ind

      ind = index(name_config,' ') - 1
      WRITE(number_file,'(i0)') 10000000+plotInt
      number_file = name_config(1:ind)//number_file
      ind = index(number_file,' ') - 1
      name_config = number_file(1:ind)//'.datbin'
   END SUBROUTINE getname
