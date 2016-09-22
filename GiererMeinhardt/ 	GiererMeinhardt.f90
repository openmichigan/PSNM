!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program solves Gierer Meinhardt equation in 3 dimensions
	! u_t + a(u_x + u_y + u_z) = d_1(u_{xx}+u_{yy}+u_{zz}) + p u^2/v - q_1 u + w_1
	! v_t + a(v_x + v_y + v_z) = d_2(v_{xx}+v_{yy}+v_{zz}) + p u^2 - q_2 u + w_2
	!
	! .. Parameters ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  Nz				= number of modes in z - power of 2 for FFT
	!  Nt				= number of timesteps to take
	!  Tmax				= maximum simulation time
	!  plotgap			= number of timesteps between plots
	!  pi = 3.14159265358979323846264338327950288419716939937510d0
	!  Lx				= width of box in x direction
	!  Ly				= width of box in y direction
	!  Lz				= width of box in z direction
	!  tol				=max error each timestep
	! .. Scalars ..
	!  pp				= order of the splitting method
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
	!  n				= loop counter for timesteps direction	
	!  allocatestatus	= error indicator during allocation
	!  start			= variable to record start time of program
	!  finish			= variable to record end time of program
	!  count_rate		= variable for clock count rate
	!  dt				= timestep
	!  modescalereal	= Number to scale after backward FFT 
	!  myid				= Process id
	!  ierr				= error code
	!  p_row			= number of rows for domain decomposition
	!  p_col			= number of columns for domain decomposition
	!  plotnum			= number of plot to save
	!  stat				= error indicator when reading inputfile
	! .. Arrays ..
	!  u 			= approximate solution u higher order
	!  v			= approximate solution v higher order
	!  ulow 			= approximate solution u lower order
	!  vlow				= approximate solution v lower order
	!  uhat 			= Fourier transform of u
	!  vhat				= Fourier transform of v
	!  savearray			= stores the realpart to save the data
	! .. Vectors ..
	!  kx				= fourier frequencies in x direction squared
	!  ky				= fourier frequencies in y direction squared
	!  kz				= fourier frequencies in z direction squared
	!  x				= x locations
	!  y				= y locations
	!  z				= z locations
	!  time				= times at which save data
	!  nameconfig		= array to store filename for data to be saved    		
	!  dpcomm
	!  intcomm
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!  sp				see http://www.2decomp.org/ for more information
	!
	! REFERENCES
	!
	! ACKNOWLEDGEMENTS
	! Modified from original version by M. Quell
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
		
		
	PROGRAM main
	USE decomp_2d
	USE decomp_2d_fft
	USE decomp_2d_io
!	USE MPI
	! Declare variables	
	IMPLICIT NONE	
	INCLUDE 'mpif.h'	
	INTEGER(kind=4)		::  Nx=512
	INTEGER(kind=4) 	::  Ny=512
	INTEGER(kind=4) 	::  Nz=1
	INTEGER(kind=4) 	::  Nt=100
	INTEGER(kind=4) 	::  plotnum=1000000000
	REAL(kind=8), PARAMETER	::&
		pi=3.14159265358979323846264338327950288419716939937510d0
	REAL(kind=8)		::  Lx=1.0,Ly=1.0,Lz=1.0
	REAL(kind=8)		::  dt=0.001
	REAL(kind=8)		::  tol1=4.0*0.1**12
	REAL(kind=8)		::  time=0.0
	REAL(kind=8)		::  ppp=0.0
	!equation specific
	REAL(kind=8)		::  a=0.03
	REAL(kind=8)		::  d1=0.001
	REAL(kind=8)		::  d2=50.0
	REAL(kind=8)		::  p=0.1
	REAL(kind=8)		::  q1=1.0
	REAL(kind=8)		::  q2=100.0
	REAL(kind=8)		::  w1=1.0
	REAL(kind=8)		::  w2=1.0
	COMPLEX(kind=8)		::  uold,vold, umean,vmean, utemp, vtemp, uhatmp
	INTEGER(kind=4) 	::  cnt=0
	INTEGER(kind=4) 	::  cntmax=100
	REAL(kind=8)		::  tol=0.1**12
	REAL(kind=8)		::  chg=1.0
	COMPLEX(kind=8),	DIMENSION(:), ALLOCATABLE	::  kx,ky
	REAL(kind=8),		DIMENSION(:), ALLOCATABLE	::  x,y 
	COMPLEX(kind=8), 	DIMENSION(:,:,:), ALLOCATABLE	::  uhigh,vhigh
	COMPLEX(kind=8), 	DIMENSION(:,:,:), ALLOCATABLE	::  ulow,vlow
	COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE	::  uhat,vhat
	REAL(kind=8), DIMENSION(:,:,:), ALLOCATABLE	:: realfield
	REAL(kind=8) 		::  modescalereal
	INTEGER(kind=4)		::  i,j,k,l,n,m,mm,ll,AllocateStatus,stat,nmax
	INTEGER(kind=4)		::  myid,numprocs,ierr
	TYPE(DECOMP_INFO)	::  decomp
	INTEGER(kind=4)		::  p_row=8, p_col=1
	INTEGER(kind=4)		::  start, finish,  starttot, finishtot, count_rate, ind
	CHARACTER*500		::  nameconfig
	CHARACTER*200		::  numberfile, InputFileName
	! splitting coeffiecents
	REAL(kind=8), DIMENSION(1:2)	::  aa1=(/0.5,0.5/)
	REAL(kind=8), DIMENSION(1:2)	::  bb1=(/1.0,0.0/)
	REAL(kind=8), DIMENSION(1:2)	::  aa2=(/0.0,1.0/)
	REAL(kind=8), DIMENSION(1:2)	::  bb2=(/1.0,0.0/)
	INTEGER(kind=4)		::  pp=2
	INTEGER(kind=4)		::  ss=2
  CHARACTER*500           ::      name
  name='./data/STRANG'

	!initialisation of MPI
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
	
	CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	IF(myid.eq.0) THEN
		PRINT *,'Nx:',Nx,'Ny:',Ny!,'Nz:',Nz
	END IF
	CALL decomp_2d_init(Nx,Ny,Nz,p_row,p_col)
	! get information about domain decomposition choosen
	CALL decomp_info_init(Nx,Ny,Nz,decomp) ! physical domain
	! initialise FFT library
	CALL decomp_2d_fft_init
	ALLOCATE(realfield(decomp%xst(1):decomp%xen(1),&
   				decomp%xst(2):decomp%xen(2),&
   				decomp%xst(3):decomp%xen(3)),&
   				uhigh(decomp%xst(1):decomp%xen(1),&
   				decomp%xst(2):decomp%xen(2),&
   				decomp%xst(3):decomp%xen(3)),&
        Vhigh(decomp%xst(1):decomp%xen(1),&
            decomp%xst(2):decomp%xen(2),&
            decomp%xst(3):decomp%xen(3)),&
	    ulow(decomp%xst(1):decomp%xen(1),&
			decomp%xst(2):decomp%xen(2),&
			decomp%xst(3):decomp%xen(3)),&
        vlow(decomp%xst(1):decomp%xen(1),&
            decomp%xst(2):decomp%xen(2),&
            decomp%xst(3):decomp%xen(3)),&
	    uhat(decomp%zst(1):decomp%zen(1),&
            decomp%zst(2):decomp%zen(2),&
            decomp%zst(3):decomp%zen(3)),&
        vhat(decomp%zst(1):decomp%zen(1),&
            decomp%zst(2):decomp%zen(2),&
            decomp%zst(3):decomp%zen(3)),&
		kx(decomp%zst(1):decomp%zen(1)), &
		ky(decomp%zst(2):decomp%zen(2)), &
		x(decomp%xst(1):decomp%xen(1)),&
		y(decomp%xst(2):decomp%xen(2)),&
		stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	modescalereal=1.0d0/REAL(Nx,KIND(0d0))
	modescalereal=modescalereal/REAL(Ny,KIND(0d0))
	modescalereal=modescalereal/REAL(Nz,KIND(0d0))
	! setup fourier frequencies and grid points
	DO i = 1,1+ Nx/2
		IF ((i.GE.decomp%zst(1)).AND.(i.LE.decomp%zen(1))) THEN
			kx(i)= cmplx(0.0d0,1.0d0)*REAL(i-1,kind(0d0))/Lx
		END IF
	END DO
	IF ((Nx/2 + 1 .GE.decomp%zst(1)).AND.(Nx/2 + 1 .LE.decomp%zen(1))) THEN
		kx( Nx/2 + 1 ) = 0.0d0 
	ENDIF
	DO i = Nx/2+2, Nx  
		IF ((i.GE.decomp%zst(1)).AND.(i.LE.decomp%zen(1))) THEN
			Kx( i) = cmplx(0.0d0,-1.0d0)*REAL(1-i+Nx,KIND(0d0))/Lx  
		ENDIF
	END DO      
	DO i=decomp%xst(1),decomp%xen(1)
		x(i)=(-1.0d0 + 2.0d0*REAL(i-1,kind(0d0))/REAL(Nx,kind(0d0)))*pi*Lx
	END DO

	DO j = 1,1+ Ny/2
		IF ((j.GE.decomp%zst(2)).AND.(j.LE.decomp%zen(2))) THEN
			ky(j)= cmplx(0.0d0,1.0d0)*REAL(j-1,kind(0d0))/Ly
		END IF
	END DO
	IF ((Ny/2 + 1 .GE.decomp%zst(2)).AND.(Ny/2 + 1 .LE.decomp%zen(2))) THEN
		ky( Ny/2 + 1 ) = 0.0d0 
	ENDIF
	DO j = Ny/2+2, Ny  
		IF ((j.GE.decomp%zst(2)).AND.(j.LE.decomp%zen(2))) THEN
			ky(j) = cmplx(0.0d0,-1.0d0)*REAL(1-j+Ny,KIND(0d0))/Ly  
		ENDIF
	END DO      
	DO j=decomp%xst(2),decomp%xen(2)
		y(j)=(-1.0d0 + 2.0d0*REAL(j-1,kind(0d0))/REAL(Ny,kind(0d0)))*pi*Ly
	END DO
	
	DO k=decomp%xst(3),decomp%xen(3)
		DO j=decomp%xst(2),decomp%xen(2)
			DO i=decomp%xst(1),decomp%xen(1)
				uhigh(i,j,k)=0.5+exp(-1.0-(x(i)**2+y(j)**2))
				vhigh(i,j,k)=0.1+exp(-1.0-(x(i)**2+y(j)**2))
			END DO
		END DO
	END DO
	plotnum=0
	! write out using 2DECOMP&FFT MPI-IO routines
	nameconfig=name
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//'u'
	WRITE(numberfile,'(i0)') plotnum+1000000000
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//numberfile
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//'.datbin'
	realfield=real(uhigh,kind(0d0))
	CALL decomp_2d_write_one(1,realfield,nameconfig)
	!same for v
	nameconfig=name
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//'v'
	WRITE(numberfile,'(i0)') plotnum+1000000000
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//numberfile
	ind=index(nameconfig,' ') -1
	nameconfig=nameconfig(1:ind)//'.datbin'
	realfield=real(vhigh,kind(0d0))
	CALL decomp_2d_write_one(1,realfield,nameconfig)
	IF (myid.eq.0) THEN
		PRINT *,'Got initial data, starting timestepping'
	END IF
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL system_clock(start,count_rate)	
		DO n=1,Nt	
		DO l=1,ss
			DO k=decomp%xst(3),decomp%xen(3)
			 DO j=decomp%xst(2),decomp%xen(2)
			  DO i=decomp%xst(1),decomp%xen(1)
						chg=1.0
						cnt=0
						uold=uhigh(i,j,k)
						vold=vhigh(i,j,k)
						DO WHILE ((chg>tol).and.(cnt<cntmax))
						cnt=cnt+1
							utemp=uhigh(i,j,k)
							vtemp=vhigh(i,j,k)
							umean=0.5*(uhigh(i,j,k)+uold)
							vmean=0.5*(vhigh(i,j,k)+vold)
							uhigh(i,j,k)=uold+dt*aa1(l)*(w1+p*umean**2/vmean)
							vhigh(i,j,k)=vold+dt*aa1(l)*(w2+p*umean**2)
							chg=abs((uhigh(i,j,k)-utemp))+abs((vhigh(i,j,k)-vtemp))	
						END DO
			END DO
			 END DO
			  END DO
	! solve linear part exactly in Fourier space
			CALL decomp_2d_fft_3d(uhigh,uhat,DECOMP_2D_FFT_FORWARD)
			CALL decomp_2d_fft_3d(vhigh,vhat,DECOMP_2D_FFT_FORWARD)
			DO k=decomp%zst(3),decomp%zen(3)
			 DO j=decomp%zst(2),decomp%zen(2)
			  DO i=decomp%zst(1),decomp%zen(1)
uhat(i,j,k)=exp(dt*bb1(l)*(-q1+d1*(kx(i)**2+ky(j)**2)-a*(kx(i)+ky(j))))*uhat(i,j,k)
vhat(i,j,k)=exp(dt*bb1(l)*(-q2+d2*(kx(i)**2+ky(j)**2)-a*(kx(i)+ky(j))))*vhat(i,j,k)
			END DO
			 END DO
			  END DO
			!dont forget to scale
			CALL decomp_2d_fft_3d(uhat,uhigh,DECOMP_2D_FFT_BACKWARD)
			CALL decomp_2d_fft_3d(vhat,vhigh,DECOMP_2D_FFT_BACKWARD)
			DO k=decomp%xst(3),decomp%xen(3)
			 DO j=decomp%xst(2),decomp%xen(2)
			  DO i=decomp%xst(1),decomp%xen(1)
				uhigh(i,j,k)=uhigh(i,j,k)*modescalereal
				vhigh(i,j,k)=vhigh(i,j,k)*modescalereal
			END DO
			 END DO
			  END DO

		END DO !l		
	if(plotnum.eq.n) then
		! write out using 2DECOMP&FFT MPI-IO routines
	  nameconfig=name
	  ind=index(nameconfig,' ') -1
	  nameconfig=nameconfig(1:ind)//'u'
	  WRITE(numberfile,'(i0)') plotnum+1000000000
	  ind=index(nameconfig,' ') -1
	  nameconfig=nameconfig(1:ind)//numberfile
	  ind=index(nameconfig,' ') -1
	  nameconfig=nameconfig(1:ind)//'.datbin'
	  realfield=real(uhigh,kind(0d0))
	  CALL decomp_2d_write_one(1,realfield,nameconfig)
	  !same for v
	  nameconfig=name
	  ind=index(nameconfig,' ') -1
	  nameconfig=nameconfig(1:ind)//'v'
	  WRITE(numberfile,'(i0)') plotnum+1000000000
	  ind=index(nameconfig,' ') -1
	  nameconfig=nameconfig(1:ind)//numberfile
	  ind=index(nameconfig,' ') -1
	  nameconfig=nameconfig(1:ind)//'.datbin'
	  realfield=real(vhigh,kind(0d0))
    CALL decomp_2d_write_one(1,realfield,nameconfig)
    plotnum=plotnum*2
	end if
END DO !Nt loop
	CALL system_clock(finish,count_rate)
		time=REAL(finish-start)/REAL(count_rate)
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)	
	! clean up 
  	CALL decomp_2d_fft_finalize
  	CALL decomp_2d_finalize
 	DEALLOCATE(realfield,uhigh,vhigh,ulow,vlow,&
		  uhat,vhat,&
   		kx,ky,x,y,&
   		stat=AllocateStatus)	
	IF (AllocateStatus .ne. 0) STOP 
	IF (myid.eq.0) THEN
	   	PRINT *,'Program execution complete'
	END IF
	CALL MPI_FINALIZE(ierr)		
	END PROGRAM main
