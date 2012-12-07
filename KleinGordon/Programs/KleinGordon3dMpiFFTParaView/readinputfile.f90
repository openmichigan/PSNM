	SUBROUTINE readinputfile(Nx,Ny,Nz,Nt,plotgap,Lx,Ly,Lz, &
								Es,DT,starttime,myid,ierr)

	!--------------------------------------------------------------------------------
	!
	!
	!  PURPOSE 
	!
	!	Read inputfile intialize parameters, which are stocked in the Input File 
	!	
	! .. INPUT ..
	!	Nx				= number of modes in the x direction
	!	Ny				= number of modes in the y direction
	!	Nz				= number of modes in the z direction
	!	Nt				= the number of timesteps
	!	plotgap			= the number of timesteps to take before plotting
	!	myid			= number of MPI process
	!	ierr			= MPI error output variable
	!	Lx				= size of the periodic domain of computation in x direction
	!	Ly				= size of the periodic domain of computation in y direction
	!	Lz				= size of the periodic domain of computation in z direction
	!	DT				= the time step
	!	starttime		= initial time of computation
	!	InputFileName	= name of the Input File
	! REFERENCES
	!
	! ACCURACY
	!		
	! ERROR INDICATORS AND WARNINGS
	!
	! FURTHER COMMENTS
	!---------------------------------------------------------------------------------
	! EXTERNAL ROUTINES REQUIRED

	IMPLICIT NONE
	INCLUDE 'mpif.h'
	! .. Scalar Arguments ..
	INTEGER(KIND=4), INTENT(IN)		::  myid	
	INTEGER(KIND=4), INTENT(OUT)	::  Nx,Ny,Nz,Nt 
	INTEGER(KIND=4), INTENT(OUT)	::  plotgap, ierr 
	REAL(KIND=8), INTENT(OUT)		::  Lx, Ly, Lz, DT, starttime, Es
	! .. Local scalars ..
	INTEGER(KIND=4)					::  stat
	! .. Local Arrays ..
	CHARACTER*40					::  InputFileName
	INTEGER(KIND=4), DIMENSION(1:5)	::  intcomm
	REAL(KIND=8), DIMENSION(1:6)	::  dpcomm

	IF(myid.eq.0) THEN
		CALL GET_ENVIRONMENT_VARIABLE(NAME="inputfile",VALUE=InputFileName, STATUS=stat)
		IF(stat.NE.0) THEN
			PRINT*,"Set environment variable inputfile to the name of the" 
			PRINT*,"file where the simulation parameters are set"
			STOP
		END IF
		OPEN(unit=11,FILE=trim(InputFileName),status="OLD") 
		REWIND(11)
		READ(11,*) intcomm(1), intcomm(2), intcomm(3), intcomm(4), intcomm(5), &
			 	dpcomm(1), dpcomm(2), dpcomm(3), dpcomm(4), dpcomm(5), dpcomm(6)
		CLOSE(11)
		PRINT *,"NX ",intcomm(1)
		PRINT *,"NY ",intcomm(2)
		PRINT *,"NZ ",intcomm(3)
		PRINT *,"NT ",intcomm(4)
		PRINT *,"plotgap ",intcomm(5)
		PRINT *,"Lx ",dpcomm(1)
		PRINT *,"Ly ",dpcomm(2)
		PRINT *,"Lz ",dpcomm(3)
		PRINT *,"Es ",dpcomm(4)		
		PRINT *,"Dt ",dpcomm(5)
		PRINT *,"strart time ",dpcomm(6)
		PRINT *,"Read inputfile"
	END IF
	CALL MPI_BCAST(dpcomm,6,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(intcomm,5,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	Nx=intcomm(1)
	Ny=intcomm(2)
	Nz=intcomm(3)
	Nt=intcomm(4)
	plotgap=intcomm(5)
	Lx=dpcomm(1)
	Ly=dpcomm(2)
	Lz=dpcomm(3)
	Es=dpcomm(4)
	DT=dpcomm(5)
	starttime=dpcomm(6)
	
	END SUBROUTINE readinputfile
