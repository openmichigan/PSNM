	PROGRAM BovCreate
	!--------------------------------------------------------------------------------
	! .. Purpose ..
	! 	BovCreate is a postprocessing program which creates header files for VisIt
	!	It uses the INPUTFILE and assumes that the filenames in the program are
	!	consistent with those in the current file.
	!	
	! .. PARAMETERS .. INITIALIZED IN INPUTFILE
	!	time		= start time of the simulation
	!	Nx			= power of two, number of modes in the x direction
	!	Ny			= power of two, number of modes in the y direction
	!	Nz			= power of two, number of modes in the z direction
	!	Nt			= the number of timesteps
	!	plotgap		= the number of timesteps to take before plotting
	!	Lx			= definition of the periodic domain of computation in x direction
	!	Ly			= definition of the periodic domain of computation in y direction
	!	Lz			= definition of the periodic domain of computation in z direction
	!	Es			= focusing or defocusing
	!	Dt			= the time step
	!
	!	REFERENCES
	! 
	!	ACCURACY
	!
	!	ERROR INDICATORS AND WARNINGS
	!
	!	FURTHER COMMENTS
	!------------------------------------------------------------------------------------
	! 	EXTERNAL ROUTINES REQUIRED
	IMPLICIT NONE
	! .. Scalar Arguments ..
	INTEGER(KIND=4)			::  Nx, Ny, Nz, Nt, plotgap 
	REAL(KIND=8)			::  Lx, Ly, Lz,  DT, time, Es
	! .. Local scalars ..
	INTEGER(KIND=4)			::  stat,plotnum,ind,n,numplots
	! .. Local Arrays ..
	CHARACTER*50			::  InputFileName, OutputFileName, OutputFileName2
	CHARACTER*10			::  number_file
	InputFileName='INPUTFILE'
	OPEN(unit=11,FILE=trim(InputFileName),status="OLD") 
	REWIND(11)
	READ(11,*) Nx, Ny, Nz, Nt, plotgap, Lx, Ly, Lz, Es, DT, time 
	CLOSE(11)	
	
	plotnum=1
	numplots=1+Nt/plotgap
	DO n=1,numplots
		OutputFileName = 'data/u'
		ind = index(OutputFileName,' ') - 1
		WRITE(number_file,'(i0)') 10000000+plotnum
		OutputFileName = OutputFileName(1:ind)//number_file
		ind = index(OutputFileName,' ') - 1
		OutputFileName = OutputFileName(1:ind)//'.bov'
		OutputFileName2='u'
		ind = index(OutputFileName2,' ') - 1
		OutputFileName2 = OutputFileName2(1:ind)//number_file
		ind = index(OutputFileName2,' ') - 1
		OutputFileName2 = OutputFileName2(1:ind)//'.datbin'
		OPEN(unit=11,FILE=trim(OutputFileName),status="UNKNOWN") 
		REWIND(11)
		WRITE(11,*) 'TIME: ',time
		WRITE(11,*) 'DATA_FILE: ',trim(OutputFileName2)
		WRITE(11,*) 'DATA_SIZE: ', Nx, Ny, Nz
		WRITE(11,*) 'DATA_FORMAT: DOUBLE'
		WRITE(11,*) 'VARIABLE: u'
		WRITE(11,*) 'DATA_ENDIAN: LITTLE'
		WRITE(11,*) 'CENTERING: ZONAL'
		WRITE(11,*) 'BRICK_ORIGIN:', -Nx/2, -Ny/2, -Nz/2
		WRITE(11,*) 'BRICK_SIZE:', Nx, Ny, Nz
		WRITE(11,*) 'DIVIDE_BRICK: true'
		WRITE(11,*) 'DATA_BRICKLETS:', Nx/2, Ny/2, Nz/2
		CLOSE(11)

		time=time+plotgap*DT
		plotnum=plotnum+1
	END DO	
	END PROGRAM BovCreate
