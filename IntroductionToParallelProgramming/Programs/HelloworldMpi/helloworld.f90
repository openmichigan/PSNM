	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program uses MPI to print hello world from all available
	! processes
	!
	! .. Parameters ..
	!
	! .. Scalars ..
	!  myid				= process id
	!  numprocs			= total number of MPI processes
	!  ierr				= error code
	!
	! .. Arrays ..
	!
	! .. Vectors ..
	!
	! REFERENCES
	! http:// en.wikipedia.org/wiki/OpenMP
	!
	! ACKNOWLEDGEMENTS
	! The program below was modified from one available at the internet
	! address in the references. This internet address was last checked
	! on 30 December 2011
	!
	! ACCURACY
	!		
	! ERROR INDICATORS AND WARNINGS
	!
	! FURTHER COMMENTS
	!
	!--------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	! MPI library
	PROGRAM hello90
	USE MPI
	IMPLICIT NONE
	INTEGER(kind=4) :: myid, numprocs, ierr
	
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
	
	PRINT*, 'Hello World from process', myid
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	IF ( myid == 0 ) THEN
		PRINT*, 'There are ', numprocs, ' MPI processes'
	END IF
	CALL MPI_FINALIZE(ierr)		
	END PROGRAM
