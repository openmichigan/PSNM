




	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program uses MPI to do a parallel monte carlo calculation of pi
	!
	! .. Parameters ..
	!  npts				= total number of Monte Carlo points
	!  xmin				= lower bound for integration region
	!  xmax             = upper bound for integration region
	! .. Scalars ..
	!  mynpts			= this processes number of Monte Carlo points
	!  myid				= process id
	!  nprocs			= total number of MPI processes
	!  ierr				= error code
	!  i                = loop counter
	!  f				= average value from summation
	!  sum              = total sum
	!  mysum            = sum on this process
	!  randnum          = random number generated from (0,1) uniform 
	!                     distribution
	!  x                = current Monte Carlo location
	!  start			= simulation start time
	!  finish			= simulation end time
	! .. Arrays ..
	!
	! .. Vectors ..
	!
	! REFERENCES
	! http://chpc.wustl.edu/mpi-fortran.html
	! Gropp, Lusk and Skjellum, "Using MPI" MIT press (1999)
	!
	! ACKNOWLEDGEMENTS
	! The program below was modified from one available at the internet
	! address in the references. This internet address was last checked
	! on 30 March 2012
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
    PROGRAM monte_carlo_mpi
    USE MPI
    IMPLICIT NONE

    INTEGER(kind=8), PARAMETER 	:: npts = 1e10
    REAL(kind=8), PARAMETER 	:: xmin=0.0d0,xmax=1.0d0
    INTEGER(kind=8) 			:: mynpts
    INTEGER(kind=4) 	        :: ierr, myid, nprocs
    INTEGER(kind=8) 	        :: i
    REAL(kind=8) 		        :: f,sum,mysum,randnum
    REAL(kind=8) 		        :: x, start, finish
    
    ! Initialize MPI
    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    start=MPI_WTIME()

	! Calculate the number of points each MPI process needs to generate
    IF (myid .eq. 0) THEN
      mynpts = npts - (nprocs-1)*(npts/nprocs)
    ELSE
      mynpts = npts/nprocs
    ENDIF
    
    ! set initial sum to zero
    mysum = 0.0d0
	! use loop on local process to generate portion of Monte Carlo integral
    DO i=1,mynpts
      CALL random_number(randnum)
      x = (xmax-xmin)*randnum + xmin
      mysum = mysum + 4.0d0/(1.0d0 + x**2)
    ENDDO

	! Do a reduction and sum the results from all processes
    CALL MPI_REDUCE(mysum,sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
          0,MPI_COMM_WORLD,ierr)
    finish=MPI_WTIME()

    ! Get one process to output the result and running time
    IF (myid .eq. 0) THEN
         f = sum/npts
         PRINT*,'PI calculated with ',npts,' points = ',f
         PRINT*,'Program took ', finish-start, ' for Time stepping'
    ENDIF

    CALL MPI_FINALIZE(ierr)

    STOP
    END PROGRAM
