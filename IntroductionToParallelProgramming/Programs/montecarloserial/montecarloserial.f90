




	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program use a monte carlo method to calculate pi
	!
	! .. Parameters ..
	!  npts				= total number of Monte Carlo points
	!  xmin				= lower bound for integration region
	!  xmax             = upper bound for integration region
	! .. Scalars ..
	!  i                = loop counter
	!  f				= average value from summation
	!  sum              = total sum
	!  randnum          = random number generated from (0,1) uniform 
	!                     distribution
	!  x                = current Monte Carlo location
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
	! None
 	PROGRAM monte_carlo
    IMPLICIT NONE

    INTEGER(kind=8), PARAMETER 	    :: npts = 1e10
    REAL(kind=8), PARAMETER 		:: xmin=0.0d0,xmax=1.0d0
    INTEGER(kind=8) 				:: i
    REAL(kind=8) 					:: f,sum, randnum,x

    DO i=1,npts
      CALL random_number(randnum)
      x = (xmax-xmin)*randnum + xmin
      sum = sum + 4.0d0/(1.0d0 + x**2)
    END DO
    f = sum/npts
    PRINT*,'PI calculated with ',npts,' points = ',f

    STOP
    END
