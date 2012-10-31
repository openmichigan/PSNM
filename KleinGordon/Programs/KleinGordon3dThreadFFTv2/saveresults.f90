	SUBROUTINE saveresults(Nt,plotgap,time,en,enstr,enkin,enpot)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine saves the energy and times stored during the
	! computation for the nonlinear Klein-Gordon equation
	!
	! INPUT
	!
	! .. Parameters ..
	!  Nt				= number of timesteps 
	!  plotgap			= number of timesteps between plots
	! .. Vectors ..
	!  time				= times at which save data
	!  en				= total energy	
	!  enstr			= strain energy
	!  enpot			= potential energy
	!  enkin			= kinetic energy
	!
	! OUTPUT
	!
	!
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  n				= loop counter 
	! .. Arrays ..
	! 	name_config		= array to hold the filename
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
	!--------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	IMPLICIT NONE					 
	! Declare variables
	INTEGER(kind=4), INTENT(IN)							:: plotgap,Nt
	REAL(KIND=8), DIMENSION(1:1+Nt/plotgap), INTENT(IN)	:: enpot, enkin
	REAL(KIND=8), DIMENSION(1:1+Nt/plotgap), INTENT(IN)	:: en,enstr,time
	INTEGER(kind=4)										:: j
	CHARACTER*100										:: name_config

	name_config = 'tdata.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,1+Nt/plotgap
		WRITE(11,*) time(j)
	END DO
	CLOSE(11)

	name_config = 'en.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,1+Nt/plotgap
		WRITE(11,*) en(j)
	END DO
	CLOSE(11)

	name_config = 'enkin.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,1+Nt/plotgap
		WRITE(11,*) enkin(j)
	END DO
	CLOSE(11)
	
	name_config = 'enpot.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,1+Nt/plotgap
		WRITE(11,*) enpot(j)
	END DO
	CLOSE(11)
	
	name_config = 'enstr.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO j=1,1+Nt/plotgap
		WRITE(11,*) enstr(j)
	END DO
	CLOSE(11)

	END SUBROUTINE saveresults