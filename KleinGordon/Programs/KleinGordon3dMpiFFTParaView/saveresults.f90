




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

	! Declare variables
	IMPLICIT NONE					 
	INTEGER(kind=4), INTENT(IN)							:: plotgap,Nt
	REAL(KIND=8), DIMENSION(1:1+Nt/plotgap), INTENT(IN)	:: enpot, enkin
	REAL(KIND=8), DIMENSION(1:1+Nt/plotgap), INTENT(IN)	:: en,enstr,time
	INTEGER(kind=4)										:: n
	CHARACTER*100										:: name_config

	name_config = 'tdata.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO n=1,1+Nt/plotgap
		WRITE(11,*) time(n)
	END DO
	CLOSE(11)

	name_config = 'en.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO n=1,1+Nt/plotgap
		WRITE(11,*) en(n)
	END DO
	CLOSE(11)

	name_config = 'enkin.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO n=1,1+Nt/plotgap
		WRITE(11,*) enkin(n)
	END DO
	CLOSE(11)
	
	name_config = 'enpot.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO n=1,1+Nt/plotgap
		WRITE(11,*) enpot(n)
	END DO
	CLOSE(11)
	
	name_config = 'enstr.dat' 
	OPEN(unit=11,FILE=name_config,status="UNKNOWN") 	
	REWIND(11)
	DO n=1,1+Nt/plotgap
		WRITE(11,*) enstr(n)
	END DO
	CLOSE(11)

	END SUBROUTINE saveresults