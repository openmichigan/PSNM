!--------------------------------------------------------------------
!
!
! PURPOSE
!
! This subroutine saves the energy and times stored during the
! computation for the nonlinear sine-Gordon equation
!
! AUTHORS
!
! B. Cloutier, B.K. Muite, P. Rigge
! 4 June 2012
!
! INPUT
!
! .. Parameters ..
!  Nx                           = number of modes in x - power of 2 for FFT
!  Ny                           = number of modes in y - power of 2 for FFT
! .. Vectors ..
!  time                         = times at which save data
!  en                           = total energy
!  enstr                        = strain energy
!  enpot                        = potential energy
!  enkin                        = kinetic energy
!
! OUTPUT
!
!
! LOCAL VARIABLES
!
! .. Scalars ..
!  n                            = loop counter
! .. Arrays ..
!       name_config             = array to hold the filename
!
SUBROUTINE saveresults(Nt,plotgap,time,en,enstr,enkin,enpot)
  IMPLICIT NONE
  INTEGER(kind=4), INTENT(IN)                             :: plotgap,Nt
  REAL(KIND=8), DIMENSION(1+Nt/plotgap), INTENT(IN)       :: time
  REAL(KIND=8), DIMENSION(1+Nt/plotgap+1), INTENT(IN)     :: enpot,enkin
  REAL(KIND=8), DIMENSION(1+Nt/plotgap+1), INTENT(IN)     :: en,enstr
  INTEGER(kind=4)                                         :: j
  CHARACTER*100                                           :: name_config
  ! time
  name_config = 'tdata.dat'
  OPEN(unit=11,FILE=name_config,status="UNKNOWN")
  REWIND(11)
  DO j=1,1+Nt/plotgap
     WRITE(11,*) time(j)
  END DO
  CLOSE(11)
  ! energy
  name_config = 'en.dat'
  OPEN(unit=11,FILE=name_config,status="UNKNOWN")
  REWIND(11)
  DO j=1,1+Nt/plotgap+1
     WRITE(11,*) en(j)
  END DO
  CLOSE(11)
  ! kinetic energy
  name_config = 'enkin.dat'
  OPEN(unit=11,FILE=name_config,status="UNKNOWN")
  REWIND(11)
  DO j=1,1+Nt/plotgap+1
     WRITE(11,*) enkin(j)
  END DO
  CLOSE(11)
  ! potential energy
  name_config = 'enpot.dat'
  OPEN(unit=11,FILE=name_config,status="UNKNOWN")
  REWIND(11)
  DO j=1,1+Nt/plotgap+1
     WRITE(11,*) enpot(j)
  END DO
  CLOSE(11)
  ! strain energy
  name_config = 'enstr.dat'
  OPEN(unit=11,FILE=name_config,status="UNKNOWN")
  REWIND(11)
  DO j=1,1+Nt/plotgap+1
     WRITE(11,*) enstr(j)
  END DO
  CLOSE(11)
END SUBROUTINE saveresults
