!--------------------------------------------------------------------
!
!
! PURPOSE
!
! This subroutine saves a two dimensional real array in binary
! format
!
! AUTHORS
!
! B. Cloutier, B.K. Muite, P. Rigge
! 4 June 2012
!
! INPUT
!
! .. Scalars ..
!  Nx                         = number of modes in x - power of 2 for FFT
!  Ny                         = number of modes in y - power of 2 for FFT
!  plotnum                    = number of plot to be made
! .. Arrays ..
!  field                      = real data to be saved
!  name_config                = root of filename to save to
!
! .. Output   ..
! plotnum                     = number of plot to be saved
!
! LOCAL VARIABLES
!
! .. Scalars ..
!  i                          = loop counter in x direction
!  j                          = loop counter in y direction
!  iol                        = size of file
!  count                      = counter
!  ind                        = index to strip data/ from the front of filepath
! .. Arrays ..
!     number_file             = datbin file name
SUBROUTINE savedata(Nx,Ny,plotnum,name_config,field)
  IMPLICIT NONE
  INTEGER(KIND=4), INTENT(IN)                           :: Nx,Ny
  INTEGER(KIND=4), INTENT(IN)                           :: plotnum
  REAL(KIND=8), DIMENSION(1:NX,1:NY), INTENT(IN)        :: field
  CHARACTER*100, INTENT(IN)                             :: name_config
  INTEGER(kind=4)                                       :: i,j,iol,count,ind
  CHARACTER*100                                         :: number_file
  ! write datbin
  ! create character array with full filename
  ind = index(name_config,' ') - 1
  WRITE(number_file,'(i0)') 10000000+plotnum
  number_file = name_config(1:ind)//number_file
  ind = index(number_file,' ') - 1
  number_file = number_file(1:ind)//'.datbin'
  INQUIRE( iolength=iol )       field(1,1)
  OPEN(unit=11,FILE=number_file,form="unformatted",&
       access="direct",recl=iol)
  count=1
  DO j=1,Ny
     DO i=1,Nx
        WRITE(11,rec=count) field(i,j)
        count=count+1
     END DO
  END DO
  CLOSE(11)

END SUBROUTINE savedata
