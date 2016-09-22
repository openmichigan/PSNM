	SUBROUTINE savedata(Nx,Ny,Nz,plotnum,name,field,u,v,decomp)
	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This subroutine saves a three dimensional real array in binary 
	! format
	!
	! INPUT
	!
	! .. Scalars ..
	!  Nx				= number of modes in x - power of 2 for FFT
	!  Ny				= number of modes in y - power of 2 for FFT
	!  Nz				= number of modes in z - power of 2 for FFT
	!  plotnum			= number of plot to be made
	! .. Arrays ..
	!  field 			= real data to be saved
	!  name_config		= root of filename to save to 
	!
	! .. Output	..	
	! plotnum			= number of plot to be saved
	! .. Special Structures ..
	!  decomp			= contains information on domain decomposition
	!					see http://www.2decomp.org/ for more information
	! LOCAL VARIABLES
	!
	! .. Scalars ..
	!  i				= loop counter in x direction
	!  j				= loop counter in y direction
	!  k				= loop counter in z direction
	!  count			= counter
	!  iol				= size of file
	! .. Arrays ..
	! 	number_file		= array to hold the number of the plot
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
	! 2DECOMP&FFT	 -- Domain decomposition and Fast Fourier Library
	!			(http://www.2decomp.org/index.html)
	! MPI library
	USE decomp_2d
	USE decomp_2d_fft
	USE decomp_2d_io
	IMPLICIT NONE					 
	INCLUDE 'mpif.h'
	! Declare variables
	INTEGER(KIND=4), INTENT(IN)						:: Nx,Ny,Nz
	INTEGER(KIND=4), INTENT(IN)						:: plotnum
	TYPE(DECOMP_INFO), INTENT(IN)					::  decomp
	REAL(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   							decomp%xst(2):decomp%xen(2),&
   							decomp%xst(3):decomp%xen(3)), &
   					                     INTENT(INOUT) :: field
	COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
   							decomp%xst(2):decomp%xen(2),&
   							decomp%xst(3):decomp%xen(3)), &
   					                     INTENT(IN) :: u,v
	CHARACTER*100, INTENT(IN)	     				:: name
	CHARACTER*100               					:: name_config
	INTEGER(kind=4)									:: i,j,k,iol,count,ind
	CHARACTER*100									:: number_file

	! create character array with full filename
		! write out using 2DECOMP&FFT MPI-IO routines
    ind=index(name,' ') -1
    name_config=name(1:ind)//'u'
    DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	 field(i,j,k)=REAL(u(i,j,k))
	END DO; END DO; END DO
	ind = index(name_config,' ') - 1
	WRITE(number_file,'(i0)') plotnum
	number_file = name_config(1:ind)//number_file
	ind = index(number_file,' ') - 1
	number_file = number_file(1:ind)//'.datbin'	
	CALL decomp_2d_write_one(1,field,number_file)

    ind=index(name,' ') -1
    name_config=name(1:ind)//'v'
    DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
	 field(i,j,k)=REAL(v(i,j,k))
	END DO; END DO; END DO
	ind = index(name_config,' ') - 1
	WRITE(number_file,'(i0)') plotnum
	number_file = name_config(1:ind)//number_file
	ind = index(number_file,' ') - 1
	number_file = number_file(1:ind)//'.datbin'	
	CALL decomp_2d_write_one(1,field,number_file)
	
	END SUBROUTINE savedata
