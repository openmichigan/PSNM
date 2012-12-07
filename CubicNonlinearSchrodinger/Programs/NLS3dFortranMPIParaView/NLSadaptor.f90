! Fortran module for interacting with the ParaView CoProcessor
! loosely based on: 
! ParaView-3.14.1-Source/CoProcessing/Adaptors/FortranAdaptors/PhastaAdaptor/phastaadaptor.f
! -- Changed for SESE

! Subroutine determines if coprocessing needed during the
! current simulation step or not, and if so, calls the coprocessor.
! Some of the subroutines are supplied by ParaView's FortranAdaptorAPI,
! others have to be supplied by the programmer.      
module NLSadaptor_module 
  use iso_c_binding
  implicit none
  public
  ! using an interface incase need arises to overload 1D, 2D versions
  interface NLSadaptor 
    module procedure NLSadaptor3D 
  end interface NLSadaptor 
contains
  subroutine NLSadaptor3D(nx, ny, nz, xst, xen, yst, yen, zst, zen, &
                                     step, time, a)
    ! nx, ny, nz -- grid dimensions
    !               used for setting whole extent
    ! xst, xen, etc -- extents of current subdomain
    ! step       -- current simulation time step
    ! time       -- current simulation time
    ! a          -- scalar array for the current time step
    integer, intent(in) :: nx, ny, nz, xst, xen, yst, yen, zst, zen, step
    real(kind=8), intent(in) :: time
    complex(kind=8), dimension(:,:,:), intent (in) :: a 
    integer :: flag

    flag = 0
    ! check if processing this time step
    ! defined in FortranAdaptorAPI.h
!   print *, "requestdatadescription"
!   print *,"step, time, flag: ", step, time, flag
    call requestdatadescription(step, time, flag)
        
    if (flag /= 0) then
       ! processing requested
       ! check if need to create grid
       ! defined in FortranAdaptorAPI.h
!      print *, "needtocreategrid"
       call needtocreategrid(flag)
       
       if (flag /= 0) then
          ! grid needed
          ! defined in adaptor.cxx
          ! takes the size of the entire grid, and the extents of the
          ! sub grid.
!         print *, "createcpimagedata"
          call createcpimagedata(nx, ny, nz, xst, xen, yst, yen, zst, zen)
       end if
          
       ! defined in adaptor.cxx
       ! call for each field of interest. Be sure to null-terminate the
       ! string for C++!
!      print *, "addfield"

       call addfield(a, "u_complex"//char(0))

       ! defined in FortranAdaptorAPI.h
!      print *,"calling coprocess()"
       call coprocess()
    end if
    
    return
  end subroutine NLSadaptor3D 
end module NLSadaptor_module
