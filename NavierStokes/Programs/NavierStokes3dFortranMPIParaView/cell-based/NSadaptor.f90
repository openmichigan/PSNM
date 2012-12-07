! Fortran module for interacting with the ParaView CoProcessor
! loosely based on: 
! ParaView-3.14.1-Source/CoProcessing/Adaptors/FortranAdaptors/PhastaAdaptor/phastaadaptor.f
! -- Changed for SESE

! Subroutine determines if coprocessing needed during the
! current simulation step or not, and if so, calls the coprocessor.
! Some of the subroutines are supplied by ParaView's FortranAdaptorAPI,
! others have to be supplied by the programmer.      
module NSadaptor_module 
  use iso_c_binding
  implicit none
  public
  interface NSadaptor 
    module procedure NSadaptor3D 
  end interface NSadaptor 
contains
  subroutine NSadaptor3D(nx, ny, nz, xst, xen, yst, yen, zst, zen, &
                                     step, time, omegax, omegay, omegaz)
    ! nx, ny, nz -- grid dimensions
    !               used for setting whole extent
    ! step       -- current simulation time step
    ! time       -- current simulation time
    ! omega*      -- scalar arrays for the current time step
    integer, intent(in) :: nx, ny, nz, xst, xen, yst, yen, zst, zen, step
    real(kind=8), intent(in) :: time
    real(kind=8), dimension(:,:,:), intent (in) :: omegax, omegay, omegaz
    integer :: flag

    ! check if processing this time step
    ! defined in FortranAdaptorAPI.h
    call requestdatadescription(step, time, flag)
        
    if (flag /= 0) then
       ! processing requested
       ! check if need to create grid
       ! defined in FortranAdaptorAPI.h
       call needtocreategrid(flag)
       
       if (flag /= 0) then
          ! grid needed
          ! defined in adaptor.cxx
          ! takes the size of the entire grid, and the extents of the
          ! sub grid.
          call createcpimagedata(nx, ny, nz, xst, xen, yst, yen, zst, zen)
       end if
          
       ! defined in adaptor.cxx
       ! call for each field of interest. Be sure to null-terminate the
       ! string for C++!
       call addfield(omegax, "realtempx"//char(0))
       call addfield(omegay, "realtempy"//char(0))
       call addfield(omegaz, "realtempz"//char(0))       


       ! defined in FortranAdaptorAPI.h
       call coprocess()
    end if
    
    return
  end subroutine NSadaptor3D 
end module NSadaptor_module
