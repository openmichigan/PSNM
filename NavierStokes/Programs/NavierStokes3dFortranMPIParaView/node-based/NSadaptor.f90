! fortran module for interacting with the ParaView CoProcessor
! loosely based on: 
! ParaView-3.14.1-Source/CoProcessing/Adaptors/FortranAdaptors/PhastaAdaptor/phastaNSadaptor.f
! -- Changed for SESE

! subroutine determines if coprocessing needed during the
! current simulation step or not, and if so, calls the coprocessor.
! some of the subroutines are supplied by ParaView's FortranAdaptorAPI,
! others have to be supplied by the programmer.      
module NSadaptor_module 
  use iso_c_binding
  implicit none
  public
  interface NSadaptor
    ! Originally also had a version that accepted 3D arrays, but in all other
    ! respects was identical. Decided this could lead to too much confusion.
    module procedure NSadaptor1D 
  end interface NSadaptor
contains
  subroutine NSadaptor1D(nx, ny, nz, xst, xen, yst, yen, zst, zen, &
                                     step, time, omegax, omegay, omegaz)
    ! nx, ny, nz     -- grid dimensions or entire mesh
    !                   used for setting whole extent
    ! xst, xen, etc. -- extents of current subdomain pencil
    ! step           -- current simulation time step
    ! time           -- current simulation time
    ! omega*          -- scalar array for the current time step
    ! flag           -- receives status from API calls
    integer, intent(in) :: nx, ny, nz, xst, xen, yst, yen, zst, zen, step
    real(kind=8), intent(in) :: time
    real(kind=8), dimension(:), intent (in) :: omegax, omegay, omegaz
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
          ! defined in VTKPointBasedDataSet.cxx
          ! takes the size of the entire grid, and the extents of the
          ! sub grid.
          call createcpimagedata(nx, ny, nz, xst, xen, yst, yen, zst, zen)
       end if
          
       ! defined in VTKPointBasedDataSet.cxx
       ! remember to null-terminate strings for C/C++
       call addfield(omegax, "realtempx"//char(0));
       call addfield(omegay, "realtempy"//char(0));
       call addfield(omegaz, "realtempz"//char(0));       

       ! defined in FortranAdaptorAPI.h
       call coprocess()
    end if
    
    return
  end subroutine NSadaptor1D 
end module NSadaptor_module
