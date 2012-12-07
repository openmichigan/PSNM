! fortran module for interacting with the ParaView CoProcessor
! loosely based on: 
! ParaView-3.14.1-Source/CoProcessing/Adaptors/FortranAdaptors/PhastaAdaptor/phastaadaptor.f
! -- Changed for SESE

! subroutine determines if coprocessing needed during the
! current simulation step or not.
      
module ns2dcnadaptor
  implicit none
  public
contains
  subroutine navierstokescoprocessor(nx, ny, nz, step, time, omeg)
    ! nx, ny, nz -- grid dimensions
    ! step       -- current simulation time step
    ! time       -- current simulation time
    ! omega      -- scalar array for the current time step
    integer, intent(in) :: nx, ny, nz, step
    real(kind=8), intent(in) :: time
    real(kind=8), dimension(:,:), intent (in) :: omeg
    ! mvm: should this be "save" ?
    integer :: flag

    ! check if processing this time step
    ! defined in FortranAdaptorAPI.h
    call requestdatadescription(step, time, flag)
        
    if (flag .ne. 0) then
       ! processing requested
       ! check if need to create grid
       ! defined in FortranAdaptorAPI.h
       call needtocreategrid(flag)
       
       if (flag .ne. 0) then
          ! grid needed
          ! defined in adaptor.cxx
          call createcpimagedata(nx, ny, nz)
       end if
          
       ! defined in adaptor.cxx
       call addfield(omeg)
       
       ! defined in FortranAdaptorAPI.h
       call coprocess()
    end if
    
    return
  end subroutine navierstokescoprocessor
end module ns2dcnadaptor
