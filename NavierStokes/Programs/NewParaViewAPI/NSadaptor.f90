! 20130816 Mark Van Moer
! Fortran module for interacting with the ParaView/Catalyst coprocessor

module nsadaptor_glue
    implicit none
    public
contains
    subroutine nsadaptor(nx, ny, nz, step, time, scalar)
        ! nx, ny, nz -- grid dimensions
        ! step       -- current time step
        ! time       -- current time
        ! scalar     -- the data
        integer, intent(in) :: nx, ny, nz, step
        real(kind=8), intent(in) :: time
        real(kind=8), dimension(:,:), intent(in) :: scalar
        integer :: flag
        
        ! check if coprocessing this step
        ! defined in adaptor header in ParaView src
        call requestdatadescription(step, time, flag)
        if (flag .ne. 0) then
            ! coprocessing requested
            ! check if grid exists
            ! also from ParaView adaptor header
            call needtocreategrid(flag)

            if (flag .ne. 0) then
                ! grid needed
                ! defined by application developer
                call createcpimagedata(nx, ny, nz)
            end if

            ! also defined by application dev
            call addfield(scalar, "test")
    
            ! from ParaView adaptor header
            call coprocess()
        end if

        return
    end subroutine nsadaptor
end module nsadaptor_glue
