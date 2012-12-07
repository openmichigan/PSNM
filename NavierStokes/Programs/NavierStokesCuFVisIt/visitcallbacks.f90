! VisIt: Functions that have to be defined even if not used
! VisIt: libsimV2 and libsimV2f are expecting these to be
! VisIt: visible, so easier not to put inside a module.
! VisIt: based on examples at wiki
  subroutine visitcommandcallback(cmd, lcmd, args, largs)
    implicit none
    include "visitfortransimV2interface.inc"
    character(len=8) :: cmd, args
    integer          :: lcmd, largs
    ! SIMSTATE common block
    integer :: runflag, simcycle
    real :: simtime
    common /SIMSTATE/ runflag, simcycle, simtime
    ! Handle commands defined in visitgetmetadata
    if (visitstrcmp(cmd, lcmd, "halt", 4) == 0) then
        runflag = 0
    elseif (visitstrcmp(cmd, lcmd, "step", 4) == 0) then
        !VisIt: this is the call that forces the simulation to
        !VisIt: be in an external function
        call simulate_one_timestep()
    elseif (visitstrcmp(cmd, lcmd, "run", 3) == 0) then
        runflag = 1
    endif
  end subroutine visitcommandcallback
  
  integer function visitbroadcastintfunction(value, sender)
    implicit none
    integer :: value, sender
    ! Replace with MPI communication if simulation is in parallel
    visitbroadcastintfunction = 0
  end function visitbroadcastintfunction

  integer function visitbroadcaststringfunction(str, lstr, sender)
    implicit none
    character(len=8) :: str
    integer :: lstr, sender
    ! Replace with MPI communication if simulation is in parallel
    visitbroadcaststringfunction = 0
  end function visitbroadcaststringfunction
 
  subroutine visitslaveprocesscallback()
    implicit none
    ! replace with MPI communication if simulation is in parallel
  end subroutine visitslaveprocesscallback

  integer function visitactivatetimestep()
    implicit none
    include "visitfortransimV2interface.inc"
    visitactivatetimestep = VISIT_OKAY
  end function visitactivatetimestep

  integer function visitgetmetadata()
    implicit none
    include "visitfortransimV2interface.inc"
    ! SIMSTATE common block
    integer :: runflag, simcycle
    real :: simtime
    common /SIMSTATE/ runflag, simcycle, simtime
    ! local variables
    integer :: metadata, mesh, vmd, cmd, err
    
    if (visitmdsimalloc(metadata) == VISIT_OKAY) then
        err = visitmdsimsetcycletime(metadata, simcycle, simtime)
        if (runflag == 1) then
            err = visitmdsimsetmode(metadata, VISIT_SIMMODE_RUNNING)
        else
            err = visitmdsimsetmode(metadata, VISIT_SIMMODE_STOPPED)
        endif
        
        ! set the mesh's properties
        if (visitmdmeshalloc(mesh) == VISIT_OKAY) then
           err = visitmdmeshsetname(mesh, "mesh", 4)
           err = visitmdmeshsetmeshtype(mesh, VISIT_MESHTYPE_RECTILINEAR)
           err = visitmdmeshsettopologicaldim(mesh, 2)
           err = visitmdmeshsetspatialdim(mesh, 2)
           err = visitmdmeshsetxlabel(mesh, "X", 1)
           err = visitmdmeshsetylabel(mesh, "Y", 1)
           err = visitmdsimaddmesh(metadata, mesh)
        endif

        ! add simulation commands
        if (visitmdcmdalloc(cmd) == VISIT_OKAY) then
           err = visitmdcmdsetname(cmd, "halt", 4)
           err = visitmdsimaddgenericcommand(metadata, cmd)
        endif
        if (visitmdcmdalloc(cmd) == VISIT_OKAY) then
           err = visitmdcmdsetname(cmd, "step", 4)
           err = visitmdsimaddgenericcommand(metadata, cmd)
        endif
        if (visitmdcmdalloc(cmd) == VISIT_OKAY) then
           err = visitmdcmdsetname(cmd, "run", 3)
           err = visitmdsimaddgenericcommand(metadata, cmd)
        endif

        ! add zonal scalar variable
        if (visitmdvaralloc(vmd) == VISIT_OKAY) then
            err = visitmdvarsetname(vmd, "nodal", 5)
            err = visitmdvarsetmeshname(vmd, "mesh", 4)
            err = visitmdvarsettype(vmd, VISIT_VARTYPE_SCALAR)  
            err = visitmdvarsetcentering(vmd, VISIT_VARCENTERING_NODE)
            err = visitmdsimaddvariable(metadata, vmd)
        endif
      endif
    visitgetmetadata = metadata
  end function visitgetmetadata

  integer function visitgetmesh(handle, domain, name, lname)
    implicit none
    character(len=8) :: name
    integer :: domain, lname, handle
    include "visitfortransimV2interface.inc"
    ! RECTMESH common block
    !VisIt: This common block is also used by
    !VisIt: simulate_one_timestep()
    integer, parameter :: nx = 1024, ny = 1024
    real, dimension(nx) :: rmx
    real, dimension(ny) :: rmy
    integer, dimension(3) :: rmdims
    integer :: rmndims
    common /RECTMESH/ rmdims, rmndims, rmx, rmy
    ! local variables
    integer :: h, x, y, z, modes, err
    
    h = VISIT_INVALID_HANDLE
    if (visitrectmeshalloc(h) == VISIT_OKAY) then
        err = visitvardataalloc(x)
        err = visitvardataalloc(y)
        err = visitvardatasetf(x, VISIT_OWNER_SIM, 1, nx, rmx)
        err = visitvardatasetf(y, VISIT_OWNER_SIM, 1, ny, rmy)
        err = visitrectmeshsetcoordsxy(h, x, y)
    endif
    visitgetmesh = h
  end function visitgetmesh

  integer function visitgetvariable(domain, name, lname)
    implicit none
    character(len=8) :: name
    integer :: domain, lname
    include "visitfortransimV2interface.inc"
    ! RECTMESH common block 
    integer :: numx = 1024, numy = 1024
    integer, dimension(3) :: rmdims
    real, dimension (:) :: rmx(1024), rmy(1024)
    real, dimension(:,:) :: nodal(1024, 1024)
    integer rmndims
    common /RECTMESH/ rmdims, rmndims, rmx, rmy, nodal
    !locals
    integer :: h, nvals, err
    rmdims = (/1024, 1024, 1/)
 
    h = VISIT_INVALID_HANDLE    
    if (visitstrcmp(name, lname, "nodal", 5) == 0) then
        if (visitvardataalloc(h) == VISIT_OKAY) then
            nvals = (rmdims(1) - 1) * (rmdims(2) - 1)
            err = visitvardatasetf(h, VISIT_OWNER_SIM, 1, nvals, nodal)
        endif
    endif

    visitgetvariable = h        
 
  end function visitgetvariable

  integer function visitgetcurve(name, lname)
    implicit none
    character(len=8) :: name
    integer :: lname
    include "visitfortransimV2interface.inc"
    visitgetcurve = VISIT_INVALID_HANDLE
  end function visitgetcurve

  integer function visitgetdomainlist(name, lname)
    implicit none
    character(len=8) :: name
    integer :: lname
    include "visitfortransimV2interface.inc"
    visitgetdomainlist = VISIT_INVALID_HANDLE
  end function visitgetdomainlist

  integer function visitgetdomainbounds(name, lname)
    implicit none
    character(len=8) :: name
    integer :: lname
    include "visitfortransimV2interface.inc"
    visitgetdomainbounds = VISIT_INVALID_HANDLE
  end function visitgetdomainbounds

  integer function visitgetdomainnesting(name, lname)
    implicit none
    character(len=8) :: name
    integer :: lname
    include "visitfortransimV2interface.inc"
    visitgetdomainnesting = VISIT_INVALID_HANDLE
  end function visitgetdomainnesting

  integer function visitgetmaterial(domain, name, lname)
    implicit none
    character(len=8) :: name
    integer :: domain, lname
    include "visitfortransimV2interface.inc"
    visitgetmaterial = VISIT_INVALID_HANDLE
  end function visitgetmaterial
