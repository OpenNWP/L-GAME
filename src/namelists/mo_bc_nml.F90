! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_bc_nml
  
  ! In this namelist, the boundary conditions are configured.
  
  use mo_definitions, only: wp
  use mo_run_nml,     only: ny,nx
  
  implicit none
  
  integer           :: n_swamp          ! thickness of the swamp layer
  logical           :: lperiodic        ! periodic boundary conditions switch
  logical           :: lfreeslip        ! free slip boundary conditions (surface) switch
  integer           :: dtime_bc         ! time step for the boundary conditions update
  character(len=64) :: bc_root_filename ! root filename of the boundary conditions
  real(wp)          :: t_latest_bc      ! latest boundary conditions update time
  
  namelist /bc/n_swamp,lperiodic,lfreeslip,dtime_bc,bc_root_filename
  
  contains
  
  subroutine bc_nml_setup
    
    ! local variables
    integer :: fileunit ! file unit of the namelist file
    
    ! default values
    n_swamp = 5
    lperiodic = .false.
    lfreeslip = .false.
    dtime_bc = 10800
    bc_root_filename = "bc"
    t_latest_bc = 0._wp
    
    ! Open and read namelist file.
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=bc,unit=fileunit)
    
    close(fileunit)
    
  end subroutine bc_nml_setup

end module mo_bc_nml
