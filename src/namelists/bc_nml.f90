! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module bc_nml

  ! In this namelist, the boundary conditions are configured.
  
  use definitions, only: wp
  use run_nml,     only: nlins,ncols
  
  implicit none
  
  integer           :: n_swamp           ! thickness of the swamp layer
  logical           :: lperiodic         ! periodic boundary conditions switch
  integer           :: dtime_bc          ! timestep for the boundary conditions update
  logical           :: lyrigid           ! switch for a rigid wall in y-direction
  character(len=64) :: bc_root_filename  ! root filename of the boundary conditions
  
  namelist /bc/n_swamp,lperiodic,dtime_bc,lyrigid,bc_root_filename
  
  contains
  
  subroutine bc_nml_setup

    ! local variables
    integer :: fileunit
    
    ! default values
    n_swamp = 5
    lperiodic = .false.
    dtime_bc = 10800._wp
    lyrigid = .false.
    bc_root_filename = "bc"
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=bc, unit=fileunit)
    
    close(fileunit)
  
  end subroutine bc_nml_setup

end module bc_nml
