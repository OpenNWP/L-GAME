! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module bc_nml

  ! In this namelist, the boundary conditions are configured.
  
  use definitions, only: wp
  use run_nml,     only: nlins,ncols
  
  implicit none
  
  integer :: n_swamp   ! thickness of the swamp layer
  logical :: lperiodic ! periodic boundary conditions switch
  
  namelist /bc/n_swamp,lperiodic
  
  contains
  
  subroutine bc_nml_setup

    ! local variables
    integer :: fileunit
    
    ! default values
    n_swamp = 5
    lperiodic = .false.
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=bc, unit=fileunit)
    
    close(fileunit)
  
  end subroutine bc_nml_setup

end module bc_nml
