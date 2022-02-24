! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module boundaries

  ! This module handles everything dealing with boundary conditions.

  use definitions, only: t_bc

  implicit none
  
  private
  
  public :: update_boundaries
  public :: rescale_tend
  
  contains
  
  subroutine update_boundaries(bc)
  
    ! updates the boundary conditions
    
    ! input arguments and output
    type(t_bc) :: bc
    
  end subroutine update_boundaries
  
  subroutine rescale_tend
  
    ! rescales the tendencies to account for boundary conditions
  
  end subroutine rescale_tend

end module boundaries








