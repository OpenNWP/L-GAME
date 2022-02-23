! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module boundaries

  ! This module handles everything dealing with boundary conditions.

  use definitions, only: t_tend

  implicit none
  
  private
  
  public :: update_boundaries
  
  contains
  
  subroutine update_boundaries(tend_bc)
  
    ! updates the boundary conditions
    
    ! input arguments and output
    type(t_tend) :: tend_bc
    
  end subroutine update_boundaries

end module boundaries








