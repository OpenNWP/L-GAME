! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module read_write_grid

  ! This module reads the grid from a file or writes the grid to a file. This is useful for efficiency.

  use definitions,       only: t_grid
  use set_initial_state, only: nc_check
  
  implicit none
  
  private
  
  public :: write_grid
  public :: read_grid
  
  contains
  
  subroutine write_grid(grid)
  
    type(t_grid), intent(in) :: grid
  
  end subroutine write_grid
  
  subroutine read_grid()
  
    ! This subroutine reads important grid properties from a file.
  
  end subroutine read_grid

end module read_write_grid








