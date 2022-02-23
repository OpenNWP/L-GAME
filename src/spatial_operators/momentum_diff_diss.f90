! This ! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module momentum_diff_diss

  ! This module handles momentum diffusion and dissipation.
  
  use definitions,          only: t_grid,t_diag,t_irrev,t_state
  use divergence_operators, only: divv_h
  use gradient_operators,   only: grad_hor
  use run_nml,              only: nlins,ncols
  
  implicit none
  
  private
  
  public :: mom_diff_h
  public :: mom_diff_v
  
  contains
  
  subroutine mom_diff_h(state,diag,irrev,grid)
  
    ! This subroutine handles horizontal momentum diffusion.
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_irrev), intent(inout) :: irrev
    type(t_grid),  intent(in)    :: grid
    
    ! calculating the divergence of the horizontal wind field
    call divv_h(state%wind_u,state%wind_v,diag%scalar_placeholder,grid)
    ! calculating the gradient of the divergence
    call grad_hor(diag%scalar_placeholder,irrev%mom_diff_tend_x,irrev%mom_diff_tend_y,irrev%mom_diff_tend_z,grid)
  
  end subroutine mom_diff_h

  subroutine mom_diff_v
  
    ! This subroutine handles vertical momentum diffusion.

  end subroutine mom_diff_v

end module momentum_diff_diss
