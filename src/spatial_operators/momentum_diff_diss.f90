! This ! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module momentum_diff_diss

  ! This module handles momentum diffusion and dissipation.
  
  use definitions,          only: t_grid,t_diag,t_state
  use divergence_operators, only: divv_h
  use gradient_operators,   only: grad_hor
  use run_nml,              only: nlins,ncols
  
  implicit none
  
  private
  
  public :: mom_diff_h
  public :: mom_diff_v
  
  contains
  
  subroutine mom_diff_h(state,diag,grid)
  
    ! This subroutine handles horizontal momentum diffusion.
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
    
    ! calculating the divergence of the horizontal wind field
    call divv_h(state%wind_u,state%wind_v,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:),grid)
    ! calculating the gradient of the divergence
    call grad_hor(diag%scalar_placeholder(2:nlins+1,2:ncols+1,:),diag%mom_diff_tend_x, &
    diag%mom_diff_tend_y,diag%mom_diff_tend_z,grid)
  
  end subroutine mom_diff_h

  subroutine mom_diff_v
  
    ! This subroutine handles vertical momentum diffusion.

  end subroutine mom_diff_v

end module momentum_diff_diss
