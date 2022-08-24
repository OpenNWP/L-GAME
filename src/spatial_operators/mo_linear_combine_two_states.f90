! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_linear_combine_two_states

  ! This module contains functionality to interpolate two states.

  use mo_definitions,   only: t_state,wp,t_grid
  use constituents_nml, only: n_condensed_constituents
  
  implicit none
  
  contains

  subroutine lin_combination(state_1,state_2,state_out,coeff_1,coeff_2,grid)
  
    ! This subroutine performs a linear combination of two states.
    
    ! input variables and output
    type(t_state), intent(in)    :: state_1,state_2 ! the states which to combine
    type(t_state), intent(inout) :: state_out       ! the state into which the result shall be written
    real(wp),      intent(in)    :: coeff_1,coeff_2 ! coefficients for the linear combination
    type(t_grid),  intent(in)    :: grid            ! grid properties (containing the background state)
    
    !$omp parallel workshare
    state_out%rho = coeff_1*state_1%rho + coeff_2*state_2%rho
    state_out%rhotheta_v = coeff_1*state_1%rhotheta_v + coeff_2*state_2%rhotheta_v
    state_out%theta_v_pert = state_out%rhotheta_v/state_out%rho(:,:,:,n_condensed_constituents+1) - grid%theta_v_bg
    state_out%exner_pert = coeff_1*state_1%exner_pert + coeff_2*state_2%exner_pert
    state_out%wind_u = coeff_1*state_1%wind_u + coeff_2*state_2%wind_u
    state_out%wind_v = coeff_1*state_1%wind_v + coeff_2*state_2%wind_v
    state_out%wind_w = coeff_1*state_1%wind_w + coeff_2*state_2%wind_w
    state_out%temperature_soil = coeff_1*state_1%temperature_soil + coeff_2*state_2%temperature_soil
    !$omp end parallel workshare
  
  end subroutine lin_combination

end module mo_linear_combine_two_states






