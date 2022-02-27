! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module linear_combine_two_states

  ! This module contains functionality to interpolate two states.

  use definitions, only: t_state,wp,t_grid
  use constituents_nml, only: no_of_condensed_constituents
  
  implicit none
  
  private
  
  public :: lin_combination
  public :: interpolation_t
  
  contains

  subroutine lin_combination(state_0,state_1,state_out,coeff_0,coeff_1,grid)
  
    ! This subroutine performs a linear combination of two states.
    
    ! input variables and output
    type(t_state), intent(in)    :: state_0,state_1 ! the states which to combine
    type(t_state), intent(inout) :: state_out       ! the state into which the result shall be written
    real(wp),      intent(in)    :: coeff_0,coeff_1 ! coefficients for the linear combination
    type(t_grid),  intent(in)    :: grid            ! grid properties (containing the background state)
    
    state_out%rho = coeff_0*state_0%rho + coeff_1*state_1%rho
    state_out%rhotheta = coeff_0*state_0%rhotheta + coeff_1*state_1%rhotheta
    state_out%theta_pert = state_out%rhotheta/state_out%rho(:,:,:,no_of_condensed_constituents+1) - grid%theta_bg
    state_out%exner_pert = coeff_0*state_0%exner_pert + coeff_1*state_1%exner_pert
    state_out%wind_u = coeff_0*state_0%wind_u + coeff_1*state_1%wind_u
    state_out%wind_v = coeff_0*state_0%wind_v + coeff_1*state_1%wind_v
    state_out%wind_w = coeff_0*state_0%wind_w + coeff_1*state_1%wind_w
    state_out%condensed_rho_t = coeff_0*state_0%condensed_rho_t + coeff_1*state_1%condensed_rho_t
    state_out%temperature_soil = coeff_0*state_0%temperature_soil + coeff_1*state_1%temperature_soil
  
  end subroutine lin_combination
  
  subroutine interpolation_t(state_0,state_1,state_out,time_old,time_new,time_write,grid)
  
    ! This subroutine performs an interpolation in time.
  
    type(t_state), intent(in)    :: state_0    ! first state
    type(t_state), intent(in)    :: state_1    ! second state
    type(t_state), intent(inout) :: state_out  ! result state
    real(wp),      intent(in)    :: time_old   ! valid time of the first state
    real(wp),      intent(in)    :: time_new   ! valid time of the second state
    real(wp),      intent(in)    :: time_write ! time at which to write out the result state
    type(t_grid),  intent(in)    :: grid       ! grid quantities
    
    ! using the lin_combination subroutine for this task
    call lin_combination(state_0,state_1,state_out,1._wp-(time_write-time_old)/(time_new-time_old), &
    (time_write-time_old)/(time_new-time_old),grid)
  
  end subroutine interpolation_t

end module linear_combine_two_states






