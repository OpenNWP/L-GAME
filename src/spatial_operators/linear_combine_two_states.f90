! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module linear_combine_two_states

  use definitions, only: t_state,wp,t_grid
  
  implicit none
  
  private
  
  public :: lin_combination
  public :: interpolation_t
  
  contains

  subroutine lin_combination(state_0,state_1,state_out,coeff_0,coeff_1,grid)
  
    ! this performs a linear combination of two states
  
    type(t_state), intent(in)    :: state_0
    type(t_state), intent(in)    :: state_1
    type(t_state), intent(inout) :: state_out
    real(wp),      intent(in)    :: coeff_0
    real(wp),      intent(in)    :: coeff_1
    type(t_grid),  intent(in)    :: grid
    
    state_out%rho        = coeff_0*state_0%rho              + coeff_1*state_1%rho
    state_out%rhotheta   = coeff_0*state_0%rhotheta         + coeff_1*state_1%rhotheta
    state_out%theta_pert = state_out%rhotheta/state_out%rho - grid%theta_bg
    state_out%exner_pert = coeff_0*state_0%exner_pert       + coeff_1*state_1%exner_pert
    state_out%wind_u     = coeff_0*state_0%wind_u           + coeff_1*state_1%wind_u
    state_out%wind_v     = coeff_0*state_0%wind_v           + coeff_1*state_1%wind_v
    state_out%wind_w     = coeff_0*state_0%wind_w           + coeff_1*state_1%wind_w
  
  end subroutine lin_combination
  
  subroutine interpolation_t(state_0,state_1,state_out,time_old,time_new,time_write,grid)
  
    ! this performs an interpolation in time
  
    type(t_state), intent(in)    :: state_0
    type(t_state), intent(in)    :: state_1
    type(t_state), intent(inout) :: state_out
    real(wp),      intent(in)    :: time_old
    real(wp),      intent(in)    :: time_new
    real(wp),      intent(in)    :: time_write
    type(t_grid),  intent(in)    :: grid
    
    call lin_combination(state_0,state_1,state_out,1._wp-(time_write-time_old)/(time_new-time_old), &
    (time_write-time_old)/(time_new-time_old),grid)
  
  end subroutine interpolation_t

end module linear_combine_two_states






