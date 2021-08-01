! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module manage_rkhevi

	use definitions, only: t_grid, t_state

	implicit none

	contains
	
	subroutine rkhevi(state_old, state_new, state_tendency, grid)
		
		type(t_state), intent(inout) :: state_old      ! the state at the old timestep
		type(t_state), intent(inout) :: state_new      ! the state at the new timestep
		type(t_state), intent(inout) :: state_tendency ! the state containing the tendency
		type(t_grid), intent(inout)  :: grid            ! the grid of the model
		
	end subroutine rkhevi

end module manage_rkhevi
