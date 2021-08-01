! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module manage_rkhevi

	! In this file, the RKHEVI time stepping is managed.

	use definitions, only: t_grid, t_state

	implicit none

	contains
	
	subroutine rkhevi(state_old, state_new, state_tendency, grid)
		
		type(t_state), intent(inout) :: state_old      ! the state at the old timestep
		type(t_state), intent(inout) :: state_new      ! the state at the new timestep
		type(t_state), intent(inout) :: state_tendency ! the state containing the tendency
		type(t_grid), intent(inout)  :: grid           ! the grid of the model
		
		! local variables
		integer :: rk_step ! index of the Runge-Kutta step
		
		do rk_step = 1,2
		
			! state_old remains unchanged the whole time.
			! At rk_step == 1, it is state_old == state_new.		
		
		enddo
		
	end subroutine rkhevi

end module manage_rkhevi
