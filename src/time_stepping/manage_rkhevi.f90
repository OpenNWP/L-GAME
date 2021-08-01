! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module manage_rkhevi

	! In this file, the RKHEVI time stepping is managed.

	use definitions,               only: t_grid, t_state, wp
	use linear_combine_two_states, only: lin_combination
	use run_nml,                   only: adv_sound_ratio, dtime

	implicit none

	contains
	
	subroutine rkhevi(state_old, state_new, state_tendency, grid, total_step_counter)
		
		type(t_state),  intent(inout) :: state_old          ! the state at the old timestep
		type(t_state),  intent(inout) :: state_new          ! the state at the new timestep
		type(t_state),  intent(inout) :: state_tendency     ! the state containing the tendency
		type(t_grid),   intent(inout) :: grid               ! the grid of the model
		integer,        intent(in)    :: total_step_counter ! time step counter
		
		! local variables
		integer :: rk_step          ! index of the Runge-Kutta step
		logical :: slow_update_bool ! switch to determine if slow terms will be updated
		real(wp):: delta_t_step     ! time step to use in this call
		
		! slow terms (momentum advection and diffusion) update switch
		slow_update_bool = .false.
		delta_t_step = dtime
		! check if slow terms have to be updated
		if (mod(total_step_counter, adv_sound_ratio) == 0) then
			! set the respective update switch to one
			slow_update_bool = .true.
			! delta_t is the large time step for the advection integration
			delta_t_step = adv_sound_ratio*dtime;
		endif
		
		do rk_step = 1,2
		
			state_new%wind_u(:,:,:) = state_old%wind_u(:,:,:) + delta_t_step*state_tendency%wind_u(:,:,:)
			state_new%wind_v(:,:,:) = state_old%wind_v(:,:,:) + delta_t_step*state_tendency%wind_v(:,:,:)
			! Horizontal velocity can be considered to be updated from now on.
		
			! state_old remains unchanged the whole time.
			! At rk_step == 1, it is state_old == state_new.		
		
		enddo
		
		! in this case, a large time step has been taken, which we modify into a small step here    
		if (slow_update_bool .and. adv_sound_ratio > 1) then
			call lin_combination(state_old, state_new, state_new, 1 - dtime/delta_t_step, dtime/delta_t_step)
		endif
		
	end subroutine rkhevi

end module manage_rkhevi




