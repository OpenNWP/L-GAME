! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module manage_rkhevi

	! In this file, the RKHEVI time stepping is managed.

	use definitions,                only: t_grid,t_state,t_diag,t_bg,wp,t_tend
	use linear_combine_two_states,  only: lin_combination
	use run_nml,                    only: adv_sound_ratio,dtime
	use pressure_gradient,          only: manage_pressure_gradient
	use explicit_vector_tendencies, only: vector_tendencies_expl
	use thermodynamics,             only: spec_heat_cap_diagnostics_v, gas_constant_diagnostics

	implicit none
	
	private
	
	public :: rkhevi

	contains
	
	subroutine rkhevi(state_old,state_new,tend,bg,grid,diag,total_step_counter)
		
		type(t_state),  intent(inout) :: state_old          ! the state at the old timestep
		type(t_state),  intent(inout) :: state_new          ! the state at the new timestep
		type(t_tend),   intent(inout) :: tend               ! the tendency
		type(t_bg),     intent(in)    :: bg                 ! background state
		type(t_grid),   intent(in)    :: grid               ! the grid of the model
		type(t_diag),   intent(inout) :: diag               ! diagnostic quantities
		integer,        intent(in)    :: total_step_counter ! time step counter
		
		! local variables
		integer  :: rk_step          ! index of the Runge-Kutta step
		logical  :: slow_update_bool ! switch to determine if slow terms will be updated
		real(wp) :: delta_t_step     ! time step to use in this call
		
		! slow terms (momentum advection and diffusion) update switch
		slow_update_bool = .false.
		delta_t_step = dtime
		! check if slow terms have to be updated
		if (mod(total_step_counter, adv_sound_ratio) == 0) then
			! set the respective update switch to one
			slow_update_bool = .true.
			! delta_t is the large time step for the advection integration
			delta_t_step = adv_sound_ratio*dtime
		endif
		
		do rk_step = 1,2
		
			! state_old remains unchanged the whole time.
			! At rk_step == 1, it is state_old == state_new.
			
			! 1.) Explicit component of the momentum equation.
			! ------------------------------------------------
			! Update of the pressure gradient.
			if (rk_step == 0) then
				call manage_pressure_gradient()
			endif
			! Only the horizontal momentum is a forward tendency.
			call vector_tendencies_expl(state_new,tend,diag,grid,slow_update_bool,rk_step,total_step_counter)
	    	! time stepping for the horizontal momentum can be directly executed
			state_new%wind_u(:,:,:) = state_old%wind_u(:,:,:) + delta_t_step*tend%wind_u(:,:,:)
			state_new%wind_v(:,:,:) = state_old%wind_v(:,:,:) + delta_t_step*tend%wind_v(:,:,:)
			! Horizontal velocity can be considered to be updated from now on.

		enddo
		
		! in this case, a large time step has been taken, which we modify into a small step here    
		if (slow_update_bool .and. adv_sound_ratio > 1) then
			call lin_combination(state_old,state_new,state_new,1._wp-dtime/delta_t_step,dtime/delta_t_step,bg)
			! this is done for self-consistency
			state_new%rhotheta(:,:,:) = state_old%rhotheta(:,:,:)*(1 + spec_heat_cap_diagnostics_v(1)/gas_constant_diagnostics(1)* &
			((bg%exner(:,:,:) + state_new%exner_pert(:,:,:))/(bg%exner(:,:,:) + state_old%exner_pert(:,:,:)) - 1))
		endif
		
	end subroutine rkhevi

end module manage_rkhevi





