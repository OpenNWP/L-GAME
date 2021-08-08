! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module explicit_vector_tendencies

	! this module manages the calculation of the explicit part of the wind tendencies

	use definitions,        only: t_grid,t_state,t_diag,t_tend,wp
	use inner_product,      only: kinetic_energy
	use gradient_operators, only: grad
	use run_nml,            only: nlins,ncols,nlays,llinear
	use vorticities,        only: calc_pot_vort
	use multiplications,    only: scalar_times_vector
	use thermodynamics,     only: gas_constant_diagnostics,spec_heat_cap_diagnostics_v, &
	                              spec_heat_cap_diagnostics_p
	use vorticity_flux,     only: calc_vorticity_flux_term
	use diff_nml,           only: lmom_diff_h
	use momentum_diff_diss, only: mom_diff_h

	implicit none
	
	private
	
	public :: vector_tendencies_expl
	
	contains

	subroutine vector_tendencies_expl(state,tend,diag,grid,slow_update_bool,rk_step,total_step_counter)
	
		type(t_state), intent(in)    :: state              ! state to use for calculating the tendencies
		type(t_tend),  intent(inout) :: tend               ! the tendency
		type(t_diag),  intent(inout) :: diag               ! diagnostic properties (e_kin is a diagnostic property)
		type(t_grid),  intent(in)    :: grid               ! grid propertiesgrid)
		logical,       intent(in)    :: slow_update_bool   ! switch to set wether the slow terms will be updated or not
		integer,       intent(in)    :: rk_step            ! Runge-Kutta step
		integer,       intent(in)    :: total_step_counter ! time step counter of the model integration
		
		! local variables
		real(wp)                     :: old_hor_pgrad_sound_weight ! old time step pressure gradient weight
		real(wp)                     :: new_hor_pgrad_sound_weight ! new time step pressure gradient weight
		
		old_hor_pgrad_sound_weight = 0.5_wp - spec_heat_cap_diagnostics_v(1)/spec_heat_cap_diagnostics_p(1)
		new_hor_pgrad_sound_weight = 1._wp - old_hor_pgrad_sound_weight
		
		! momentum advection
		if (((slow_update_bool .and. rk_step == 2) .or. total_step_counter == 0) .and. .not. llinear) then
			! calculating the mass flux density
			call scalar_times_vector(state%rho,state%wind_u,state%wind_v,state%wind_w, &
			diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
			! calculating the potential vorticity
			call calc_pot_vort(state,diag,grid)
			! calculating the potential voritcity flux term
			call calc_vorticity_flux_term(diag,grid)
			
			! Kinetic energy is prepared for the gradient term of the Lamb transformation.
			call kinetic_energy(state,diag,grid)
			! taking the gradient of the kinetic energy
			call grad(diag%e_kin,diag%e_kin_grad_x,diag%e_kin_grad_y,diag%e_kin_grad_z,grid)
		endif
		
		! momentum diffusion and dissipation (only updated at the first RK step and if advection is updated as well)
    	if (rk_step == 1 .and. slow_update_bool) then
			! horizontal momentum diffusion
			if (lmom_diff_h) then
				call mom_diff_h(state,diag,grid)
			endif
		endif
		
		! adding up the explicit wind tendencies
		! x-direction
		tend%wind_u(:,:,:) = diag%p_grad_acc_l_u(:,:,:) + diag%p_grad_acc_nl_u(:,:,:) &
		! momentum advection
		- diag%e_kin_grad_x(:,:,:) + diag%pot_vort_tend_x(:,:,:) & 
		! momentum diffusion
		+ diag%mom_diff_tend_x(:,:,:)
		! y-direction
		tend%wind_v(:,:,:) = diag%p_grad_acc_l_v(:,:,:) + diag%p_grad_acc_nl_v(:,:,:) &
		! momentum advection
		- diag%e_kin_grad_y(:,:,:) + diag%pot_vort_tend_y(:,:,:) &
		! momentum diffusion
		+ diag%mom_diff_tend_y(:,:,:)
		! z-direction
		tend%wind_w(:,:,:) = gas_constant_diagnostics(1)/spec_heat_cap_diagnostics_p(1)*diag%p_grad_acc_l_w(:,:,:) + &
		gas_constant_diagnostics(1)/spec_heat_cap_diagnostics_p(1)*diag%p_grad_acc_nl_w(:,:,:) &
		! momentum advection
		- diag%e_kin_grad_z(:,:,:) + diag%pot_vort_tend_z(:,:,:) &
		! momentum diffusion
		+ diag%mom_diff_tend_z(:,:,:)
	
	end subroutine vector_tendencies_expl

end module explicit_vector_tendencies





