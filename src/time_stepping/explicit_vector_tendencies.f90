! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module explicit_vector_tendencies

	! this module manages the calculation of the explicit part of the wind tendencies

	use definitions,        only: t_grid,t_state,t_diag,t_tend
	use inner_product,      only: kinetic_energy
	use gradient_operators, only: grad
	use run_nml,            only: nlins,ncols,nlays
	use vorticities,        only: calc_pot_vort
	use multiplications,    only: scalar_times_vector_h

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
		
		! momentum advection
		if ((slow_update_bool .and. rk_step == 2) .or. total_step_counter == 0) then
			! calculating the mass flux density
			call scalar_times_vector_h(state%rho,state%wind_u,state%wind_v, &
			diag%u_placeholder,diag%v_placeholder)
			! calculating the potential vorticity
			call calc_pot_vort(state,diag,grid)
			! Kinetic energy is prepared for the gradient term of the Lamb transformation.
			call kinetic_energy(state,diag,grid)
			! taking the gradient of the kinetic energy
			call grad(diag%e_kin,diag%e_kin_grad_x,diag%e_kin_grad_y,diag%e_kin_grad_z,grid)
		endif
		
		! adding up the explicit wind tendencies
		tend%wind_u(:,:,:) = -diag%e_kin_grad_x(:,:,:) + diag%pot_vort_tend_x(:,:,:) + diag%mom_diff_tend_x(:,:,:)
		tend%wind_v(:,:,:) = -diag%e_kin_grad_y(:,:,:) + diag%pot_vort_tend_y(:,:,:) + diag%mom_diff_tend_y(:,:,:)
		tend%wind_w(:,:,:) = -diag%e_kin_grad_z(:,:,:) + diag%pot_vort_tend_z(:,:,:) + diag%mom_diff_tend_z(:,:,:)
	
	end subroutine vector_tendencies_expl

end module explicit_vector_tendencies





