! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module explicit_vector_tendencies

	use definitions,   only: t_grid,t_state,t_diag
	use inner_product, only: kinetic_energy

	implicit none
	
	private
	
	public :: vector_tendencies_expl
	
	contains

	subroutine vector_tendencies_expl(state,diag,grid,slow_update_bool,rk_step,total_step_counter)
	
		type(t_state), intent(in)    :: state              ! state to use for calculating e_kin
		type(t_diag),  intent(inout) :: diag               ! diagnostic properties (e_kin is a diagnostic property)
		type(t_grid),  intent(in)    :: grid               ! grid propertiesgrid)
		logical,       intent(in)    :: slow_update_bool   ! switch to set wether the slow terms will be updated or not
		integer,       intent(in)    :: rk_step            ! Runge-Kutta step
		integer,       intent(in)    :: total_step_counter ! time step counter of the model integration
		
		! momentum advection
		if ((slow_update_bool .and. rk_step == 2) .or. total_step_counter == 0) then
			! Kinetic energy is prepared for the gradient term of the Lamb transformation.
			call kinetic_energy(state,diag,grid)
		endif
	
	end subroutine vector_tendencies_expl

end module explicit_vector_tendencies
