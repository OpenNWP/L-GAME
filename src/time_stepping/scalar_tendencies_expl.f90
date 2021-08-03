! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module scalar_tendencies_expl

	! This module manages the calculation of the explicit component of the scalar tendencies.

	use definitions,          only: wp,t_grid,t_state,t_diag,t_tend
	use multiplications,      only: scalar_times_vector_scalar_h
	use divergence_operators, only: divv_h

	implicit none
	
	private
	
	public :: expl_scalar_tend
	
	contains
	
	subroutine expl_scalar_tend(grid,state,tend,diag)
	
		type(t_grid),  intent(in)    :: grid
		type(t_state), intent(in)    :: state
		type(t_tend),  intent(inout) :: tend
		type(t_diag),  intent(inout) :: diag
	
		! explicit mass density integration
		! calculating the mass density flux
		call scalar_times_vector_scalar_h(state%rho,state%wind_u,state%wind_v,diag%u_placeholder,diag%v_placeholder)
		! calculating the divergence of the mass density flux
		call divv_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder,grid)
		
		tend%rho(:,:,:) = diag%scalar_placeholder(:,:,:)
		
		! explicit potential temperature density integration
		! calculating the potential temperature density flux
		call scalar_times_vector_scalar_h(state%theta,diag%u_placeholder,diag%v_placeholder,diag%u_placeholder,diag%v_placeholder)
		! calculating the divergence of the potential temperature density flux
		call divv_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder,grid)
		
		tend%rhotheta(:,:,:) = diag%scalar_placeholder(:,:,:)
	
	end subroutine

end module scalar_tendencies_expl
