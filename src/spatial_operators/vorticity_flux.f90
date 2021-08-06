! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module vorticity_flux

	! This module computes the voriticity flux term.
	
	use definitions, only: t_grid,t_diag,wp
	use run_nml,     only: nlins,ncols,nlays
	
	implicit none
	
	private
	
	public :: calc_vorticity_flux_term
	
	contains
	
	subroutine calc_vorticity_flux_term(diag,grid)

		type(t_diag), intent(inout) :: diag ! diagnostic quantities
		type(t_grid), intent(in)    :: grid ! model grid
		
		! local variables
		integer                     :: ji,jk,jl ! loop indices
		
		! horizontal velocity tendency due to vertical vorticity and horizontal wind (TRSK)
		! u
		do ji=1,nlins
			do jk=1,ncols-1
				diag%pot_vort_tend_x(ji,jk,:) = &
				grid%trsk_weights_u(ji,jk,1)*diag%v_placeholder(ji+1,jk+1,:)*0.25_wp* &
				(diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk,:)) &
				+ grid%trsk_weights_u(ji,jk,2)*diag%u_placeholder(ji+2,jk,:)*0.25_wp* &
				(diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji,jk,:)) &
				+ grid%trsk_weights_u(ji,jk,3)*diag%v_placeholder(ji,jk+1,:)*0.25_wp* &
				(diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji,jk,:)+diag%z_eta_z(ji,jk+1,:)) &
				+ grid%trsk_weights_u(ji,jk,4)*diag%v_placeholder(ji,jk+2,:)*0.25_wp* &
				(diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji,jk+2,:)) &
				+ grid%trsk_weights_u(ji,jk,5)*diag%u_placeholder(ji+1,jk+2,:)*0.25_wp* &
				(diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji,jk+2,:)+diag%z_eta_z(ji+1,jk+2,:)) &
				+ grid%trsk_weights_u(ji,jk,6)*diag%v_placeholder(ji+1,jk+2,:)*0.25_wp* &
				(diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk+2,:)+diag%z_eta_z(ji+1,jk+1,:))
			enddo
		enddo
		! v
		do ji=1,nlins-1
			do jk=1,ncols
				diag%pot_vort_tend_y(ji,jk,:) = &
				grid%trsk_weights_v(ji,jk,1)*diag%u_placeholder(ji+1,jk,:)*0.25_wp* &
				(diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji,jk,:)) &
				+ grid%trsk_weights_v(ji,jk,2)*diag%u_placeholder(ji+1,jk+1,:)*0.25_wp* &
				(diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji,jk+1,:)) &
				+ grid%trsk_weights_v(ji,jk,3)*diag%u_placeholder(ji+2,jk+1,:)*0.25_wp* &
				(diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+2,jk+1,:)) &
				+ grid%trsk_weights_v(ji,jk,4)*diag%u_placeholder(ji+2,jk,:)*0.25_wp* &
				(diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+2,jk,:))
			enddo
		enddo
		
		! horizontal velocity tendency due to horizontal vorticity and vertical wind
		
		! vertical velocity tendency due to horizontal vorticity and horizontal wind
	
	end subroutine calc_vorticity_flux_term

end module vorticity_flux








