! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module inner_product

	use run_nml,     only: nlins,ncols
	use definitions, only: t_state,t_diag,t_grid,wp
	
	implicit none
	
	private
	
	public :: kinetic_energy
	
	contains

	subroutine kinetic_energy(state,diag,grid)
	
		! This subroutine calculates the specific kinetic energy.
		
		type(t_state), intent(in)    :: state ! state to use for calculating e_kin
		type(t_diag),  intent(inout) :: diag  ! diagnostic properties (e_kin is a diagnostic property)
		type(t_grid),  intent(in)    :: grid  ! grid properties
		! local variables
		integer :: ji,jk                     ! loop indices
		
		do ji=1,nlins
			do jk=1,ncols
				diag%e_kin(ji,jk,:) = ( &
				state%wind_u(ji  ,jk+1,:)**2*grid%area_x(ji  ,jk+1,:)*grid%dx(ji  ,jk+1,:) + &
				state%wind_v(ji  ,jk  ,:)**2*grid%area_y(ji  ,jk  ,:)*grid%dy(ji  ,jk  ,:) + &
				state%wind_u(ji  ,jk  ,:)**2*grid%area_x(ji  ,jk  ,:)*grid%dx(ji  ,jk  ,:) + &
				state%wind_v(ji+1,jk  ,:)**2*grid%area_y(ji+1,jk  ,:)*grid%dy(ji+1,jk  ,:))/ &
				(4._wp*grid%volume(ji,jk,:))
			enddo
		enddo
	
	end subroutine kinetic_energy

end module inner_product
