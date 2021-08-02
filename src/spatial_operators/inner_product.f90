! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module inner_product

	use run_nml,     only: nlins,ncols,nlays
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
		integer :: ji,jk,jl                   ! loop indices
		
		do ji=1,nlins
			do jk=1,ncols
				do jl=1,nlays
					diag%e_kin(ji,jk,jl) = ( &
					state%wind_u(ji  ,jk+1,jl  )**2*grid%area_x(ji  ,jk+1,jl  )*grid%dx(ji  ,jk+1,jl  ) + &
					state%wind_v(ji  ,jk  ,jl  )**2*grid%area_y(ji  ,jk  ,jl  )*grid%dy(ji  ,jk  ,jl  ) + &
					state%wind_u(ji  ,jk  ,jl  )**2*grid%area_x(ji  ,jk  ,jl  )*grid%dx(ji  ,jk  ,jl  ) + &
					state%wind_v(ji+1,jk  ,jl  )**2*grid%area_y(ji+1,jk  ,jl  )*grid%dy(ji+1,jk  ,jl  ) + &
					state%wind_w(ji  ,jk  ,jl  )**2*grid%area_z(ji  ,jk  ,jl  )*grid%dz(ji  ,jk  ,jl  ) + &
					state%wind_w(ji  ,jk  ,jl+1)**2*grid%area_z(ji  ,jk  ,jl+1)*grid%dz(ji  ,jk  ,jl+1))/ &
					(4._wp*grid%volume(ji,jk,jl))
				enddo
			enddo
		enddo
	
	end subroutine kinetic_energy

end module inner_product
