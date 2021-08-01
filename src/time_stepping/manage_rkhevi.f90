! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module manage_rkhevi

	use grid_generator, only: t_grid

	implicit none

	contains
	
	subroutine rkhevi(grid)
		
		type(t_grid), intent(inout) :: grid
		
	end subroutine rkhevi

end module manage_rkhevi
