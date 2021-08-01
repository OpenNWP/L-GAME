! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file contains the calculation of the grid properties as well as some type definitions.

module grid_generator

	use definitions, only: wp, t_grid

	implicit none
	
	private
	
	public :: grid_setup
	
	contains
	
	subroutine grid_setup(grid)
	
		type(t_grid), intent(inout) :: grid
		! local variables
		real(wp) :: re ! Earth radius
	
	end subroutine grid_setup
	
end module grid_generator






