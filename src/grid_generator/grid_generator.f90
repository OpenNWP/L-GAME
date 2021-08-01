! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file contains the calculation of the grid properties as well as some type definitions.

module grid_generator

	use definitions, only: wp,t_grid
	use run_nml,     only: nlins,ncols,nlevs,dy

	implicit none
	
	private
	
	public :: grid_setup
	
	contains
	
	subroutine grid_setup(grid)
	
		type(t_grid), intent(inout) :: grid
		! local variables
		real(wp) :: semimajor  ! Earth radius
		real(wp) :: semiminor  ! Earth radius
		real(wp) :: re         ! Earth radius
		integer  :: ji, jk, jl ! loop indices
		
		semiminor=6356752.314_wp
		semimajor=6378137.0_wp
		re = (semimajor*semimajor*semiminor)**(1._wp/3._wp)
		
		! setting the dy of the model grid
		do ji=1,ncols+1
			do jk=1,nlins
				do jl=1,nlevs
					grid%dy(ji,jk,jl) = dy
				enddo
			enddo
		enddo
	
	end subroutine grid_setup
	
end module grid_generator






