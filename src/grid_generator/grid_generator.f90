! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file contains the calculation of the grid properties as well as some type definitions.

module grid_generator

	use definitions, only: wp,t_grid
	use run_nml,     only: nlins,ncols,nlevs,dy,dx

	implicit none
	
	private
	
	public :: grid_setup
	
	contains
	
	subroutine grid_setup(grid)
	
		type(t_grid), intent(inout) :: grid
		! local variables
		real(wp) :: semimajor      ! Earth radius
		real(wp) :: semiminor      ! Earth radius
		real(wp) :: re             ! Earth radius
		real(wp) :: lat_left_lower ! latitude coordinate of lower left corner
		real(wp) :: lon_left_lower ! longitude coordinate of lower left corner
		real(wp) :: dlat           ! mesh size in y direction as angle
		real(wp) :: dlon           ! mesh size in x direction as angle
		integer  :: ji, jk, jl ! loop indices
		
		semiminor=6356752.314_wp
		semimajor=6378137.0_wp
		re = (semimajor*semimajor*semiminor)**(1._wp/3._wp)
		
		! setting the latitude and longitude coordinates of the scalar grid points
		! setting the dy of the model grid
		dlat = dy/re
		dlon = dx/re
		lat_left_lower = -(nlins - 1)/2*dlat
		lon_left_lower = -(ncols - 1)/2*dlon
		do ji=1,nlins
			grid%lat_scalar(ji) = lat_left_lower + dlat*(ji - 1)
		enddo
		do ji=1,ncols
			grid%lon_scalar(ji) = lon_left_lower + dlon*(ji - 1)
		enddo
		
		! setting the dy of the model grid
		do ji=1,nlins+1
			do jk=1,ncols
				do jl=1,nlevs
					grid%dy(ji,jk,jl) = dy
				enddo
			enddo
		enddo
	
	end subroutine grid_setup
	
end module grid_generator






