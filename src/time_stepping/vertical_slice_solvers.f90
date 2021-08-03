! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module vertical_slice_solvers

	! This module contains the implicit vertical routines (implicit part of the HEVI scheme).

	use run_nml,     only: nlins,ncols
	use definitions, only: t_grid

	implicit none
	
	private
	
	public :: three_band_solver_ver
	
	contains
	
	subroutine three_band_solver_ver(grid)

		type(t_grid), intent(in) :: grid  ! model grid

		! local variables
		integer                  :: ji,jk ! loop variables

		do ji=1,nlins
			do jk=1,ncols
		
				
		
			enddo
		enddo

	end subroutine three_band_solver_ver

end module vertical_slice_solvers
