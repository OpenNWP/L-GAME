! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module run_nml

	use io, only: wp
	
	implicit none

	integer  :: nlins ! number of lines
	integer  :: ncols ! number of columns
	integer  :: nlev  ! number of levels
	real(wp) :: dy    ! mesh size in y direction
	real(wp) :: dtime ! time step
	
	namelist /run/nlins,ncols,nlev,dy,dtime

	contains

	subroutine run_nml_setup
	
		nlins = 101
		ncols = 101
		nlev  = 101
		dy    = 500._wp
		dtime = 5._wp
		
	end subroutine run_nml_setup
	
end module run_nml
