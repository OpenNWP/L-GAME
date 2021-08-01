! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module run_nml

	use io, only: wp
	
	implicit none

	integer  :: nlins       ! number of lines
	integer  :: ncols       ! number of columns
	integer  :: nlev        ! number of levels
	real(wp) :: dy          ! mesh size in y direction
	real(wp) :: dtime       ! time steprun_span
	integer  :: run_span_hr ! run span in hours
	real     :: t_init      ! epoch time stamp of the initialization
	
	namelist /run/nlins,ncols,nlev,dy,dtime

	contains

	subroutine run_nml_setup
	
		nlins       = 101
		ncols       = 101
		nlev        = 101
		dy          = 500._wp
		dtime       = 5._wp
		run_span_hr = 1
		t_init      = 0._wp
		
		write(*,*) "Time step: ", dtime, " s."
		
	end subroutine run_nml_setup
	
end module run_nml
