! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module run_nml

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
	use definitions, only: wp
	
	implicit none

	integer  :: nlins              ! number of lines
	integer  :: ncols              ! number of columns
	integer  :: nlev               ! number of levels
	real(wp) :: dy                 ! mesh size in y direction
	real(wp) :: dtime              ! time steprun_span
	integer  :: run_span_hr        ! run span in hours
	real     :: t_init             ! epoch time stamp of the initialization
	character(len = 100) :: run_id ! the ID of the run
	
	namelist /run/run_id,nlins,ncols,nlev,dy,dtime,run_span_hr

	contains

	subroutine run_nml_setup
	
		! local vaiables
		integer :: rc, fileunit
		
		nlins       = 101
		ncols       = 101
		nlev        = 101
		dy          = 500._wp
		dtime       = 5._wp
		run_span_hr = 1
		t_init      = 0._wp
		run_id      = "ideal"
		
        ! Open and read Namelist file.
        open(action="read", file="namelist.nml", iostat=rc, newunit=fileunit)
        read(nml=run, iostat=rc, unit=fileunit)
        
        close(fileunit)
		
		write(*,*) "Time step: ", dtime, " s."
		
	end subroutine run_nml_setup
	
end module run_nml





