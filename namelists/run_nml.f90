! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module run_nml

    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
	use definitions, only: wp
	
	implicit none

	integer  :: nlins              ! number of lines
	integer  :: ncols              ! number of columns
	integer  :: nlevs              ! number of levels
	real(wp) :: dy                 ! mesh size in y direction at sea level
	real(wp) :: dx                 ! mesh size in x direction at sea level at the equator
	real(wp) :: dtime              ! time steprun_span
	integer  :: run_span_hr        ! run span in hours
	real     :: t_init             ! epoch time stamp of the initialization
	character(len = 100) :: run_id ! the ID of the run
	integer  :: adv_sound_ratio    ! ratio of advective to sound time step
	
	namelist /run/run_id,nlins,ncols,nlevs,dy,dx,dtime,run_span_hr,adv_sound_ratio

	contains

	subroutine run_nml_setup
	
		! local vaiables
		integer :: fileunit
		
		nlins           = 101
		ncols           = 101
		nlevs           = 80
		dy              = 800._wp
		dx              = 850._wp
		! this calculates the time step using the CFL criterion
		dtime           = 0.5_wp*dy/350._wp
		run_span_hr     = 63
		t_init          = 0._wp
		run_id          = "ideal"
		adv_sound_ratio = 5
		
        ! Open and read Namelist file.
        open(action="read", file="namelist.nml", newunit=fileunit)
        read(nml=run, unit=fileunit)
        
        close(fileunit)
        
        ! checking input data for correctness
        if (mod(nlins, 2) == 0) then
        	write(*,*) "Error: nlins must be odd. Aborting."
        	call exit(1)
        endif
        if (mod(ncols, 2) == 0) then
        	write(*,*) "Error: ncols must be odd. Aborting."
        	call exit(1)
        endif
		
		write(*,*) "Time step: ", dtime, " s."
		
	end subroutine run_nml_setup
	
end module run_nml









