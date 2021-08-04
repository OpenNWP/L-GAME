! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module run_nml

	use definitions, only: wp
	
	implicit none

	integer           :: nlins              ! number of lines
	integer           :: ncols              ! number of columns
	integer		      :: nlays              ! number of levels
	integer           :: nlays_oro          ! number of levels following the orography
	real(wp)          :: dy                 ! mesh size in y direction at sea level
	real(wp)          :: dx                 ! mesh size in x direction at sea level at the equator
	real(wp)          :: dtime              ! time step
	real(wp)          :: toa                ! top of atmosphere
	real(wp)          :: sigma              ! vertical grid stretching parameter
	integer           :: run_span_hr        ! run span in hours
	real              :: t_init             ! epoch time stamp of the initialization
	integer           :: adv_sound_ratio    ! ratio of advective to sound time step
	real(wp)          :: semimajor          ! large halfaxis of the Earth
	real(wp)          :: semiminor          ! small halfaxis of the Earth
	real(wp)          :: re                 ! Earth radius
	integer           :: dt_write_min       ! output interval in minutes
	real(wp)          :: dt_write           ! output interval in seconds
	logical           :: lrestart           ! switch for restart runs
	logical           :: lideal             ! switch for analytic test cases
	logical           :: l3dvar             ! switch for 3d-Var
	logical           :: l4dvar             ! switch for 4d-Var
	character(len=64) :: scenario           ! scenario for ideal runs
	real(wp)          :: p0                 ! reference pressure
	real(wp)          :: omega              ! angular frequency of Earth rotation
	
	namelist /run/nlins,ncols,nlays,dy,dx,run_span_hr, &
	adv_sound_ratio,toa,dt_write_min,scenario

	contains

	subroutine run_nml_setup
	
		! local vaiables
		integer :: fileunit
		
		nlins           = 101
		ncols           = 101
		nlays           = 80
		dy              = 800._wp
		dx              = 850._wp
		run_span_hr     = 63
		t_init          = 0._wp
		adv_sound_ratio = 4
		toa             = 40000._wp
		sigma           = 1.3_wp
		semiminor       = 6356752.314_wp
		semimajor       = 6378137.0_wp
		re              = (semimajor*semimajor*semiminor)**(1._wp/3._wp)
		dt_write_min    = 60
		lrestart        = .false.
		lideal          = .true.
		l3dvar          = .false.
		l4dvar          = .false.
		scenario        = "standard"
		p0              = 100000._wp
		omega           = 7.292115e-5
		
        ! Open and read Namelist file.
        open(action="read", file="namelist.nml", newunit=fileunit)
        read(nml=run, unit=fileunit)
        
        close(fileunit)
        
		! this calculates the time step using the CFL criterion
		dtime           = 0.4_wp*dy/350._wp
		! calculating the output timestep in seconds
		dt_write        = 60._wp*dt_write_min
		! number of layers following orography
		nlays_oro       = int(0.66*nlays)
        
        ! checking input data for correctness
        if (mod(nlins, 2) == 0) then
        	write(*,*) "Error: nlins must be odd. Aborting."
        	call exit(1)
        endif
        if (mod(ncols, 2) == 0) then
        	write(*,*) "Error: ncols must be odd. Aborting."
        	call exit(1)
        endif
        if (toa <= 0) then
        	write(*,*) "Error: TOA must be positive."
        	call exit(1)
        endif
		
		write(*,*) "Time step: ", dtime, " s."
		
	end subroutine run_nml_setup
	
end module run_nml









