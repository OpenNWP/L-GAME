! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

program control

	! This controls the model run from the beginning to the end.

	use run_nml,        only: run_nml_setup, run_span_hr, dtime, &
	                          t_init, nlins, ncols, nlev
	use grid_generator, only: grid_setup, t_grid
	use io,             only: read_init, wp
	use manage_rkhevi,  only: rkhevi
	
	implicit none

	integer       :: run_span, time_step_counter
	real(wp)     :: t_0
	type(t_grid) :: grid

	! reading in all namelists so that we know what we have to do
	write(*,*) "Reading in run namelist ..."
	call run_nml_setup
	write(*,*) "... run namelist read."

	! allocating the grid
	write(*,*) "Allocating the grid ..."
	allocate(grid%z_geo_scal(nlins,ncols,nlev))
	allocate(grid%z_agl_scal(nlins,ncols,nlev))
	write(*,*) "... finished."

	! firstly, the grid generator needs to be called to calculate the grid properties
	write(*,*) "Setting up the grid ..."
	call grid_setup(grid)
	write(*,*) "... grid set up."

	! reading the initial state
	write(*,*) "Reading the initial state..."
	call read_init()
	write(*,*) "... initial state read."

	! the loop over the time steps
	t_0 = t_init
	run_span = 3600*run_span_hr
	time_step_counter = 0
	do while (t_0 < t_init + run_span + 300)
		
		! this is the RKHEVI routine performing the time stepping
		call rkhevi(grid)
		
        t_0 = t_0 + dtime
		time_step_counter = time_step_counter + 1
		write(*,*) "Step ", time_step_counter, " completed."
		
	enddo
  
end program control






