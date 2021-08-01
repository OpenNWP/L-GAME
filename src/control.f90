! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

program control

	! This controls the model run from the beginning to the end.

	use run_nml,        only: run_nml_setup, run_span_hr, dtime, &
	                          t_init, nlins, ncols, nlev
	use definitions,    only: t_grid, t_state, wp
	use grid_generator, only: grid_setup
	use io,             only: read_init
	use manage_rkhevi,  only: rkhevi
	
	implicit none

	integer       :: run_span, time_step_counter
	real(wp)      :: t_0
	type(t_grid)  :: grid
	type(t_state) :: state_old, state_new, state_tendency, state_write

	! reading in all namelists so that we know what we have to do
	write(*,*) "Reading in run namelist ..."
	call run_nml_setup
	write(*,*) "... run namelist read."

	! allocating the grid
	write(*,*) "Allocating memory ..."
	allocate(grid%z_geo_scal(nlins,ncols,nlev))
	allocate(grid%z_agl_scal(nlins,ncols,nlev))
	! state at the old time step
	allocate(state_old%rho(nlins,ncols,nlev))
	allocate(state_old%rhotheta(nlins,ncols,nlev))
	allocate(state_old%exner_bg(nlins,ncols,nlev))
	allocate(state_old%theta_bg(nlins,ncols,nlev))
	allocate(state_old%theta_pert(nlins,ncols,nlev))
	allocate(state_old%exner_pert(nlins,ncols,nlev))
	write(*,*) "here"
	allocate(state_old%wind_h%x(nlins,ncols+1,nlev))
	write(*,*) "here"
	allocate(state_old%wind_h%y(nlins+1,ncols,nlev))
	allocate(state_old%wind_v(nlins,ncols,nlev))
	! state at the new time step
	allocate(state_new%rho(nlins,ncols,nlev))
	allocate(state_new%rhotheta(nlins,ncols,nlev))
	allocate(state_new%exner_bg(nlins,ncols,nlev))
	allocate(state_new%theta_bg(nlins,ncols,nlev))
	allocate(state_new%theta_pert(nlins,ncols,nlev))
	allocate(state_new%exner_pert(nlins,ncols,nlev))
	allocate(state_new%wind_h%x(nlins,ncols+1,nlev))
	allocate(state_new%wind_h%y(nlins+1,ncols,nlev))
	allocate(state_new%wind_v(nlins,ncols,nlev))
	! state containing the tendency
	allocate(state_tendency%rho(nlins,ncols,nlev))
	allocate(state_tendency%rhotheta(nlins,ncols,nlev))
	allocate(state_tendency%exner_bg(nlins,ncols,nlev))
	allocate(state_tendency%theta_bg(nlins,ncols,nlev))
	allocate(state_tendency%theta_pert(nlins,ncols,nlev))
	allocate(state_tendency%exner_pert(nlins,ncols,nlev))
	allocate(state_tendency%wind_h%x(nlins,ncols+1,nlev))
	allocate(state_tendency%wind_h%y(nlins+1,ncols,nlev))
	allocate(state_tendency%wind_v(nlins,ncols,nlev))
	! state to be written out
	allocate(state_write%rho(nlins,ncols,nlev))
	allocate(state_write%rhotheta(nlins,ncols,nlev))
	allocate(state_write%exner_bg(nlins,ncols,nlev))
	allocate(state_write%theta_bg(nlins,ncols,nlev))
	allocate(state_write%theta_pert(nlins,ncols,nlev))
	allocate(state_write%exner_pert(nlins,ncols,nlev))
	allocate(state_write%wind_h%x(nlins,ncols+1,nlev))
	allocate(state_write%wind_h%y(nlins+1,ncols,nlev))
	allocate(state_write%wind_v(nlins,ncols,nlev))
	write(*,*) "... finished."

	! firstly, the grid generator needs to be called to calculate the grid properties
	write(*,*) "Setting up the grid ..."
	call grid_setup(grid)
	write(*,*) "... grid set up."

	! reading the initial state
	write(*,*) "Reading the initial state..."
	call read_init(state_old)
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






