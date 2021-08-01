! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

program control

	! This controls the model run from the beginning to the end.

	use run_nml,        only: run_nml_setup, run_span_hr, dtime, &
	                          t_init, nlins, ncols, nlevs
	use definitions,    only: t_grid, t_state, wp
	use grid_generator, only: grid_setup
	use io,             only: read_init
	use manage_rkhevi,  only: rkhevi
	
	implicit none

	integer       :: time_step_counter
	real(wp)      :: t_0, run_span
	type(t_grid)  :: grid
	type(t_state) :: state_old, state_new, state_tendency, state_write

	! reading in all namelists so that we know what we have to do
	write(*,*) "Reading in run namelist ..."
	call run_nml_setup
	write(*,*) "... run namelist read."

	! allocating memory
	write(*,*) "Allocating memory ..."
	allocate(grid%z_geo_scal(nlins,ncols,nlevs))
	allocate(grid%z_agl_scal(nlins,ncols,nlevs))
	allocate(grid%dy(nlins+1,ncols,nlevs))
	allocate(grid%dx(nlins,ncols+1,nlevs))
	! state at the old time step
	allocate(state_old%rho(nlins,ncols,nlevs))
	allocate(state_old%rhotheta(nlins,ncols,nlevs))
	allocate(state_old%exner_bg(nlins,ncols,nlevs))
	allocate(state_old%theta_bg(nlins,ncols,nlevs))
	allocate(state_old%theta_pert(nlins,ncols,nlevs))
	allocate(state_old%exner_pert(nlins,ncols,nlevs))
	allocate(state_old%wind_u(nlins,ncols+1,nlevs))
	allocate(state_old%wind_v(nlins+1,ncols,nlevs))
	allocate(state_old%wind_w(nlins,ncols,nlevs))
	! state at the new time step
	allocate(state_new%rho(nlins,ncols,nlevs))
	allocate(state_new%rhotheta(nlins,ncols,nlevs))
	allocate(state_new%exner_bg(nlins,ncols,nlevs))
	allocate(state_new%theta_bg(nlins,ncols,nlevs))
	allocate(state_new%theta_pert(nlins,ncols,nlevs))
	allocate(state_new%exner_pert(nlins,ncols,nlevs))
	allocate(state_new%wind_u(nlins,ncols+1,nlevs))
	allocate(state_new%wind_v(nlins+1,ncols,nlevs))
	allocate(state_new%wind_w(nlins,ncols,nlevs))
	! state containing the tendency
	allocate(state_tendency%rho(nlins,ncols,nlevs))
	allocate(state_tendency%rhotheta(nlins,ncols,nlevs))
	allocate(state_tendency%exner_bg(nlins,ncols,nlevs))
	allocate(state_tendency%theta_bg(nlins,ncols,nlevs))
	allocate(state_tendency%theta_pert(nlins,ncols,nlevs))
	allocate(state_tendency%exner_pert(nlins,ncols,nlevs))
	allocate(state_tendency%wind_u(nlins,ncols+1,nlevs))
	allocate(state_tendency%wind_v(nlins+1,ncols,nlevs))
	allocate(state_tendency%wind_w(nlins,ncols,nlevs))
	! state to be written out
	allocate(state_write%rho(nlins,ncols,nlevs))
	allocate(state_write%rhotheta(nlins,ncols,nlevs))
	allocate(state_write%exner_bg(nlins,ncols,nlevs))
	allocate(state_write%theta_bg(nlins,ncols,nlevs))
	allocate(state_write%theta_pert(nlins,ncols,nlevs))
	allocate(state_write%exner_pert(nlins,ncols,nlevs))
	allocate(state_write%wind_u(nlins,ncols+1,nlevs))
	allocate(state_write%wind_v(nlins+1,ncols,nlevs))
	allocate(state_write%wind_w(nlins,ncols,nlevs))
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
	run_span = 3600._wp*run_span_hr
	time_step_counter = 0
	do while (t_0 < t_init + run_span + 300)
		
		! this is the RKHEVI routine performing the time stepping
		call rkhevi(state_old, state_new, state_tendency, grid)
		
        t_0 = t_0 + dtime
		time_step_counter = time_step_counter + 1
		write(*,*) "Step ", time_step_counter, " completed."
		
	enddo
	
	! deallocating the memory
	write(*,*) "Deallocating memory ..."
	deallocate(grid%z_geo_scal)
	deallocate(grid%z_agl_scal)
	! state at the old time step
	deallocate(state_old%rho)
	deallocate(state_old%rhotheta)
	deallocate(state_old%exner_bg)
	deallocate(state_old%theta_bg)
	deallocate(state_old%theta_pert)
	deallocate(state_old%exner_pert)
	deallocate(state_old%wind_u)
	deallocate(state_old%wind_v)
	deallocate(state_old%wind_w)
	! state at the new time step
	deallocate(state_new%rho)
	deallocate(state_new%rhotheta)
	deallocate(state_new%exner_bg)
	deallocate(state_new%theta_bg)
	deallocate(state_new%theta_pert)
	deallocate(state_new%exner_pert)
	deallocate(state_new%wind_u)
	deallocate(state_new%wind_v)
	deallocate(state_new%wind_w)
	! state containing the tendency
	deallocate(state_tendency%rho)
	deallocate(state_tendency%rhotheta)
	deallocate(state_tendency%exner_bg)
	deallocate(state_tendency%theta_bg)
	deallocate(state_tendency%theta_pert)
	deallocate(state_tendency%exner_pert)
	deallocate(state_tendency%wind_u)
	deallocate(state_tendency%wind_v)
	deallocate(state_tendency%wind_w)
	! state to be written out
	deallocate(state_write%rho)
	deallocate(state_write%rhotheta)
	deallocate(state_write%exner_bg)
	deallocate(state_write%theta_bg)
	deallocate(state_write%theta_pert)
	deallocate(state_write%exner_pert)
	deallocate(state_write%wind_u)
	deallocate(state_write%wind_v)
	deallocate(state_write%wind_w)
	write(*,*) "... finished."
  
end program control





