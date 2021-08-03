! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

program control

	! This controls the model run from the beginning to the end.

	use run_nml,        only: run_nml_setup,run_span_hr,dtime, &
	                          t_init,nlins,ncols,nlays
	use definitions,    only: t_grid,t_state,wp,t_diag,t_bg,t_tend
	use grid_generator, only: grid_setup
	use io,             only: read_init
	use manage_rkhevi,  only: rkhevi
	use linear_combine_two_states, only: lin_combination
	
	implicit none

	integer       :: time_step_counter
	real(wp)      :: t_0, run_span
	type(t_grid)  :: grid
	type(t_state) :: state_old, state_new, state_write
	type(t_diag)  :: diag
	type(t_bg)    :: bg
	type(t_tend)  :: tend

	! reading in all namelists so that we know what we have to do
	write(*,*) "Reading in run namelist ..."
	call run_nml_setup
	write(*,*) "... run namelist read."

	! allocating memory
	write(*,*) "Allocating memory ..."
	allocate(grid%lat_scalar(nlins))
	allocate(grid%lon_scalar(ncols))
	allocate(grid%z_geo_scal(nlins,ncols,nlays))
	allocate(grid%z_agl_scal(nlins,ncols,nlays))
	allocate(grid%z_geo_w(nlins,ncols,nlays+1))
	allocate(grid%dy(nlins+1,ncols,nlays))
	allocate(grid%dx(nlins,ncols+1,nlays))
	allocate(grid%dz(nlins,ncols,nlays+1))
	allocate(grid%slope_x(nlins,ncols+1,nlays))
	allocate(grid%slope_y(nlins+1,ncols,nlays))
	allocate(grid%volume(nlins,ncols,nlays))
	allocate(grid%area_x(nlins,ncols+1,nlays))
	allocate(grid%area_y(nlins+1,ncols,nlays))
	allocate(grid%area_z(nlins,ncols,nlays+1))
	allocate(grid%inner_product_weights(nlins,ncols,nlays,6))
	! state at the old time step
	allocate(state_old%rho(nlins,ncols,nlays))
	allocate(state_old%rhotheta(nlins,ncols,nlays))
	allocate(state_old%theta_pert(nlins,ncols,nlays))
	allocate(state_old%theta(nlins,ncols,nlays))
	allocate(state_old%exner_pert(nlins,ncols,nlays))
	allocate(state_old%wind_u(nlins,ncols+1,nlays))
	allocate(state_old%wind_v(nlins+1,ncols,nlays))
	allocate(state_old%wind_w(nlins,ncols,nlays))
	! state at the new time step
	allocate(state_new%rho(nlins,ncols,nlays))
	allocate(state_new%rhotheta(nlins,ncols,nlays))
	allocate(state_new%theta_pert(nlins,ncols,nlays))
	allocate(state_new%theta(nlins,ncols,nlays))
	allocate(state_new%exner_pert(nlins,ncols,nlays))
	allocate(state_new%wind_u(nlins,ncols+1,nlays))
	allocate(state_new%wind_v(nlins+1,ncols,nlays))
	allocate(state_new%wind_w(nlins,ncols,nlays))
	! state containing the tendency
	allocate(tend%rho(nlins,ncols,nlays))
	allocate(tend%rhotheta(nlins,ncols,nlays))
	allocate(tend%wind_u(nlins,ncols+1,nlays))
	allocate(tend%wind_v(nlins+1,ncols,nlays))
	allocate(tend%wind_w(nlins,ncols,nlays+1))
	! state to be written out
	allocate(state_write%rho(nlins,ncols,nlays))
	allocate(state_write%rhotheta(nlins,ncols,nlays))
	allocate(state_write%theta_pert(nlins,ncols,nlays))
	allocate(state_write%theta(nlins,ncols,nlays))
	allocate(state_write%exner_pert(nlins,ncols,nlays))
	allocate(state_write%wind_u(nlins,ncols+1,nlays))
	allocate(state_write%wind_v(nlins+1,ncols,nlays))
	allocate(state_write%wind_w(nlins,ncols,nlays))
	! background state
	allocate(bg%theta(nlins,ncols,nlays))
	allocate(bg%exner(nlins,ncols,nlays))
	! state containing diagnostic quantities
	allocate(diag%e_kin(nlins,ncols,nlays))
	allocate(diag%e_kin_grad_x(nlins,ncols+1,nlays))
	allocate(diag%e_kin_grad_y(nlins+1,ncols,nlays))
	allocate(diag%e_kin_grad_z(nlins,ncols,nlays+1))
	allocate(diag%scalar_placeholder(nlins,ncols,nlays))
	allocate(diag%u_placeholder(nlins,ncols+1,nlays))
	allocate(diag%v_placeholder(nlins+1,ncols,nlays))
	write(*,*) "... finished."

	! firstly, the grid generator needs to be called to calculate the grid properties
	write(*,*) "Setting up the grid ..."
	call grid_setup(grid)
	write(*,*) "... grid set up."

	call lin_combination(state_old,state_old,state_new,1._wp,0._wp,bg)
    
	! reading the initial state
	write(*,*) "Reading the initial state..."
	call read_init(state_old)
	write(*,*) "... initial state read."

	! the loop over the time steps
	t_0 = t_init
	run_span = 3600._wp*run_span_hr
	time_step_counter = 0
	do while (t_0 < t_init + run_span + 300)
		
    	call lin_combination(state_new,state_new,state_old,1._wp,0._wp,bg)
    	
		! this is the RKHEVI routine performing the time stepping
		call rkhevi(state_old,state_new,tend,bg,grid,diag,time_step_counter)
		
        t_0 = t_0 + dtime
		time_step_counter = time_step_counter + 1
		write(*,*) "Step ", time_step_counter, " completed."
		
	enddo
	
	! deallocating the memory
	write(*,*) "Deallocating memory ..."
	deallocate(grid%lat_scalar)
	deallocate(grid%lon_scalar)
	deallocate(grid%z_geo_scal)
	deallocate(grid%z_agl_scal)
	deallocate(grid%z_geo_w)
	deallocate(grid%dy)
	deallocate(grid%dx)
	deallocate(grid%dz)
	deallocate(grid%slope_x)
	deallocate(grid%slope_y)
	deallocate(grid%volume)
	deallocate(grid%area_x)
	deallocate(grid%area_y)
	deallocate(grid%area_z)
	deallocate(grid%inner_product_weights)
	! state at the old time step
	deallocate(state_old%rho)
	deallocate(state_old%rhotheta)
	deallocate(state_old%theta_pert)
	deallocate(state_old%theta)
	deallocate(state_old%exner_pert)
	deallocate(state_old%wind_u)
	deallocate(state_old%wind_v)
	deallocate(state_old%wind_w)
	! state at the new time step
	deallocate(state_new%rho)
	deallocate(state_new%rhotheta)
	deallocate(state_new%theta_pert)
	deallocate(state_new%theta)
	deallocate(state_new%exner_pert)
	deallocate(state_new%wind_u)
	deallocate(state_new%wind_v)
	deallocate(state_new%wind_w)
	! state containing the tendency
	deallocate(tend%rho)
	deallocate(tend%rhotheta)
	deallocate(tend%wind_u)
	deallocate(tend%wind_v)
	deallocate(tend%wind_w)
	! state to be written out
	deallocate(state_write%rho)
	deallocate(state_write%rhotheta)
	deallocate(state_write%theta_pert)
	deallocate(state_write%theta)
	deallocate(state_write%exner_pert)
	deallocate(state_write%wind_u)
	deallocate(state_write%wind_v)
	deallocate(state_write%wind_w)
	! background state
	deallocate(bg%theta)
	deallocate(bg%exner)
	! state containing diagnostic quantities
	deallocate(diag%e_kin)
	deallocate(diag%e_kin_grad_x)
	deallocate(diag%e_kin_grad_y)
	deallocate(diag%e_kin_grad_z)
	deallocate(diag%scalar_placeholder)
	deallocate(diag%u_placeholder)
	deallocate(diag%v_placeholder)
	write(*,*) "... finished."
  
end program control





