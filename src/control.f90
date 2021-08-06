! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

program control

	! This controls the model run from the beginning to the end.

	use run_nml,                   only: run_nml_setup,run_span_hr,dtime, &
	                                     t_init,nlins,ncols,nlays,dt_write, &
	                                     lrestart,lideal,l3dvar,l4dvar
	use definitions,               only: t_grid,t_state,wp,t_diag,t_bg,t_tend
	use grid_generator,            only: grid_setup,bg_setup
	use io,                        only: restart,ideal,var_3d,var_4d,write_output
	use manage_rkhevi,             only: rkhevi
	use linear_combine_two_states, only: lin_combination,interpolation_t
	
	implicit none

	integer       :: time_step_counter
	real(wp)      :: t_0,run_span,t_write
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
	allocate(grid%lat_scalar(nlins+2))
	allocate(grid%lon_scalar(ncols+2))
	allocate(grid%z_geo_scal(nlins+2,ncols+2,nlays))
	allocate(grid%dx(nlins+2,ncols+1,nlays))
	allocate(grid%dy(nlins+1,ncols+2,nlays))
	allocate(grid%dz(nlins+2,ncols+2,nlays+1))
	allocate(grid%z_geo_u(nlins+2,ncols+1,nlays))
	allocate(grid%z_geo_v(nlins+1,ncols+2,nlays))
	allocate(grid%z_geo_w(nlins+2,ncols+2,nlays+1))
	allocate(grid%slope_x(nlins+2,ncols+1,nlays))
	allocate(grid%slope_y(nlins+1,ncols+2,nlays))
	allocate(grid%volume(nlins,ncols,nlays))
	allocate(grid%area_x(nlins,ncols+1,nlays))
	allocate(grid%area_y(nlins+1,ncols,nlays))
	allocate(grid%area_z(nlins,ncols,nlays+1))
	allocate(grid%area_dual_x(nlins+1,ncols,nlays+1))
	allocate(grid%area_dual_y(nlins,ncols+1,nlays+1))
	allocate(grid%area_dual_z(nlins+1,ncols+1,nlays))
	allocate(grid%z_geo_area_dual_z(nlins+1,ncols+1,nlays))
	allocate(grid%inner_product_weights(nlins,ncols,nlays,6))
	allocate(grid%fvec_x(nlins+1,ncols))
	allocate(grid%fvec_y(nlins,ncols+1))
	allocate(grid%fvec_z(nlins+1,ncols+1))
	allocate(grid%trsk_weights_u(nlins,ncols-1,6))
	allocate(grid%trsk_weights_v(nlins-1,ncols,4))
	! state at the old time step
	allocate(state_old%rho(nlins+2,ncols+2,nlays))
	allocate(state_old%rhotheta(nlins+2,ncols+2,nlays))
	allocate(state_old%theta_pert(nlins+2,ncols+2,nlays))
	allocate(state_old%exner_pert(nlins+2,ncols+2,nlays))
	allocate(state_old%wind_u(nlins+2,ncols+1,nlays))
	allocate(state_old%wind_v(nlins+1,ncols+2,nlays))
	allocate(state_old%wind_w(nlins+2,ncols+2,nlays+1))
	! state at the new time step
	allocate(state_new%rho(nlins+2,ncols+2,nlays))
	allocate(state_new%rhotheta(nlins+2,ncols+2,nlays))
	allocate(state_new%theta_pert(nlins+2,ncols+2,nlays))
	allocate(state_new%exner_pert(nlins+2,ncols+2,nlays))
	allocate(state_new%wind_u(nlins+2,ncols+1,nlays))
	allocate(state_new%wind_v(nlins+1,ncols+2,nlays))
	allocate(state_new%wind_w(nlins+2,ncols+2,nlays+1))
	! state containing the tendency
	allocate(tend%rho(nlins,ncols,nlays))
	allocate(tend%rhotheta(nlins,ncols,nlays))
	allocate(tend%wind_u(nlins,ncols-1,nlays))
	allocate(tend%wind_v(nlins-1,ncols,nlays))
	allocate(tend%wind_w(nlins,ncols,nlays+1))
	! state to be written out
	allocate(state_write%rho(nlins+2,ncols+2,nlays))
	allocate(state_write%rhotheta(nlins+2,ncols+2,nlays))
	allocate(state_write%theta_pert(nlins+2,ncols+2,nlays))
	allocate(state_write%exner_pert(nlins+2,ncols+2,nlays))
	allocate(state_write%wind_u(nlins+2,ncols+1,nlays))
	allocate(state_write%wind_v(nlins+1,ncols+2,nlays))
	allocate(state_write%wind_w(nlins+2,ncols+2,nlays+1))
	! background state
	allocate(bg%theta(nlins+2,ncols+2,nlays))
	allocate(bg%exner(nlins+2,ncols+2,nlays))
	! state containing diagnostic quantities
	allocate(diag%e_kin(nlins,ncols,nlays))
	allocate(diag%e_kin_grad_x(nlins,ncols-1,nlays))
	allocate(diag%e_kin_grad_y(nlins-1,ncols,nlays))
	allocate(diag%e_kin_grad_z(nlins,ncols,nlays+1))
	allocate(diag%p_grad_acc_l_u(nlins,ncols-1,nlays))
	allocate(diag%p_grad_acc_l_v(nlins-1,ncols,nlays))
	allocate(diag%p_grad_acc_l_w(nlins,ncols,nlays+1))
	allocate(diag%p_grad_acc_nl_u(nlins,ncols-1,nlays))
	allocate(diag%p_grad_acc_nl_v(nlins-1,ncols,nlays))
	allocate(diag%p_grad_acc_nl_w(nlins,ncols,nlays+1))
	allocate(diag%p_grad_acc_old_u(nlins,ncols-1,nlays))
	allocate(diag%p_grad_acc_old_v(nlins-1,ncols,nlays))
	allocate(diag%p_grad_acc_old_w(nlins,ncols,nlays+1))
	allocate(diag%pot_vort_tend_x(nlins,ncols-1,nlays))
	allocate(diag%pot_vort_tend_y(nlins-1,ncols,nlays))
	allocate(diag%pot_vort_tend_z(nlins,ncols,nlays+1))
	allocate(diag%mom_diff_tend_x(nlins,ncols-1,nlays))
	allocate(diag%mom_diff_tend_y(nlins-1,ncols,nlays))
	allocate(diag%mom_diff_tend_z(nlins,ncols,nlays+1))
	allocate(diag%scalar_placeholder(nlins+2,ncols+2,nlays))
	allocate(diag%u_placeholder(nlins+2,ncols+1,nlays))
	allocate(diag%v_placeholder(nlins+1,ncols+2,nlays))
	allocate(diag%u_10(nlins,ncols))
	allocate(diag%v_10(nlins,ncols))
	allocate(diag%mslp(nlins,ncols))
	allocate(diag%t_2(nlins,ncols))
	allocate(diag%z_eta_x(nlins+1,ncols,nlays+1))
	allocate(diag%z_eta_y(nlins,ncols+1,nlays+1))
	allocate(diag%z_eta_z(nlins+1,ncols+1,nlays))
	write(*,*) "... finished."

	! firstly, the grid generator needs to be called to calculate the grid properties
	write(*,*) "Setting up the grid ..."
	call grid_setup(grid)
	write(*,*) "... grid set up."

	! copying the old state to the new state
	call lin_combination(state_old,state_old,state_new,1._wp,0._wp,bg)
	
	! setting up the background state
	call bg_setup(grid,bg)
    
	! setting the initial state
	write(*,*) "Setting the initial state..."
	if (lrestart) then
		call restart(state_old,diag,bg,grid)
	elseif (lideal) then
		call ideal(state_old,diag,bg,grid)
	elseif (l3dvar) then
		call var_3d(state_old,diag,bg,grid)
	elseif (l4dvar) then
		call var_4d(state_old,diag,bg,grid)
	endif
	write(*,*) "... initial state set."
	
	call write_output(state_old,diag,0,grid,bg)

	! the loop over the time steps
	t_0 = t_init
	t_write = t_0 + dt_write
	run_span = 3600._wp*run_span_hr
	time_step_counter = 0
	do while (t_0 < t_init + run_span + 300)
		
    	call lin_combination(state_new,state_new,state_old,1._wp,0._wp,bg)
    	
		! this is the RKHEVI routine performing the time stepping
		call rkhevi(state_old,state_new,tend,bg,grid,diag,time_step_counter)
		
		! managing the calls to the output routine
		if (t_0 + dtime >= t_write) then
		
            call interpolation_t(state_old,state_new,state_write,t_0,t_0+dtime,t_write,bg)
			call write_output(state_write,diag,int((t_write-t_0)/60._wp),grid,bg)
		
			t_write = t_write + dt_write
		
		endif
		
        t_0 = t_0 + dtime
		time_step_counter = time_step_counter + 1
		write(*,*) "Step ", time_step_counter, " completed."
		
	enddo
	
	! deallocating the memory
	write(*,*) "Deallocating memory ..."
	deallocate(grid%lat_scalar)
	deallocate(grid%lon_scalar)
	deallocate(grid%z_geo_scal)
	deallocate(grid%dy)
	deallocate(grid%dx)
	deallocate(grid%dz)
	deallocate(grid%z_geo_u)
	deallocate(grid%z_geo_v)
	deallocate(grid%z_geo_w)
	deallocate(grid%slope_x)
	deallocate(grid%slope_y)
	deallocate(grid%volume)
	deallocate(grid%area_x)
	deallocate(grid%area_y)
	deallocate(grid%area_z)
	deallocate(grid%area_dual_x)
	deallocate(grid%area_dual_y)
	deallocate(grid%area_dual_z)
	deallocate(grid%z_geo_area_dual_z)
	deallocate(grid%inner_product_weights)
	deallocate(grid%fvec_x)
	deallocate(grid%fvec_y)
	deallocate(grid%fvec_z)
	deallocate(grid%trsk_weights_u)
	deallocate(grid%trsk_weights_v)
	! state at the old time step
	deallocate(state_old%rho)
	deallocate(state_old%rhotheta)
	deallocate(state_old%theta_pert)
	deallocate(state_old%exner_pert)
	deallocate(state_old%wind_u)
	deallocate(state_old%wind_v)
	deallocate(state_old%wind_w)
	! state at the new time step
	deallocate(state_new%rho)
	deallocate(state_new%rhotheta)
	deallocate(state_new%theta_pert)
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
	deallocate(state_write%exner_pert)
	deallocate(state_write%wind_u)
	deallocate(state_write%wind_v)
	deallocate(state_write%wind_w)
	! background state
	deallocate(bg%theta)
	deallocate(bg%exner)
	! state containing diagnostic quantities
	deallocate(diag%e_kin)
	deallocate(diag%p_grad_acc_l_u)
	deallocate(diag%p_grad_acc_l_v)
	deallocate(diag%p_grad_acc_l_w)
	deallocate(diag%p_grad_acc_nl_u)
	deallocate(diag%p_grad_acc_nl_v)
	deallocate(diag%p_grad_acc_nl_w)
	deallocate(diag%p_grad_acc_old_u)
	deallocate(diag%p_grad_acc_old_v)
	deallocate(diag%p_grad_acc_old_w)
	deallocate(diag%e_kin_grad_x)
	deallocate(diag%e_kin_grad_y)
	deallocate(diag%e_kin_grad_z)
	deallocate(diag%pot_vort_tend_x)
	deallocate(diag%pot_vort_tend_y)
	deallocate(diag%pot_vort_tend_z)
	deallocate(diag%mom_diff_tend_x)
	deallocate(diag%mom_diff_tend_y)
	deallocate(diag%mom_diff_tend_z)
	deallocate(diag%scalar_placeholder)
	deallocate(diag%u_placeholder)
	deallocate(diag%v_placeholder)
	deallocate(diag%u_10)
	deallocate(diag%v_10)
	deallocate(diag%mslp)
	deallocate(diag%t_2)
	deallocate(diag%z_eta_x)
	deallocate(diag%z_eta_y)
	deallocate(diag%z_eta_z)
	write(*,*) "... finished."
  
end program control





