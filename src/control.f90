! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

program control

  ! This controls the model run from the beginning to the end.

  use run_nml,                   only: run_nml_setup,run_span_hr,dtime, &
                                       t_init,nlins,ncols,nlays,lrestart, &
                                       lideal
  use io_nml,                    only: io_nml_setup,dt_write
  use constituents_nml,          only: constituents_nml_setup,no_of_condensed_constituents,no_of_constituents
  use diff_nml,                  only: diff_nml_setup
  use surface_nml,               only: surface_nml_setup,nsoillays
  use definitions,               only: t_grid,t_state,wp,t_diag,t_tend,t_irrev
  use grid_generator,            only: grid_setup,bg_setup
  use set_initial_state,         only: restart,ideal
  use write_out,                 only: write_output
  use manage_rkhevi,             only: rkhevi
  use linear_combine_two_states, only: lin_combination,interpolation_t
  use bc_nml,                    only: bc_nml_setup
  
  implicit none

  integer        :: time_step_counter
  real(wp)       :: t_0,run_span,t_write
  type(t_grid)   :: grid
  type(t_state)  :: state_old, state_new, state_write
  type(t_diag)   :: diag
  type(t_tend)   :: tend
  type(t_irrev)  :: irrev

  ! reading in all namelists so that we know what we have to do
  write(*,*) "Reading in run namelist ..."
  call run_nml_setup
  write(*,*) "... run namelist read."
  write(*,*) "Reading in diff namelist ..."
  call diff_nml_setup
  write(*,*) "... diff namelist read."
  write(*,*) "Reading in I/O namelist ..."
  call io_nml_setup
  write(*,*) "... I/O namelist read."
  write(*,*) "Reading in constituents namelist ..."
  call constituents_nml_setup()
  write(*,*) "... constituents namelist read."
  write(*,*) "Reading in surface namelist ..."
  call surface_nml_setup()
  write(*,*) "... surface namelist read."
  write(*,*) "Reading in boundary conditions namelist ..."
  call bc_nml_setup()
  write(*,*) "... boundary conditions namelist read."
  
  ! allocating memory
  write(*,*) "Allocating memory ..."
  allocate(grid%lat_scalar(nlins))
  allocate(grid%lon_scalar(ncols))
  allocate(grid%z_geo_scal(nlins,ncols,nlays))
  allocate(grid%dx(nlins,ncols+1,nlays))
  allocate(grid%dy(nlins+1,ncols,nlays))
  allocate(grid%dz(nlins,ncols,nlays+1))
  allocate(grid%z_geo_u(nlins,ncols+1,nlays))
  allocate(grid%z_geo_v(nlins+1,ncols,nlays))
  allocate(grid%z_geo_w(nlins,ncols,nlays+1))
  allocate(grid%slope_x(nlins,ncols+1,nlays))
  allocate(grid%slope_y(nlins+1,ncols,nlays))
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
  allocate(grid%theta_bg(nlins,ncols,nlays))
  allocate(grid%exner_bg(nlins,ncols,nlays))
  allocate(grid%exner_bg_grad_u(nlins,ncols-1,nlays))
  allocate(grid%exner_bg_grad_v(nlins-1,ncols,nlays))
  allocate(grid%exner_bg_grad_w(nlins,ncols,nlays+1))
  allocate(grid%sfc_albedo(nlins,ncols))
  allocate(grid%sfc_rho_c(nlins,ncols))
  allocate(grid%t_conduc_soil(nlins,ncols))
  allocate(grid%roughness_length(nlins,ncols))
  allocate(grid%is_land(nlins,ncols))
  allocate(grid%z_soil_interface(nsoillays+1))
  allocate(grid%z_soil_center(nsoillays))
  allocate(grid%t_const_soil(nlins,ncols))
  ! state at the old time step
  allocate(state_old%rho(nlins,ncols,nlays,no_of_constituents))
  allocate(state_old%rhotheta(nlins,ncols,nlays))
  allocate(state_old%theta_pert(nlins,ncols,nlays))
  allocate(state_old%exner_pert(nlins,ncols,nlays))
  allocate(state_old%condensed_rho_t(nlins,ncols,nlays,no_of_condensed_constituents))
  allocate(state_old%wind_u(nlins,ncols+1,nlays))
  allocate(state_old%wind_v(nlins+1,ncols,nlays))
  allocate(state_old%wind_w(nlins,ncols,nlays+1))
  allocate(state_old%temperature_soil(nlins,ncols,nsoillays))
  ! state at the new time step
  allocate(state_new%rho(nlins,ncols,nlays,no_of_constituents))
  allocate(state_new%rhotheta(nlins,ncols,nlays))
  allocate(state_new%theta_pert(nlins,ncols,nlays))
  allocate(state_new%exner_pert(nlins,ncols,nlays))
  allocate(state_new%condensed_rho_t(nlins,ncols,nlays,no_of_condensed_constituents))
  allocate(state_new%wind_u(nlins,ncols+1,nlays))
  allocate(state_new%wind_v(nlins+1,ncols,nlays))
  allocate(state_new%wind_w(nlins,ncols,nlays+1))
  allocate(state_new%temperature_soil(nlins,ncols,nsoillays))
  ! state containing the tendency
  allocate(tend%rho(nlins,ncols,nlays,no_of_constituents))
  allocate(tend%rhotheta(nlins,ncols,nlays))
  allocate(tend%wind_u(nlins,ncols+1,nlays))
  allocate(tend%wind_v(nlins+1,ncols,nlays))
  allocate(tend%wind_w(nlins,ncols,nlays+1))
  allocate(tend%condensed_rho_t(nlins,ncols,nlays,no_of_condensed_constituents))
  ! state to be written out
  allocate(state_write%rho(nlins,ncols,nlays,no_of_constituents))
  allocate(state_write%rhotheta(nlins,ncols,nlays))
  allocate(state_write%theta_pert(nlins,ncols,nlays))
  allocate(state_write%exner_pert(nlins,ncols,nlays))
  allocate(state_write%condensed_rho_t(nlins,ncols,nlays,no_of_condensed_constituents))
  allocate(state_write%wind_u(nlins,ncols+1,nlays))
  allocate(state_write%wind_v(nlins+1,ncols,nlays))
  allocate(state_write%wind_w(nlins,ncols,nlays+1))
  allocate(state_write%temperature_soil(nlins,ncols,nsoillays))
  ! type containing diagnostic quantities
  allocate(diag%e_kin(nlins,ncols,nlays))
  allocate(diag%e_kin_grad_x(nlins,ncols-1,nlays))
  allocate(diag%e_kin_grad_y(nlins-1,ncols,nlays))
  allocate(diag%e_kin_grad_z(nlins,ncols,nlays+1))
  allocate(diag%p_grad_acc_neg_l_u(nlins,ncols-1,nlays))
  allocate(diag%p_grad_acc_neg_l_v(nlins-1,ncols,nlays))
  allocate(diag%p_grad_acc_neg_l_w(nlins,ncols,nlays+1))
  allocate(diag%p_grad_acc_neg_nl_u(nlins,ncols-1,nlays))
  allocate(diag%p_grad_acc_neg_nl_v(nlins-1,ncols,nlays))
  allocate(diag%p_grad_acc_neg_nl_w(nlins,ncols,nlays+1))
  allocate(diag%p_grad_acc_old_u(nlins,ncols-1,nlays))
  allocate(diag%p_grad_acc_old_v(nlins-1,ncols,nlays))
  allocate(diag%p_grad_acc_old_w(nlins,ncols,nlays+1))
  allocate(diag%pot_vort_tend_x(nlins,ncols-1,nlays))
  allocate(diag%pot_vort_tend_y(nlins-1,ncols,nlays))
  allocate(diag%pot_vort_tend_z(nlins,ncols,nlays+1))
  allocate(diag%scalar_placeholder(nlins,ncols,nlays))
  allocate(diag%temperature_gas(nlins,ncols,nlays))
  allocate(diag%u_placeholder(nlins,ncols+1,nlays))
  allocate(diag%v_placeholder(nlins+1,ncols,nlays))
  allocate(diag%w_placeholder(nlins,ncols,nlays+1))
  allocate(diag%u_10(nlins,ncols))
  allocate(diag%v_10(nlins,ncols))
  allocate(diag%mslp(nlins,ncols))
  allocate(diag%t_2(nlins,ncols))
  allocate(diag%z_eta_x(nlins+1,ncols,nlays+1))
  allocate(diag%z_eta_y(nlins,ncols+1,nlays+1))
  allocate(diag%z_eta_z(nlins+1,ncols+1,nlays))
  allocate(diag%scalar_flux_resistance(nlins,ncols))
  allocate(diag%monin_obukhov_length(nlins,ncols))
  allocate(diag%power_flux_density_sensible(nlins,ncols))
  allocate(diag%power_flux_density_latent(nlins,ncols))
  allocate(diag%sfc_sw_in(nlins,ncols))
  allocate(diag%sfc_lw_out(nlins,ncols))
  allocate(diag%roughness_velocity(nlins,ncols))
  ! type containing irreversible quantities
  allocate(irrev%tke(nlins,ncols,nlays))
  allocate(irrev%mom_diff_tend_x(nlins,ncols-1,nlays))
  allocate(irrev%mom_diff_tend_y(nlins-1,ncols,nlays))
  allocate(irrev%mom_diff_tend_z(nlins,ncols,nlays+1))
  allocate(irrev%heating_diss(nlins,ncols,nlays))
  allocate(irrev%mass_source_rates(nlins,ncols,nlays,no_of_condensed_constituents+1))
  allocate(irrev%heat_source_rates(nlins,ncols,nlays,no_of_condensed_constituents))
  write(*,*) "... finished."

  ! firstly, the grid generator needs to be called to calculate the grid properties
  write(*,*) "Setting up the grid ..."
  call grid_setup(grid)
  write(*,*) "... grid set up."

  ! setting up the background state
  call bg_setup(grid)
  
  ! setting the initial state
  write(*,*) "Setting the initial state..."
  if (lrestart) then
    call restart(state_old,diag,grid)
  elseif (lideal) then
    call ideal(state_old,diag,grid)
  endif
  write(*,*) "... initial state set."
  
  ! writing out the initial state
  call write_output(state_old,diag,0,grid)
  
  ! initializing the wind tendencies with zero
  diag%e_kin_grad_x = 0._wp
  diag%e_kin_grad_y = 0._wp
  diag%e_kin_grad_z = 0._wp
  diag%pot_vort_tend_x = 0._wp
  diag%pot_vort_tend_y = 0._wp
  diag%pot_vort_tend_z = 0._wp
  irrev%mom_diff_tend_x = 0._wp
  irrev%mom_diff_tend_y = 0._wp
  irrev%mom_diff_tend_z = 0._wp
  tend%wind_u = 0._wp
  tend%wind_v = 0._wp
  tend%wind_w = 0._wp
  tend%condensed_rho_t = 0._wp
  
  ! copying the new state to the old state
  state_new = state_old
  
  ! the loop over the time steps
  t_0 = t_init
  t_write = t_0 + dt_write
  run_span = 3600._wp*run_span_hr
  time_step_counter = 0
  
  do while (t_0 < t_init + run_span + 300._wp)
    call lin_combination(state_new,state_new,state_old,0._wp,1._wp,grid)
      
    ! this is the RKHEVI routine performing the time stepping
    call rkhevi(state_old,state_new,tend,grid,diag,irrev,time_step_counter)
    
    ! managing the calls to the output routine
    if (t_0 + dtime >= t_write) then
    
      call interpolation_t(state_old,state_new,state_write,t_0,t_0+dtime,t_write,grid)
      call write_output(state_write,diag,int((t_write-t_init)/60._wp),grid)
    
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
  deallocate(grid%theta_bg)
  deallocate(grid%exner_bg)
  deallocate(grid%exner_bg_grad_u)
  deallocate(grid%exner_bg_grad_v)
  deallocate(grid%exner_bg_grad_w)
  deallocate(grid%sfc_albedo)
  deallocate(grid%sfc_rho_c)
  deallocate(grid%t_conduc_soil)
  deallocate(grid%roughness_length)
  deallocate(grid%is_land)
  deallocate(grid%z_soil_interface)
  deallocate(grid%z_soil_center)
  deallocate(grid%t_const_soil)
  ! state at the old time step
  deallocate(state_old%rho)
  deallocate(state_old%rhotheta)
  deallocate(state_old%theta_pert)
  deallocate(state_old%exner_pert)
  deallocate(state_old%condensed_rho_t)
  deallocate(state_old%wind_u)
  deallocate(state_old%wind_v)
  deallocate(state_old%wind_w)
  deallocate(state_old%temperature_soil)
  ! state at the new time step
  deallocate(state_new%rho)
  deallocate(state_new%rhotheta)
  deallocate(state_new%theta_pert)
  deallocate(state_new%exner_pert)
  deallocate(state_new%condensed_rho_t)
  deallocate(state_new%wind_u)
  deallocate(state_new%wind_v)
  deallocate(state_new%wind_w)
  deallocate(state_new%temperature_soil)
  ! state containing the tendency
  deallocate(tend%rho)
  deallocate(tend%rhotheta)
  deallocate(tend%wind_u)
  deallocate(tend%wind_v)
  deallocate(tend%wind_w)
  deallocate(tend%condensed_rho_t)
  ! state to be written out
  deallocate(state_write%rho)
  deallocate(state_write%rhotheta)
  deallocate(state_write%theta_pert)
  deallocate(state_write%exner_pert)
  deallocate(state_write%condensed_rho_t)
  deallocate(state_write%wind_u)
  deallocate(state_write%wind_v)
  deallocate(state_write%wind_w)
  deallocate(state_write%temperature_soil)
  ! type containing diagnostic quantities
  deallocate(diag%e_kin)
  deallocate(diag%p_grad_acc_neg_l_u)
  deallocate(diag%p_grad_acc_neg_l_v)
  deallocate(diag%p_grad_acc_neg_l_w)
  deallocate(diag%p_grad_acc_neg_nl_u)
  deallocate(diag%p_grad_acc_neg_nl_v)
  deallocate(diag%p_grad_acc_neg_nl_w)
  deallocate(diag%p_grad_acc_old_u)
  deallocate(diag%p_grad_acc_old_v)
  deallocate(diag%p_grad_acc_old_w)
  deallocate(diag%e_kin_grad_x)
  deallocate(diag%e_kin_grad_y)
  deallocate(diag%e_kin_grad_z)
  deallocate(diag%pot_vort_tend_x)
  deallocate(diag%pot_vort_tend_y)
  deallocate(diag%pot_vort_tend_z)
  deallocate(diag%scalar_placeholder)
  deallocate(diag%temperature_gas)
  deallocate(diag%u_placeholder)
  deallocate(diag%v_placeholder)
  deallocate(diag%w_placeholder)
  deallocate(diag%u_10)
  deallocate(diag%v_10)
  deallocate(diag%mslp)
  deallocate(diag%t_2)
  deallocate(diag%z_eta_x)
  deallocate(diag%z_eta_y)
  deallocate(diag%z_eta_z)
  deallocate(diag%scalar_flux_resistance)
  deallocate(diag%monin_obukhov_length)
  deallocate(diag%power_flux_density_sensible)
  deallocate(diag%power_flux_density_latent)
  deallocate(diag%sfc_sw_in)
  deallocate(diag%sfc_lw_out)
  deallocate(diag%roughness_velocity)
  ! type containing irreversible quantities
  deallocate(irrev%tke)
  deallocate(irrev%mom_diff_tend_x)
  deallocate(irrev%mom_diff_tend_y)
  deallocate(irrev%mom_diff_tend_z)
  deallocate(irrev%heating_diss)
  deallocate(irrev%mass_source_rates)
  deallocate(irrev%heat_source_rates)
  write(*,*) "... finished."
  
end program control








