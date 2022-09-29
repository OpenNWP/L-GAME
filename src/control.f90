! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

program control

  ! This controls the model run from the beginning to the end.

  use mo_run_nml,                   only: run_nml_setup,run_span_min,dtime,t_init,ny,nx,n_layers,lrestart, &
                                          lideal,n_levels
  use mo_io_nml,                    only: io_nml_setup,dt_write
  use mo_constituents_nml,          only: constituents_nml_setup,n_condensed_constituents,n_constituents, &
                                          snow_velocity,rain_velocity
  use mo_diff_nml,                  only: diff_nml_setup
  use mo_surface_nml,               only: surface_nml_setup,nsoillays
  use mo_definitions,               only: t_grid,t_state,wp,t_diag,t_tend,t_bc
  use mo_grid_generator,            only: grid_setup,bg_setup
  use mo_set_initial_state,         only: restart,ideal_init
  use mo_write_out,                 only: write_output
  use mo_manage_pchevi,             only: pchevi
  use mo_linear_combine_two_states, only: lin_combination
  use mo_bc_nml,                    only: bc_nml_setup,lperiodic,t_latest_bc,dtime_bc
  use mo_rad_nml,                   only: rad_nml_setup,lrad,dtime_rad
  use mo_manage_radiation_calls,    only: update_rad_fluxes
  use mo_boundaries,                only: setup_bc_factor,read_boundaries
  use mo_rrtmgp_coupler,            only: radiation_init
  use mo_derived,                   only: temperature_diagnostics
  
  implicit none

  ! local variables
  integer           :: timestep_counter                ! counter of the timestep
  real(wp)          :: t_0,run_span,t_write            ! time information
  type(t_grid)      :: grid                            ! grid properties
  type(t_state)     :: state_1,state_2,state_write ! states at different time steps
  type(t_diag)      :: diag                            ! diagnostic quantities
  type(t_tend)      :: tend                            ! state containing the tendency
  type(t_bc)        :: bc                              ! boundary conditions
  real(wp)          :: normal_dist_min_vert            ! minimum vertical gridpoint distance
  logical           :: lrad_update                     ! radiation update switch
  real(wp)          :: t_rad_update                    ! radiation update time
  character(len=82) :: stars                           ! character containing stars

  stars = "**********************************************************************************"
  write(*,*) stars
  write(*,*) "*                                                                                *"
  write(*,*) "*                                 This is L-GAME                                 *"
  write(*,*) "*               Limited-area Geophysical Fluids Modeling Framework               *"
  write(*,*) "*                                                                                *"
  write(*,*) "*                         Released under the MIT license.                        *"
  write(*,*) "*         Visit https://github.com/OpenNWP/L-GAME for more information.          *"
  write(*,*) "*                                                                                *"
  write(*,*) stars

  ! reading in all namelists so that we know what we have to do
  write(*,*) "Reading in run namelist ..."
  call run_nml_setup()
  write(*,*) "... run namelist read."
  
  write(*,*) "Reading in diff namelist ..."
  call diff_nml_setup()
  write(*,*) "... diff namelist read."
  
  write(*,*) "Reading in I/O namelist ..."
  call io_nml_setup()
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
  
  write(*,*) "Reading in radiation namelist ..."
  call rad_nml_setup()
  write(*,*) "... radiation namelist read."
  
  ! sanity check
  if (dtime>dt_write) then
    write(*,*) "Error: It must be dtime <= dt_write."
    call exit(1)
  endif
  
  ! allocating memory
  write(*,*) "Allocating memory ..."
  allocate(grid%lat_scalar(ny))
  allocate(grid%lon_scalar(nx))
  allocate(grid%lat_geo_scalar(ny,nx))
  allocate(grid%lon_geo_scalar(ny,nx))
  allocate(grid%lat_geo_u(ny,nx+1))
  allocate(grid%lon_geo_u(ny,nx+1))
  allocate(grid%dir_geo_u(ny,nx+1))
  allocate(grid%lat_geo_v(ny+1,nx))
  allocate(grid%lon_geo_v(ny+1,nx))
  allocate(grid%dir_geo_v(ny+1,nx))
  allocate(grid%dir_geo_u_scalar(ny,nx))
  allocate(grid%z_scalar(ny,nx,n_layers))
  allocate(grid%dx(ny,nx+1,n_layers))
  allocate(grid%dy(ny+1,nx,n_layers))
  allocate(grid%dz(ny,nx,n_levels))
  allocate(grid%layer_thickness(ny,nx,n_layers))
  allocate(grid%dx_dual(ny+1,nx,n_layers))
  allocate(grid%dy_dual(ny,nx+1,n_layers))
  allocate(grid%z_u(ny,nx+1,n_layers))
  allocate(grid%z_v(ny+1,nx,n_layers))
  allocate(grid%z_w(ny,nx,n_levels))
  allocate(grid%slope_x(ny,nx+1,n_layers))
  allocate(grid%slope_y(ny+1,nx,n_layers))
  allocate(grid%volume(ny,nx,n_layers))
  allocate(grid%area_x(ny,nx+1,n_layers))
  allocate(grid%area_y(ny+1,nx,n_layers))
  allocate(grid%area_z(ny,nx,n_levels))
  allocate(grid%area_dual_x(ny+1,nx,n_levels))
  allocate(grid%area_dual_y(ny,nx+1,n_levels))
  allocate(grid%area_dual_z(ny+1,nx+1,n_layers))
  allocate(grid%z_area_dual_z(ny+1,nx+1,n_layers))
  allocate(grid%inner_product_weights(ny,nx,n_layers,6))
  allocate(grid%fvec_x(ny+1,nx))
  allocate(grid%fvec_y(ny,nx+1))
  allocate(grid%fvec_z(ny+1,nx+1))
  allocate(grid%trsk_weights_u(ny,6))
  allocate(grid%trsk_weights_v(ny+1,4))
  allocate(grid%theta_v_bg(ny,nx,n_layers))
  allocate(grid%exner_bg(ny,nx,n_layers))
  allocate(grid%exner_bg_grad_u(ny,nx+1,n_layers))
  allocate(grid%exner_bg_grad_v(ny+1,nx,n_layers))
  allocate(grid%exner_bg_grad_w(ny,nx,n_levels))
  allocate(grid%sfc_albedo(ny,nx))
  allocate(grid%sfc_rho_c(ny,nx))
  allocate(grid%t_conduc_soil(ny,nx))
  allocate(grid%roughness_length(ny,nx))
  allocate(grid%is_land(ny,nx))
  grid%is_land = 0
  allocate(grid%z_soil_interface(nsoillays+1))
  allocate(grid%z_soil_center(nsoillays))
  allocate(grid%t_const_soil(ny,nx))
  ! state at the old time step
  allocate(state_1%rho(ny,nx,n_layers,n_constituents))
  allocate(state_1%rhotheta_v(ny,nx,n_layers))
  allocate(state_1%theta_v_pert(ny,nx,n_layers))
  allocate(state_1%exner_pert(ny,nx,n_layers))
  allocate(state_1%wind_u(ny,nx+1,n_layers))
  allocate(state_1%wind_v(ny+1,nx,n_layers))
  allocate(state_1%wind_w(ny,nx,n_levels))
  allocate(state_1%temperature_soil(ny,nx,nsoillays))
  ! state at the new time step
  allocate(state_2%rho(ny,nx,n_layers,n_constituents))
  allocate(state_2%rhotheta_v(ny,nx,n_layers))
  allocate(state_2%theta_v_pert(ny,nx,n_layers))
  allocate(state_2%exner_pert(ny,nx,n_layers))
  allocate(state_2%wind_u(ny,nx+1,n_layers))
  allocate(state_2%wind_v(ny+1,nx,n_layers))
  allocate(state_2%wind_w(ny,nx,n_levels))
  allocate(state_2%temperature_soil(ny,nx,nsoillays))
  ! state containing the tendency
  allocate(tend%rho(ny,nx,n_layers,n_constituents))
  allocate(tend%rhotheta_v(ny,nx,n_layers))
  allocate(tend%wind_u(ny,nx+1,n_layers))
  allocate(tend%wind_v(ny+1,nx,n_layers))
  allocate(tend%wind_w(ny,nx,n_levels))
  ! state containing the tendency of the boundary conditions
  allocate(bc%rho(ny,nx,n_layers,n_constituents,2))
  allocate(bc%rhotheta_v(ny,nx,n_layers,2))
  allocate(bc%wind_u(ny,nx+1,n_layers,2))
  allocate(bc%wind_v(ny+1,nx,n_layers,2))
  allocate(bc%wind_w(ny,nx,n_levels,2))
  allocate(bc%scalar_bc_factor(ny,nx))
  allocate(bc%u_bc_factor(ny,nx+1))
  allocate(bc%v_bc_factor(ny+1,nx))
  bc%index_old = 1
  bc%index_new = 2
  ! state to be written out
  allocate(state_write%rho(ny,nx,n_layers,n_constituents))
  allocate(state_write%rhotheta_v(ny,nx,n_layers))
  allocate(state_write%theta_v_pert(ny,nx,n_layers))
  allocate(state_write%exner_pert(ny,nx,n_layers))
  allocate(state_write%wind_u(ny,nx+1,n_layers))
  allocate(state_write%wind_v(ny+1,nx,n_layers))
  allocate(state_write%wind_w(ny,nx,n_levels))
  allocate(state_write%temperature_soil(ny,nx,nsoillays))
  ! type containing diagnostic quantities
  allocate(diag%v_squared(ny,nx,n_layers))
  allocate(diag%v_squared_grad_x(ny,nx+1,n_layers))
  allocate(diag%v_squared_grad_y(ny+1,nx,n_layers))
  allocate(diag%v_squared_grad_z(ny,nx,n_levels))
  allocate(diag%p_grad_acc_neg_l_u(ny,nx+1,n_layers))
  allocate(diag%p_grad_acc_neg_l_v(ny+1,nx,n_layers))
  allocate(diag%p_grad_acc_neg_l_w(ny,nx,n_levels))
  allocate(diag%p_grad_acc_neg_nl_u(ny,nx+1,n_layers))
  allocate(diag%p_grad_acc_neg_nl_v(ny+1,nx,n_layers))
  allocate(diag%p_grad_acc_neg_nl_w(ny,nx,n_levels))
  allocate(diag%p_grad_acc_old_u(ny,nx+1,n_layers))
  allocate(diag%p_grad_acc_old_v(ny+1,nx,n_layers))
  allocate(diag%p_grad_acc_old_w(ny,nx,n_levels))
  allocate(diag%pot_vort_tend_x(ny,nx+1,n_layers))
  allocate(diag%pot_vort_tend_y(ny+1,nx,n_layers))
  allocate(diag%pot_vort_tend_z(ny,nx,n_levels))
  allocate(diag%scalar_placeholder(ny,nx,n_layers))
  allocate(diag%temperature(ny,nx,n_layers))
  allocate(diag%u_placeholder(ny,nx+1,n_layers))
  allocate(diag%v_placeholder(ny+1,nx,n_layers))
  allocate(diag%w_placeholder(ny,nx,n_levels))
  allocate(diag%u_10(ny,nx))
  allocate(diag%v_10(ny,nx))
  allocate(diag%gust(ny,nx))
  allocate(diag%mslp(ny,nx))
  allocate(diag%t_2(ny,nx))
  allocate(diag%zeta_x(ny+1,nx,n_levels))
  allocate(diag%zeta_y(ny,nx+1,n_levels))
  allocate(diag%zeta_z(ny+1,nx+1,n_layers))
  allocate(diag%eta_x(ny+1,nx,n_levels))
  allocate(diag%eta_y(ny,nx+1,n_levels))
  allocate(diag%eta_z(ny+1,nx+1,n_layers))
  allocate(diag%radiation_tendency(ny,nx,n_layers))
  allocate(diag%scalar_flux_resistance(ny,nx))
  allocate(diag%monin_obukhov_length(ny,nx))
  allocate(diag%power_flux_density_sensible(ny,nx))
  allocate(diag%power_flux_density_latent(ny,nx))
  allocate(diag%sfc_sw_in(ny,nx))
  allocate(diag%sfc_lw_out(ny,nx))
  allocate(diag%roughness_velocity(ny,nx))
  allocate(diag%flux_density_u(ny,nx+1,n_layers))
  allocate(diag%flux_density_v(ny+1,nx,n_layers))
  allocate(diag%flux_density_w(ny,nx,n_levels))
  allocate(diag%flux_density_div(ny,nx,n_layers))
  allocate(diag%du_dz(ny,nx+1,n_levels))
  allocate(diag%dv_dz(ny+1,nx,n_levels))
  allocate(diag%n_squared(ny,nx,n_layers))
  allocate(diag%tke(ny,nx,n_layers))
  allocate(diag%viscosity_molecular(ny,nx,n_layers))
  allocate(diag%viscosity_coeff_div(ny,nx,n_layers))
  allocate(diag%viscosity_coeff_curl(ny,nx,n_layers))
  allocate(diag%viscosity_coeff_curl_dual(ny+1,nx+1,n_layers))
  allocate(diag%vert_hor_viscosity_u(ny,nx+1,n_levels))
  allocate(diag%vert_hor_viscosity_v(ny+1,nx,n_levels))
  allocate(diag%scalar_diff_coeff_h(ny,nx,n_layers))
  allocate(diag%scalar_diff_coeff_v(ny,nx,n_layers))
  allocate(diag%pressure_gradient_decel_factor(ny,nx,n_layers))
  allocate(diag%mom_diff_tend_x(ny,nx+1,n_layers))
  allocate(diag%mom_diff_tend_y(ny+1,nx,n_layers))
  allocate(diag%mom_diff_tend_z(ny,nx,n_levels))
  allocate(diag%heating_diss(ny,nx,n_layers))
  allocate(diag%phase_trans_rates(ny,nx,n_layers,n_condensed_constituents+1))
  allocate(diag%phase_trans_heating_rate(ny,nx,n_layers))
  allocate(diag%temp_diff_heating(ny,nx,n_layers))
  allocate(diag%condensates_sediment_heat(ny,nx,n_layers))
  allocate(diag%mass_diff_tendency(ny,nx,n_layers,n_constituents))
  write(*,*) "... finished."
  
  ! firstly, the grid generator needs to be called to calculate the grid properties
  write(*,*) "Setting up the grid ..."
  call grid_setup(grid)
  write(*,*) "... grid set up."
  
  if (.not. lperiodic) then
    call setup_bc_factor(bc)
  endif
  
  ! limitting the hydrometeor sedimentation velocities for stability reasons
  normal_dist_min_vert = minval(grid%dz(:,:,n_levels))
  rain_velocity = min(0.8_wp*normal_dist_min_vert/dtime,rain_velocity)
  snow_velocity = min(0.8_wp*normal_dist_min_vert/dtime,snow_velocity)
  write(*,*) "Rain falling velocity set to", rain_velocity, "m/s."
  write(*,*) "Snow falling velocity set to", snow_velocity, "m/s."
  
  ! setting up the background state
  call bg_setup(grid)
  
  ! setting the initial state
  write(*,*) "Setting the initial state..."
  if (lrestart) then
    call restart(state_1,grid)
    call read_boundaries(bc,real(t_init,wp),1)
    call read_boundaries(bc,real(t_init+dtime_bc,wp),2)
    t_latest_bc = t_0
  elseif (lideal) then
    call ideal_init(state_1,diag,grid)
  endif
  write(*,*) "... initial state set."
  
  ! updating radiation for the first time if nescessary
  t_0 = t_init
  t_rad_update = t_0
  ! calculating the temperature of the gas phase
  call temperature_diagnostics(state_1,diag,grid)
  if (lrad) then
    call radiation_init()
    call update_rad_fluxes(state_1,grid,diag,t_0)
  endif
  ! setting the next time for the radiation update
  t_rad_update = t_rad_update+dtime_rad
  
  ! writing out the initial state
  call write_output(state_1,diag,0,grid)
  
  ! copying the first state to the second state
  state_2 = state_1
  
  ! the loop over the time steps
  t_write = t_0 + dt_write
  run_span = 3600._wp*run_span_min
  timestep_counter = 0
  do while (t_0<t_init+run_span+300._wp .and. run_span/=0)
    
    ! writing the new state into the old state
    call lin_combination(state_2,state_2,state_1,0._wp,1._wp,grid)

    if (lrad .and. t_0<=t_rad_update .and. t_0+dtime>=t_rad_update) then
      lrad_update = .true.
      t_rad_update = t_rad_update+dtime_rad
    else
      lrad_update = .false.
    endif

    ! this is the RKHEVI routine performing the time stepping
    if (mod(timestep_counter,2)==0) then
      call pchevi(state_1,state_2,tend,bc,grid,diag,timestep_counter,lrad_update,t_0)
    else
      call pchevi(state_2,state_1,tend,bc,grid,diag,timestep_counter,lrad_update,t_0)
    endif
    
    ! managing the calls to the output routine
    if (t_0+dtime>=t_write) then
      call lin_combination(state_1,state_2,state_write,1._wp-(t_write-t_0)/dtime,(t_write-t_0)/dtime,grid)
      call write_output(state_write,diag,int((t_write-t_init)/60._wp),grid)
    
      t_write = t_write+dt_write
    
    endif
    
    t_0 = t_0+dtime
    timestep_counter = timestep_counter+1
    write(*,*) "Step", timestep_counter, " completed."
    
  enddo
  
  ! deallocating the memory
  write(*,*) "Deallocating memory ..."
  deallocate(grid%lat_scalar)
  deallocate(grid%lon_scalar)
  deallocate(grid%lat_geo_scalar)
  deallocate(grid%lon_geo_scalar)
  deallocate(grid%lat_geo_u)
  deallocate(grid%lon_geo_u)
  deallocate(grid%dir_geo_u)
  deallocate(grid%lat_geo_v)
  deallocate(grid%lon_geo_v)
  deallocate(grid%dir_geo_v)
  deallocate(grid%dir_geo_u_scalar)
  deallocate(grid%z_scalar)
  deallocate(grid%dy)
  deallocate(grid%dx)
  deallocate(grid%dz)
  deallocate(grid%layer_thickness)
  deallocate(grid%dx_dual)
  deallocate(grid%dy_dual)
  deallocate(grid%z_u)
  deallocate(grid%z_v)
  deallocate(grid%z_w)
  deallocate(grid%slope_x)
  deallocate(grid%slope_y)
  deallocate(grid%volume)
  deallocate(grid%area_x)
  deallocate(grid%area_y)
  deallocate(grid%area_z)
  deallocate(grid%area_dual_x)
  deallocate(grid%area_dual_y)
  deallocate(grid%area_dual_z)
  deallocate(grid%z_area_dual_z)
  deallocate(grid%inner_product_weights)
  deallocate(grid%fvec_x)
  deallocate(grid%fvec_y)
  deallocate(grid%fvec_z)
  deallocate(grid%trsk_weights_u)
  deallocate(grid%trsk_weights_v)
  deallocate(grid%theta_v_bg)
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
  deallocate(state_1%rho)
  deallocate(state_1%rhotheta_v)
  deallocate(state_1%theta_v_pert)
  deallocate(state_1%exner_pert)
  deallocate(state_1%wind_u)
  deallocate(state_1%wind_v)
  deallocate(state_1%wind_w)
  deallocate(state_1%temperature_soil)
  ! state at the new time step
  deallocate(state_2%rho)
  deallocate(state_2%rhotheta_v)
  deallocate(state_2%theta_v_pert)
  deallocate(state_2%exner_pert)
  deallocate(state_2%wind_u)
  deallocate(state_2%wind_v)
  deallocate(state_2%wind_w)
  deallocate(state_2%temperature_soil)
  ! state containing the tendency
  deallocate(tend%rho)
  deallocate(tend%rhotheta_v)
  deallocate(tend%wind_u)
  deallocate(tend%wind_v)
  deallocate(tend%wind_w)
  ! state containing the tendency of the boundary conditions
  deallocate(bc%rho)
  deallocate(bc%rhotheta_v)
  deallocate(bc%wind_u)
  deallocate(bc%wind_v)
  deallocate(bc%wind_w)
  deallocate(bc%scalar_bc_factor)
  deallocate(bc%u_bc_factor)
  deallocate(bc%v_bc_factor)
  ! state to be written out
  deallocate(state_write%rho)
  deallocate(state_write%rhotheta_v)
  deallocate(state_write%theta_v_pert)
  deallocate(state_write%exner_pert)
  deallocate(state_write%wind_u)
  deallocate(state_write%wind_v)
  deallocate(state_write%wind_w)
  deallocate(state_write%temperature_soil)
  ! type containing diagnostic quantities
  deallocate(diag%v_squared)
  deallocate(diag%p_grad_acc_neg_l_u)
  deallocate(diag%p_grad_acc_neg_l_v)
  deallocate(diag%p_grad_acc_neg_l_w)
  deallocate(diag%p_grad_acc_neg_nl_u)
  deallocate(diag%p_grad_acc_neg_nl_v)
  deallocate(diag%p_grad_acc_neg_nl_w)
  deallocate(diag%p_grad_acc_old_u)
  deallocate(diag%p_grad_acc_old_v)
  deallocate(diag%p_grad_acc_old_w)
  deallocate(diag%v_squared_grad_x)
  deallocate(diag%v_squared_grad_y)
  deallocate(diag%v_squared_grad_z)
  deallocate(diag%pot_vort_tend_x)
  deallocate(diag%pot_vort_tend_y)
  deallocate(diag%pot_vort_tend_z)
  deallocate(diag%scalar_placeholder)
  deallocate(diag%temperature)
  deallocate(diag%u_placeholder)
  deallocate(diag%v_placeholder)
  deallocate(diag%w_placeholder)
  deallocate(diag%u_10)
  deallocate(diag%v_10)
  deallocate(diag%gust)
  deallocate(diag%mslp)
  deallocate(diag%t_2)
  deallocate(diag%zeta_x)
  deallocate(diag%zeta_y)
  deallocate(diag%zeta_z)
  deallocate(diag%eta_x)
  deallocate(diag%eta_y)
  deallocate(diag%eta_z)
  deallocate(diag%radiation_tendency)
  deallocate(diag%scalar_flux_resistance)
  deallocate(diag%monin_obukhov_length)
  deallocate(diag%power_flux_density_sensible)
  deallocate(diag%power_flux_density_latent)
  deallocate(diag%sfc_sw_in)
  deallocate(diag%sfc_lw_out)
  deallocate(diag%roughness_velocity)
  deallocate(diag%flux_density_u)
  deallocate(diag%flux_density_v)
  deallocate(diag%flux_density_w)
  deallocate(diag%flux_density_div)
  deallocate(diag%du_dz)
  deallocate(diag%dv_dz)
  deallocate(diag%n_squared)
  deallocate(diag%tke)
  deallocate(diag%viscosity_molecular)
  deallocate(diag%viscosity_coeff_div)
  deallocate(diag%viscosity_coeff_curl)
  deallocate(diag%viscosity_coeff_curl_dual)
  deallocate(diag%vert_hor_viscosity_u)
  deallocate(diag%vert_hor_viscosity_v)
  deallocate(diag%scalar_diff_coeff_h)
  deallocate(diag%scalar_diff_coeff_v)
  deallocate(diag%mom_diff_tend_x)
  deallocate(diag%mom_diff_tend_y)
  deallocate(diag%mom_diff_tend_z)
  deallocate(diag%pressure_gradient_decel_factor)
  deallocate(diag%heating_diss)
  deallocate(diag%phase_trans_rates)
  deallocate(diag%phase_trans_heating_rate)
  deallocate(diag%temp_diff_heating)
  deallocate(diag%condensates_sediment_heat)
  deallocate(diag%mass_diff_tendency)
  write(*,*) "... finished."
  write(*,*) "L-GAME over."
  
end program control







