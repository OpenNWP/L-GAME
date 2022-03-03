! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

program control

  ! This controls the model run from the beginning to the end.

  use run_nml,                   only: run_nml_setup,run_span_hr,dtime, &
                                       t_init,nlins,ncols,nlays,lrestart, &
                                       lideal
  use io_nml,                    only: io_nml_setup,dt_write
  use constituents_nml,          only: constituents_nml_setup,no_of_condensed_constituents,no_of_constituents, &
                                       snow_velocity,rain_velocity,cloud_droplets_velocity
  use diff_nml,                  only: diff_nml_setup
  use surface_nml,               only: surface_nml_setup,nsoillays
  use definitions,               only: t_grid,t_state,wp,t_diag,t_tend,t_bc,t_irrev
  use grid_generator,            only: grid_setup,bg_setup
  use set_initial_state,         only: restart,ideal_init
  use write_out,                 only: write_output
  use manage_rkhevi,             only: rkhevi
  use linear_combine_two_states, only: lin_combination,interpolation_t
  use bc_nml,                    only: bc_nml_setup,lperiodic
  use rad_nml,                   only: rad_nml_setup,lrad,dtime_rad
  use manage_radiation_calls,    only: call_radiation
  use boundaries,                only: setup_bc_factor
  use radiation,                 only: radiation_init
  use derived_quantities,        only: temperature_diagnostics
  
  implicit none

  ! local variables
  integer           :: timestep_counter                ! counter of the timestep
  real(wp)          :: t_0,run_span,t_write            ! time information
  type(t_grid)      :: grid                            ! grid properties
  type(t_state)     :: state_old,state_new,state_write ! states at different time steps
  type(t_diag)      :: diag                            ! diagnostic quantities
  type(t_tend)      :: tend                            ! state containing the tendency
  type(t_bc)        :: bc                              ! boundary conditions
  type(t_irrev)     :: irrev                           ! irreversible quantities
  real(wp)          :: normal_dist_min_vert            ! minimum vertical gridpoint distance
  logical           :: lrad_update                     ! radiation update switch
  real(wp)          :: t_rad_update                    ! radiation update time
  character(len=82) :: stars                           ! character containing stars

  stars = "**********************************************************************************"
  write(*,*) stars
  write(*,*) "*                                                                                *"
  write(*,*) "*                                 This is L-GAME                                 *"
  write(*,*) "*                  Regional Geophysical Fluids Modeling Framework                *"
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
  if (dtime > dt_write) then
    write(*,*) "Error: It must be dtime <= dt_write."
    call exit(1)
  endif
  
  ! allocating memory
  write(*,*) "Allocating memory ..."
  allocate(grid%lat_scalar(nlins))
  grid%lat_scalar = 0._wp
  allocate(grid%lon_scalar(ncols))
  grid%lon_scalar = 0._wp
  allocate(grid%lat_geo_scalar(nlins,ncols))
  grid%lat_geo_scalar = 0._wp
  allocate(grid%lon_geo_scalar(nlins,ncols))
  grid%lon_geo_scalar = 0._wp
  allocate(grid%lat_geo_u(nlins,ncols+1))
  grid%lat_geo_u = 0._wp
  allocate(grid%lon_geo_u(nlins,ncols+1))
  grid%lon_geo_u = 0._wp
  allocate(grid%dir_geo_u(nlins,ncols+1))
  grid%dir_geo_u = 0._wp
  allocate(grid%lat_geo_v(nlins+1,ncols))
  grid%lat_geo_v = 0._wp
  allocate(grid%lon_geo_v(nlins+1,ncols))
  grid%lon_geo_v = 0._wp
  allocate(grid%dir_geo_v(nlins+1,ncols))
  grid%dir_geo_v = 0._wp
  allocate(grid%dir_geo_u_scalar(nlins,ncols))
  grid%dir_geo_u_scalar = 0._wp
  allocate(grid%z_scalar(nlins,ncols,nlays))
  grid%z_scalar = 0._wp
  allocate(grid%dx(nlins,ncols+1,nlays))
  grid%dx = 0._wp
  allocate(grid%dy(nlins+1,ncols,nlays))
  grid%dy = 0._wp
  allocate(grid%dz(nlins,ncols,nlays+1))
  grid%dz = 0._wp
  allocate(grid%dx_dual(nlins+1,ncols,nlays))
  grid%dx_dual = 0._wp
  allocate(grid%dy_dual(nlins,ncols+1,nlays))
  grid%dy_dual = 0._wp
  allocate(grid%z_u(nlins,ncols+1,nlays))
  grid%z_u = 0._wp
  allocate(grid%z_v(nlins+1,ncols,nlays))
  grid%z_v = 0._wp
  allocate(grid%z_w(nlins,ncols,nlays+1))
  grid%z_w = 0._wp
  allocate(grid%slope_x(nlins,ncols+1,nlays))
  grid%slope_x = 0._wp
  allocate(grid%slope_y(nlins+1,ncols,nlays))
  grid%slope_y = 0._wp
  allocate(grid%volume(nlins,ncols,nlays))
  grid%volume = 0._wp
  allocate(grid%area_x(nlins,ncols+1,nlays))
  grid%area_x = 0._wp
  allocate(grid%area_y(nlins+1,ncols,nlays))
  grid%area_y = 0._wp
  allocate(grid%area_z(nlins,ncols,nlays+1))
  grid%area_z = 0._wp
  allocate(grid%area_dual_x(nlins+1,ncols,nlays+1))
  grid%area_dual_x = 0._wp
  allocate(grid%area_dual_y(nlins,ncols+1,nlays+1))
  grid%area_dual_y = 0._wp
  allocate(grid%area_dual_z(nlins+1,ncols+1,nlays))
  grid%area_dual_z = 0._wp
  allocate(grid%z_area_dual_z(nlins+1,ncols+1,nlays))
  grid%z_area_dual_z = 0._wp
  allocate(grid%inner_product_weights(nlins,ncols,nlays,6))
  grid%inner_product_weights = 0._wp
  allocate(grid%fvec_x(nlins+1,ncols))
  grid%fvec_x = 0._wp
  allocate(grid%fvec_y(nlins,ncols+1))
  grid%fvec_y = 0._wp
  allocate(grid%fvec_z(nlins+1,ncols+1))
  grid%fvec_z = 0._wp
  allocate(grid%trsk_weights_u(nlins,6))
  grid%trsk_weights_u = 0._wp
  allocate(grid%trsk_weights_v(nlins+1,4))
  grid%trsk_weights_v = 0._wp
  allocate(grid%theta_bg(nlins,ncols,nlays))
  grid%theta_bg = 0._wp
  allocate(grid%exner_bg(nlins,ncols,nlays))
  grid%exner_bg = 0._wp
  allocate(grid%exner_bg_grad_u(nlins,ncols+1,nlays))
  grid%exner_bg_grad_u = 0._wp
  allocate(grid%exner_bg_grad_v(nlins+1,ncols,nlays))
  grid%exner_bg_grad_v = 0._wp
  allocate(grid%exner_bg_grad_w(nlins,ncols,nlays+1))
  grid%exner_bg_grad_w = 0._wp
  allocate(grid%sfc_albedo(nlins,ncols))
  grid%sfc_albedo = 0._wp
  allocate(grid%sfc_rho_c(nlins,ncols))
  grid%sfc_rho_c = 0._wp
  allocate(grid%t_conduc_soil(nlins,ncols))
  grid%t_conduc_soil = 0._wp
  allocate(grid%roughness_length(nlins,ncols))
  grid%roughness_length = 0._wp
  allocate(grid%is_land(nlins,ncols))
  grid%is_land = 0
  allocate(grid%z_soil_interface(nsoillays+1))
  grid%z_soil_interface = 0._wp
  allocate(grid%z_soil_center(nsoillays))
  grid%z_soil_center = 0._wp
  allocate(grid%t_const_soil(nlins,ncols))
  grid%t_const_soil = 0._wp
  ! state at the old time step
  allocate(state_old%rho(nlins,ncols,nlays,no_of_constituents))
  state_old%rho = 0._wp
  allocate(state_old%rhotheta(nlins,ncols,nlays))
  state_old%rhotheta = 0._wp
  allocate(state_old%theta_pert(nlins,ncols,nlays))
  state_old%theta_pert = 0._wp
  allocate(state_old%exner_pert(nlins,ncols,nlays))
  state_old%exner_pert = 0._wp
  allocate(state_old%condensed_rho_t(nlins,ncols,nlays,no_of_condensed_constituents))
  state_old%condensed_rho_t = 0._wp
  allocate(state_old%wind_u(nlins,ncols+1,nlays))
  state_old%wind_u = 0._wp
  allocate(state_old%wind_v(nlins+1,ncols,nlays))
  state_old%wind_v = 0._wp
  allocate(state_old%wind_w(nlins,ncols,nlays+1))
  state_old%wind_w = 0._wp
  allocate(state_old%temperature_soil(nlins,ncols,nsoillays))
  state_old%temperature_soil = 0._wp
  ! state at the new time step
  allocate(state_new%rho(nlins,ncols,nlays,no_of_constituents))
  state_new%rho = 0._wp
  allocate(state_new%rhotheta(nlins,ncols,nlays))
  state_new%rhotheta = 0._wp
  allocate(state_new%theta_pert(nlins,ncols,nlays))
  state_new%theta_pert = 0._wp
  allocate(state_new%exner_pert(nlins,ncols,nlays))
  state_new%exner_pert = 0._wp
  allocate(state_new%condensed_rho_t(nlins,ncols,nlays,no_of_condensed_constituents))
  state_new%condensed_rho_t = 0._wp
  allocate(state_new%wind_u(nlins,ncols+1,nlays))
  state_new%wind_u = 0._wp
  allocate(state_new%wind_v(nlins+1,ncols,nlays))
  state_new%wind_v = 0._wp
  allocate(state_new%wind_w(nlins,ncols,nlays+1))
  state_new%wind_w = 0._wp
  allocate(state_new%temperature_soil(nlins,ncols,nsoillays))
  state_new%temperature_soil = 0._wp
  ! state containing the tendency
  allocate(tend%rho(nlins,ncols,nlays,no_of_constituents))
  tend%rho = 0._wp
  allocate(tend%rhotheta(nlins,ncols,nlays))
  tend%rhotheta = 0._wp
  allocate(tend%wind_u(nlins,ncols+1,nlays))
  tend%wind_u = 0._wp
  allocate(tend%wind_v(nlins+1,ncols,nlays))
  tend%wind_v = 0._wp
  allocate(tend%wind_w(nlins,ncols,nlays+1))
  tend%wind_w = 0._wp
  allocate(tend%condensed_rho_t(nlins,ncols,nlays,no_of_condensed_constituents))
  tend%condensed_rho_t = 0._wp
  ! state containing the tendency of the boundary conditions
  allocate(bc%rho_old(nlins,ncols,nlays,no_of_constituents))
  bc%rho_old = 0._wp
  allocate(bc%rhotheta_old(nlins,ncols,nlays))
  bc%rhotheta_old = 0._wp
  allocate(bc%wind_u_old(nlins,ncols+1,nlays))
  bc%wind_u_old = 0._wp
  allocate(bc%wind_v_old(nlins+1,ncols,nlays))
  bc%wind_v_old = 0._wp
  allocate(bc%wind_w_old(nlins,ncols,nlays+1))
  bc%wind_w_old = 0._wp
  allocate(bc%condensed_rho_t_old(nlins,ncols,nlays,no_of_condensed_constituents))
  bc%condensed_rho_t_old = 0._wp
  allocate(bc%rho_new(nlins,ncols,nlays,no_of_constituents))
  bc%rho_new = 0._wp
  allocate(bc%rhotheta_new(nlins,ncols,nlays))
  bc%rhotheta_new = 0._wp
  allocate(bc%wind_u_new(nlins,ncols+1,nlays))
  bc%wind_u_new = 0._wp
  allocate(bc%wind_v_new(nlins+1,ncols,nlays))
  bc%wind_v_new = 0._wp
  allocate(bc%wind_w_new(nlins,ncols,nlays+1))
  bc%wind_w_new = 0._wp
  allocate(bc%condensed_rho_t_new(nlins,ncols,nlays,no_of_condensed_constituents))
  bc%condensed_rho_t_new = 0._wp
  allocate(bc%scalar_bc_factor(nlins,ncols))
  bc%scalar_bc_factor = 0._wp
  allocate(bc%u_bc_factor(nlins,ncols+1))
  bc%u_bc_factor = 0._wp
  allocate(bc%v_bc_factor(nlins+1,ncols))
  bc%v_bc_factor = 0._wp
  ! state to be written out
  allocate(state_write%rho(nlins,ncols,nlays,no_of_constituents))
  state_write%rho = 0._wp
  allocate(state_write%rhotheta(nlins,ncols,nlays))
  state_write%rhotheta = 0._wp
  allocate(state_write%theta_pert(nlins,ncols,nlays))
  state_write%theta_pert = 0._wp
  allocate(state_write%exner_pert(nlins,ncols,nlays))
  state_write%exner_pert = 0._wp
  allocate(state_write%condensed_rho_t(nlins,ncols,nlays,no_of_condensed_constituents))
  state_write%condensed_rho_t = 0._wp
  allocate(state_write%wind_u(nlins,ncols+1,nlays))
  state_write%wind_u = 0._wp
  allocate(state_write%wind_v(nlins+1,ncols,nlays))
  state_write%wind_v = 0._wp
  allocate(state_write%wind_w(nlins,ncols,nlays+1))
  state_write%wind_w = 0._wp
  allocate(state_write%temperature_soil(nlins,ncols,nsoillays))
  state_write%temperature_soil = 0._wp
  ! type containing diagnostic quantities
  allocate(diag%v_squared(nlins,ncols,nlays))
  diag%v_squared = 0._wp
  allocate(diag%v_squared_grad_x(nlins,ncols+1,nlays))
  diag%v_squared_grad_x = 0._wp
  allocate(diag%v_squared_grad_y(nlins+1,ncols,nlays))
  diag%v_squared_grad_y = 0._wp
  allocate(diag%v_squared_grad_z(nlins,ncols,nlays+1))
  diag%v_squared_grad_z = 0._wp
  allocate(diag%p_grad_acc_neg_l_u(nlins,ncols+1,nlays))
  diag%p_grad_acc_neg_l_u = 0._wp
  allocate(diag%p_grad_acc_neg_l_v(nlins+1,ncols,nlays))
  diag%p_grad_acc_neg_l_v = 0._wp
  allocate(diag%p_grad_acc_neg_l_w(nlins,ncols,nlays+1))
  diag%p_grad_acc_neg_l_w = 0._wp
  allocate(diag%p_grad_acc_neg_nl_u(nlins,ncols+1,nlays))
  diag%p_grad_acc_neg_nl_u = 0._wp
  allocate(diag%p_grad_acc_neg_nl_v(nlins+1,ncols,nlays))
  diag%p_grad_acc_neg_nl_v = 0._wp
  allocate(diag%p_grad_acc_neg_nl_w(nlins,ncols,nlays+1))
  diag%p_grad_acc_neg_nl_w = 0._wp
  allocate(diag%p_grad_acc_old_u(nlins,ncols+1,nlays))
  diag%p_grad_acc_old_u = 0._wp
  allocate(diag%p_grad_acc_old_v(nlins+1,ncols,nlays))
  diag%p_grad_acc_old_v = 0._wp
  allocate(diag%p_grad_acc_old_w(nlins,ncols,nlays+1))
  diag%p_grad_acc_old_w = 0._wp
  allocate(diag%pot_vort_tend_x(nlins,ncols+1,nlays))
  diag%pot_vort_tend_x = 0._wp
  allocate(diag%pot_vort_tend_y(nlins+1,ncols,nlays))
  diag%pot_vort_tend_y = 0._wp
  allocate(diag%pot_vort_tend_z(nlins,ncols,nlays+1))
  diag%pot_vort_tend_z = 0._wp
  allocate(diag%scalar_placeholder(nlins,ncols,nlays))
  diag%scalar_placeholder = 0._wp
  allocate(diag%temperature_gas(nlins,ncols,nlays))
  diag%temperature_gas = 0._wp
  allocate(diag%u_placeholder(nlins,ncols+1,nlays))
  diag%u_placeholder = 0._wp
  allocate(diag%v_placeholder(nlins+1,ncols,nlays))
  diag%v_placeholder = 0._wp
  allocate(diag%w_placeholder(nlins,ncols,nlays+1))
  diag%w_placeholder = 0._wp
  allocate(diag%u_10(nlins,ncols))
  diag%u_10 = 0._wp
  allocate(diag%v_10(nlins,ncols))
  diag%v_10 = 0._wp
  allocate(diag%gust(nlins,ncols))
  diag%gust = 0._wp
  allocate(diag%mslp(nlins,ncols))
  diag%mslp = 0._wp
  allocate(diag%t_2(nlins,ncols))
  diag%t_2 = 0._wp
  allocate(diag%zeta_x(nlins+1,ncols,nlays+1))
  diag%zeta_x = 0._wp
  allocate(diag%zeta_y(nlins,ncols+1,nlays+1))
  diag%zeta_y = 0._wp
  allocate(diag%zeta_z(nlins+1,ncols+1,nlays))
  diag%zeta_z = 0._wp
  allocate(diag%eta_x(nlins+1,ncols,nlays+1))
  diag%eta_x = 0._wp
  allocate(diag%eta_y(nlins,ncols+1,nlays+1))
  diag%eta_y = 0._wp
  allocate(diag%eta_z(nlins+1,ncols+1,nlays))
  diag%eta_z = 0._wp
  allocate(diag%radiation_tendency(nlins,ncols,nlays))
  diag%radiation_tendency = 0._wp
  allocate(diag%scalar_flux_resistance(nlins,ncols))
  diag%scalar_flux_resistance = 0._wp
  allocate(diag%monin_obukhov_length(nlins,ncols))
  diag%monin_obukhov_length = 0._wp
  allocate(diag%power_flux_density_sensible(nlins,ncols))
  diag%power_flux_density_sensible = 0._wp
  allocate(diag%power_flux_density_latent(nlins,ncols))
  diag%power_flux_density_latent = 0._wp
  allocate(diag%sfc_sw_in(nlins,ncols))
  diag%sfc_sw_in = 0._wp
  allocate(diag%sfc_lw_out(nlins,ncols))
  diag%sfc_lw_out = 0._wp
  allocate(diag%roughness_velocity(nlins,ncols))
  diag%roughness_velocity = 0._wp
  allocate(diag%flux_density_u(nlins,ncols+1,nlays))
  diag%flux_density_u = 0._wp
  allocate(diag%flux_density_v(nlins+1,ncols,nlays))
  diag%flux_density_v = 0._wp
  allocate(diag%flux_density_w(nlins,ncols,nlays+1))
  diag%flux_density_w = 0._wp
  allocate(diag%flux_density_div(nlins,ncols,nlays))
  diag%flux_density_div = 0._wp
  allocate(diag%du_dz(nlins,ncols+1,nlays+1))
  diag%du_dz = 0._wp
  allocate(diag%dv_dz(nlins+1,ncols,nlays+1))
  diag%dv_dz = 0._wp
  ! type containing irreversible quantities
  allocate(irrev%tke(nlins,ncols,nlays))
  irrev%tke = 0._wp
  allocate(irrev%viscosity_molecular(nlins,ncols,nlays))
  irrev%viscosity_molecular = 0._wp
  allocate(irrev%viscosity_coeff_div(nlins,ncols,nlays))
  irrev%viscosity_coeff_div = 0._wp
  allocate(irrev%viscosity_coeff_curl(nlins,ncols,nlays))
  irrev%viscosity_coeff_curl = 0._wp
  allocate(irrev%viscosity_coeff_curl_dual(nlins+1,ncols+1,nlays))
  irrev%viscosity_coeff_curl_dual = 0._wp
  allocate(irrev%vert_hor_viscosity_u(nlins,ncols+1,nlays+1))
  irrev%vert_hor_viscosity_u = 0._wp
  allocate(irrev%vert_hor_viscosity_v(nlins+1,ncols,nlays+1))
  irrev%vert_hor_viscosity_v = 0._wp
  allocate(irrev%scalar_diff_coeff_h(nlins,ncols,nlays))
  irrev%scalar_diff_coeff_h = 0._wp
  allocate(irrev%scalar_diff_coeff_v(nlins,ncols,nlays))
  irrev%scalar_diff_coeff_v = 0._wp
  allocate(irrev%mom_diff_tend_x(nlins,ncols+1,nlays))
  irrev%mom_diff_tend_x = 0._wp
  allocate(irrev%mom_diff_tend_y(nlins+1,ncols,nlays))
  irrev%mom_diff_tend_y = 0._wp
  allocate(irrev%mom_diff_tend_z(nlins,ncols,nlays+1))
  irrev%mom_diff_tend_z = 0._wp
  allocate(irrev%heating_diss(nlins,ncols,nlays))
  irrev%heating_diss = 0._wp
  allocate(irrev%mass_source_rates(nlins,ncols,nlays,no_of_condensed_constituents+1))
  irrev%mass_source_rates = 0._wp
  allocate(irrev%heat_source_rates(nlins,ncols,nlays,no_of_condensed_constituents))
  irrev%heat_source_rates = 0._wp
  allocate(irrev%temp_diff_heating(nlins,ncols,nlays))
  irrev%temp_diff_heating = 0._wp
  write(*,*) "... finished."
  
  ! firstly, the grid generator needs to be called to calculate the grid properties
  write(*,*) "Setting up the grid ..."
  call grid_setup(grid)
  write(*,*) "... grid set up."
  
  if (.not. lperiodic) then
    call setup_bc_factor(bc)
  endif
  
  ! limitting the hydrometeor sedimentation velocities for stability reasons
  normal_dist_min_vert = minval(grid%dz(:,:,nlays+1))
  rain_velocity = min(0.8_wp*normal_dist_min_vert/dtime,rain_velocity)
  snow_velocity = min(0.8_wp*normal_dist_min_vert/dtime,snow_velocity)
  write(*,*) "Rain falling velocity set to", rain_velocity, "m/s."
  write(*,*) "Snow falling velocity set to", snow_velocity, "m/s."
  
  ! maximum horizontal diffusion coefficient  
  irrev%max_diff_h_coeff_turb = 0.125_wp*grid%mean_velocity_area/dtime
  ! setting up the background state
  call bg_setup(grid)
  
  ! setting the initial state
  write(*,*) "Setting the initial state..."
  if (lrestart) then
    call restart(state_old,grid)
  elseif (lideal) then
    call ideal_init(state_old,diag,grid)
  endif
  write(*,*) "... initial state set."
  
  ! updating radiation for the first time if nescessary
  t_0 = t_init
  t_rad_update = t_0
  ! calculating the temperature of the gas phase
  call temperature_diagnostics(state_old,diag,grid)
  if (lrad) then
    call radiation_init()
    call call_radiation(state_old,grid,diag,t_0)
  endif
  ! setting the next time for the radiation update
  t_rad_update = t_rad_update+dtime_rad
  
  ! writing out the initial state
  call write_output(state_old,diag,0,grid)
  
  ! copying the new state to the old state
  state_new = state_old
  
  ! the loop over the time steps
  t_write = t_0 + dt_write
  run_span = 3600._wp*run_span_hr
  timestep_counter = 0
  do while (t_0<t_init+run_span+300._wp .and. run_span/=0)
    
    ! writing the new state into the old state
    call lin_combination(state_new,state_new,state_old,0._wp,1._wp,grid)

    if (lrad .and. t_0<=t_rad_update .and. t_0+dtime>=t_rad_update) then
      lrad_update = .true.
      t_rad_update = t_rad_update+dtime_rad
    else
      lrad_update = .false.
    endif

    ! this is the RKHEVI routine performing the time stepping
    call rkhevi(state_old,state_new,tend,bc,grid,diag,irrev,timestep_counter,lrad_update,t_0)
    
    ! managing the calls to the output routine
    if (t_0+dtime>=t_write) then
      call interpolation_t(state_old,state_new,state_write,t_0,t_0+dtime,t_write,grid)
      call write_output(state_write,diag,int((t_write-t_init)/60._wp),grid)
    
      t_write = t_write+dt_write
    
    endif
    
    t_0 = t_0+dtime
    timestep_counter = timestep_counter+1
    write(*,*) "Step ", timestep_counter, " completed."
    
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
  ! state containing the tendency of the boundary conditions
  deallocate(bc%rho_old)
  deallocate(bc%rhotheta_old)
  deallocate(bc%wind_u_old)
  deallocate(bc%wind_v_old)
  deallocate(bc%wind_w_old)
  deallocate(bc%condensed_rho_t_old)
  deallocate(bc%rho_new)
  deallocate(bc%rhotheta_new)
  deallocate(bc%wind_u_new)
  deallocate(bc%wind_v_new)
  deallocate(bc%wind_w_new)
  deallocate(bc%condensed_rho_t_new)
  deallocate(bc%scalar_bc_factor)
  deallocate(bc%u_bc_factor)
  deallocate(bc%v_bc_factor)
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
  deallocate(diag%temperature_gas)
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
  ! type containing irreversible quantities
  deallocate(irrev%tke)
  deallocate(irrev%viscosity_molecular)
  deallocate(irrev%viscosity_coeff_div)
  deallocate(irrev%viscosity_coeff_curl)
  deallocate(irrev%viscosity_coeff_curl_dual)
  deallocate(irrev%vert_hor_viscosity_u)
  deallocate(irrev%vert_hor_viscosity_v)
  deallocate(irrev%scalar_diff_coeff_h)
  deallocate(irrev%scalar_diff_coeff_v)
  deallocate(irrev%mom_diff_tend_x)
  deallocate(irrev%mom_diff_tend_y)
  deallocate(irrev%mom_diff_tend_z)
  deallocate(irrev%heating_diss)
  deallocate(irrev%mass_source_rates)
  deallocate(irrev%heat_source_rates)
  deallocate(irrev%temp_diff_heating)
  write(*,*) "... finished."
  write(*,*) "L-GAME over."
  
end program control







