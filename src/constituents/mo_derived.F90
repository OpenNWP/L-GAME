! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_derived

  ! In this module more complex thermodynamic quantities are being calculated.
  
  use mo_definitions,      only: wp,t_grid,t_state,t_diag
  use mo_run_nml,          only: ny,nx,n_layers
  use mo_constants,        only: k_b,M_PI,m_d,n_a,r_d,r_v,c_d_p,c_v_p,c_d_v,c_v_v,t_0,m_v,rho_h2o,gravity
  use mo_constituents_nml, only: n_condensed_constituents,n_gaseous_constituents,n_constituents,lmoist
  use mo_dictionary,       only: saturation_pressure_over_ice,saturation_pressure_over_water,c_p_cond
  
  implicit none
  
  contains

  subroutine temperature_diagnostics(state,diag,grid)
    
    ! This function diagnoses the temperature of the gas phase.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities (needed for the background state)
    
    if (.not. lmoist) then
      !$omp parallel workshare
      diag%temperature = (grid%theta_v_bg + state%theta_v_pert)*(grid%exner_bg + state%exner_pert)
      !$omp end parallel workshare
    else
      !$omp parallel workshare
      diag%temperature = (grid%theta_v_bg + state%theta_v_pert)*(grid%exner_bg + state%exner_pert) &
      /(1._wp+state%rho(:,:,:,n_condensed_constituents+2) &
      /state%rho(:,:,:,n_condensed_constituents+1)*(m_d/m_v-1._wp))
      !$omp end parallel workshare
    endif
    
  end subroutine temperature_diagnostics

  function spec_heat_cap_diagnostics_v(state,ji,jk,jl)
    
    ! This function calculates the specific heat capacity at constant volume of the air at a certain gridpoint.
    
    type(t_state), intent(in) :: state                       ! state with which to calculate c_v
    integer,       intent(in) :: ji                          ! horizontal index of the gridpoint
    integer,       intent(in) :: jk                          ! horizontal index of the gridpoint
    integer,       intent(in) :: jl                          ! vertical index of the gridpoint
    real(wp)                  :: spec_heat_cap_diagnostics_v ! result
    
    ! local variables
    integer  :: n_relevant_constituents       ! the number of relevant constituents for this calculation
    integer  :: jc                            ! constituent index
    real(wp) :: rho_g                         ! gas density
    real(wp) :: spec_heat_capacities_v_gas(2) ! vector containing specific heat capacities
    
    n_relevant_constituents = 1
    rho_g = state%rho(ji,jk,jl,n_condensed_constituents+1)
    
    spec_heat_capacities_v_gas(1) = c_d_v
    spec_heat_capacities_v_gas(2) = c_v_v
    
    spec_heat_cap_diagnostics_v = 0._wp
    do jc=1,n_relevant_constituents
      spec_heat_cap_diagnostics_v = spec_heat_cap_diagnostics_v +  state%rho(ji,jk,jl,n_condensed_constituents+jc) &
      /rho_g*spec_heat_capacities_v_gas(jc)
    enddo
    
  end function spec_heat_cap_diagnostics_v

  function spec_heat_cap_diagnostics_p(state,ji,jk,jl)
    
    ! This function calculates the specific heat capacity at constant pressure of the air at a certain gridpoint.
    
    type(t_state), intent(in) :: state                       ! state with which to calculate c_p
    integer,       intent(in) :: ji                          ! horizontal index of the gridpoint
    integer,       intent(in) :: jk                          ! horizontal index of the gridpoint
    integer,       intent(in) :: jl                          ! vertical index of the gridpoint
    real(wp)                  :: spec_heat_cap_diagnostics_p ! result
    
    ! local variables
    integer  :: n_relevant_constituents       ! the number of relevant constituents for this calculation
    integer  :: jc                            ! constituent index
    real(wp) :: rho_g                         ! gas density
    real(wp) :: spec_heat_capacities_p_gas(2) ! vector containing specific heat capacities
    
    n_relevant_constituents = 1
    rho_g = state%rho(ji,jk,jl,n_condensed_constituents+1)
    
    spec_heat_capacities_p_gas(1) = c_d_p
    spec_heat_capacities_p_gas(2) = c_v_p
    
    spec_heat_cap_diagnostics_p = 0._wp
    do jc=1,n_relevant_constituents
      spec_heat_cap_diagnostics_p = spec_heat_cap_diagnostics_p + state%rho(ji,jk,jl,n_condensed_constituents+jc) &
      /rho_g*spec_heat_capacities_p_gas(jc)
    enddo
    
  end function spec_heat_cap_diagnostics_p
    
  function gas_constant_diagnostics(state,ji,jk,jl)
    
    ! This function calculates the specific gas constant at a certain gridpoint.
    
    type(t_state), intent(in) :: state                    ! state with which to calculate the gas constant
    integer,       intent(in) :: ji                       ! horizontal index of the gridpoint
    integer,       intent(in) :: jk                       ! horizontal index of the gridpoint
    integer,       intent(in) :: jl                       ! vertical index of the gridpoint
    real(wp)                  :: gas_constant_diagnostics ! result
    
    ! local variables
    integer  :: n_relevant_constituents   ! the number of relevant constituents for this calculation
    integer  :: jc                        ! constituent index
    real(wp) :: rho_g                     ! gas density
    real(wp) :: specific_gas_constants(2) ! vector containing specific gas constants
    
    n_relevant_constituents = 1
    rho_g = state%rho(ji,jk,jl,n_condensed_constituents+1)
    
    gas_constant_diagnostics = 0._wp
    
    specific_gas_constants(1) = r_d
    specific_gas_constants(2) = r_v
    
    do jc=1,n_relevant_constituents
      gas_constant_diagnostics = gas_constant_diagnostics + state%rho(ji,jk,jl,n_condensed_constituents+jc) &
      /rho_g*specific_gas_constants(jc)
    enddo
    
  end function gas_constant_diagnostics

  function rel_humidity(abs_humidity,temperature)
    
    ! This function returns the relative humidity as a function of the absolute humidity in kg/m^3 and the temperature in K.
    
    real(wp), intent(in) :: abs_humidity ! absolute humidity (mass density of water vapour)
    real(wp), intent(in) :: temperature  ! temperature
    real(wp)             :: rel_humidity ! relative humidity
    
    ! local variables
    real(wp) :: vapour_pressure     ! actual water vapour pressure
    real(wp) :: saturation_pressure ! saturation water vapour pressure
    
    ! calculation of the water vapour pressure according to the equation of state
    vapour_pressure = abs_humidity*r_v*temperature
    
    if (temperature>t_0) then
      saturation_pressure = saturation_pressure_over_water(temperature)
    endif
    if (temperature<=t_0) then
      saturation_pressure = saturation_pressure_over_ice(temperature)
    endif
    
    rel_humidity = vapour_pressure/saturation_pressure
    
  end function rel_humidity
  
  function c_v_mass_weighted_air(rho,temperature,ji,jk,jl)
    
    ! This function calculates the mass-weighted c_v of the air.
    
    real(wp), intent(in) :: rho(ny,nx,n_layers,n_constituents) ! mass densities
    real(wp), intent(in) :: temperature(ny,nx,n_layers)        ! temperature
    integer,  intent(in) :: ji                                 ! horizontal index
    integer,  intent(in) :: jk                                 ! horizontal index
    integer,  intent(in) :: jl                                 ! layer index
    real(wp)             :: c_v_mass_weighted_air              ! result
    
    ! local variables
    integer :: jc ! constituent index
    
    c_v_mass_weighted_air = 0._wp
    do jc=1,n_condensed_constituents
      ! It is correct to use c_p here because the compression of the condensates has almost no effect on the air pressure.
      c_v_mass_weighted_air = c_v_mass_weighted_air + rho(ji,jk,jl,jc)*c_p_cond(jc,temperature(ji,jk,jl))
    enddo
    if (lmoist) then
      ! dry air
      c_v_mass_weighted_air = c_v_mass_weighted_air &
                              + (rho(ji,jk,jl,n_condensed_constituents+1) - rho(ji,jk,jl,n_condensed_constituents+2))*c_d_v
      ! water vapour
      c_v_mass_weighted_air = c_v_mass_weighted_air + rho(ji,jk,jl,n_condensed_constituents+2)*c_v_v
    else
      ! dry air
      c_v_mass_weighted_air = rho(ji,jk,jl,n_condensed_constituents+1)*c_d_v
    endif
    
  end function c_v_mass_weighted_air

  function calc_diff_coeff(temperature,density)
    
    ! This function calculates the molecular diffusion coefficient.
    
    real(wp) :: temperature     ! temperature
    real(wp) :: density         ! mass density
    real(wp) :: calc_diff_coeff ! result
    
    ! local variables
    real(wp) :: particle_radius  ! radius of the particles
    real(wp) :: particle_mass    ! mass of the particles
    real(wp) :: thermal_velocity ! thermal velocity of the particles
    real(wp) :: particle_density ! numeric particle density
    real(wp) :: cross_section    ! cross section of the particles
    real(wp) :: mean_free_path   ! mean free path of the particles

    ! these things are hardly ever modified
    particle_radius = 130e-12_wp
    particle_mass = m_d/n_a
    
    ! actual calculation
    thermal_velocity = sqrt(8._wp*k_b*temperature/(M_PI*particle_mass))
    particle_density = density/particle_mass
    cross_section = 4._wp*M_PI*particle_radius**2
    mean_free_path = 1._wp/(sqrt(2._wp)*particle_density*cross_section)
    calc_diff_coeff = 1._wp/3._wp*thermal_velocity*mean_free_path
    
  end function calc_diff_coeff
  
  function v_fall_liquid(state,diag,grid,radius,ji,jk,jl)
    
    ! This function returns the fall velocity of water droplets as a function of the radius of the droplets and the air density.
    
    type(t_state), intent(in)    :: state         ! state variables
    type(t_diag),  intent(inout) :: diag          ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid          ! grid quantities
    real(wp),      intent(in)    :: radius        ! radius of the droplet
    integer,       intent(in)    :: ji            ! horizontal index
    integer,       intent(in)    :: jk            ! horizontal index
    integer,       intent(in)    :: jl            ! layer index
    real(wp)                     :: v_fall_liquid ! result
    
    ! local variables
    real(wp) :: kinematic_viscosity ! kinematic viscosity
    real(wp) :: reynolds_number     ! Reynolds number of the flow around the particle
    real(wp) :: reynolds_crit       ! critical Reynolds number, where turbulent flow takes over
    real(wp) :: gravity_local       ! local gravity acceleration
    real(wp) :: c_w                 ! c_w value of a liquid particle in turbulent flow
    
    reynolds_crit = 10._wp
    
    ! computing the local gravity acceleration
    if (jl==1) then
      gravity_local = grid%gravity_m_v(ji,jk,jl+1)
    elseif (jl==n_layers) then
      gravity_local = grid%gravity_m_v(ji,jk,jl)
    else
      gravity_local &
      = grid%inner_product_weights(5,ji,jk,jl)*grid%gravity_m_v(ji,jk,jl) &
      + grid%inner_product_weights(6,ji,jk,jl)*grid%gravity_m_v(ji,jk,jl+1)
    endif
    
    kinematic_viscosity = calc_diff_coeff(diag%temperature(ji,jk,jl),state%rho(ji,jk,jl,n_condensed_constituents+1))
    
    v_fall_liquid = 2._wp*radius**2*(rho_h2o-state%rho(ji,jk,jl,n_condensed_constituents+1))*gravity_local &
                    /(9._wp*state%rho(ji,jk,jl,n_condensed_constituents+1)*kinematic_viscosity)
    
    reynolds_number = v_fall_liquid*radius/kinematic_viscosity
    
    if (reynolds_number>reynolds_crit) then
      c_w = .8_wp
      v_fall_liquid = sqrt(8._wp*radius*rho_h2o*gravity_local/(3._wp*state%rho(ji,jk,jl,n_condensed_constituents+1)*c_w))
    endif
    
  end function v_fall_liquid
  
  function v_fall_solid(state,diag,radius,ji,jk,jl)
    
    ! This function returns the fall velocity of ice particles (including snow) as a function of the radius of the droplets and the air density.
    
    type(t_state), intent(in)    :: state        ! state variables
    type(t_diag),  intent(inout) :: diag         ! diagnostic quantities
    real(wp),      intent(in)    :: radius       ! radius of the droplet
    integer,       intent(in)    :: ji           ! horizontal index
    integer,       intent(in)    :: jk           ! horizontal index
    integer,       intent(in)    :: jl           ! layer index
    real(wp)                     :: v_fall_solid ! result
    
    v_fall_solid = 0.001_wp
    
    if (radius>250.e-6) then
      v_fall_solid = 1._wp
    endif
    
  end function v_fall_solid

end module mo_derived











