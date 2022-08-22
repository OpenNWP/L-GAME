! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_derived

  ! In this module more complex thermodynamic quantities are being calculated.
  
  use mo_definitions,   only: wp,t_grid,t_state,t_diag
  use run_nml,          only: ny,nx,nlays
  use constants,        only: k_b,M_PI,m_d,n_a,r_d,r_v,c_d_p,c_v_p,c_d_v,c_v_v,t_0
  use constituents_nml, only: n_condensed_constituents,n_gaseous_constituents,n_constituents
  use dictionary,       only: saturation_pressure_over_ice,saturation_pressure_over_water
  
  implicit none
  
  contains

  subroutine temperature_diagnostics(state,diag,grid)
    
    ! This function diagnoses the temperature of the gas phase.
    
    ! input arguments
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx
        do jl=1,nlays
          diag%temperature(ji,jk,jl) = (grid%theta_v_bg(ji,jk,jl) + state%theta_v_pert(ji,jk,jl)) &
          *(grid%exner_bg(ji,jk,jl) + state%exner_pert(ji,jk,jl))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine temperature_diagnostics

  function spec_heat_cap_diagnostics_v(state,ji,jk,jl)
  
    ! This function calculates the specific heat capacity at constant volume of the air at a certain gridpoint.
  
    ! input arguments
    type(t_state), intent(in) :: state    ! state with which to calculate c_v
    integer, intent(in)       :: ji,jk,jl ! spatial indices of the gridpoint
    
    ! output
    real(wp) :: spec_heat_cap_diagnostics_v
    
    ! local variables
    integer  :: n_relevant_constituents ! the number of relevant constituents for this calculation
    integer  :: jc               ! constituent index
    real(wp) :: rho_g                       ! gas density
    real(wp) :: spec_heat_capacities_v_gas(2)
    
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
  
    ! input arguments
    type(t_state), intent(in) :: state    ! state with which to calculate c_p
    integer,       intent(in) :: ji,jk,jl ! spatial indices of the gridpoint
    
    ! output
    real(wp) :: spec_heat_cap_diagnostics_p
    
    ! local variables
    integer  :: n_relevant_constituents ! the number of relevant constituents for this calculation
    integer  :: jc               ! constituent index
    real(wp) :: rho_g                       ! gas density
    real(wp) :: spec_heat_capacities_p_gas(2)
    
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
  
    ! input arguments
    type(t_state), intent(in) :: state    ! state with which to calculate the gas constant
    integer,       intent(in) :: ji,jk,jl ! spatial indices of the gridpoint
    ! output
    real(wp)                  :: gas_constant_diagnostics
    
    ! local variables
    integer  :: n_relevant_constituents ! the number of relevant constituents for this calculation
    integer  :: jc               ! constituent index
    real(wp) :: rho_g                       ! gas density
    real(wp) :: specific_gas_constants(2)
    
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

  function density_gas(state,ji,jk,jl)
    
    ! This function calculates the density of the gas phase at a certain gridpoint.
    
    ! input arguments
    type(t_state), intent(in) :: state       ! state with which to calculate the result
    integer,       intent(in) :: ji,jk,jl    ! indices of the gridpoint
    ! output
    real(wp)                  :: density_gas ! the result
    
    ! local variables
    integer :: jc ! index of the constituent
    
    density_gas = 0._wp
    
    do jc=1,n_gaseous_constituents
      density_gas = density_gas + state%rho(ji,jk,jl,n_condensed_constituents+jc)
    enddo
    
  end function density_gas

  function rel_humidity(abs_humidity,temperature)
    
    ! This function returns the relative humidity as a function of the absolute humidity in kg/m^3 and the temperature in K.
    
    real(wp), intent(in) :: abs_humidity,temperature
    real(wp)             :: rel_humidity
    
    ! local variables
    real(wp)             :: vapour_pressure     ! actual water vapour pressure
    real(wp)             :: saturation_pressure ! saturation water vapour pressure
    
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

  function calc_diffusion_coeff(temperature,density)
  
    ! This function calculates the molecular diffusion coefficient.
  
    real(wp) :: temperature,density
    real(wp) :: calc_diffusion_coeff
    
    ! local variables
    real(wp) :: particle_radius,particle_mass,thermal_velocity,particle_density,cross_section,mean_free_path

    ! these things are hardly ever modified
    particle_radius = 130e-12_wp
    particle_mass = m_d/n_a
    
    ! actual calculation
    thermal_velocity = sqrt(8.0_wp*k_b*temperature/(M_PI*particle_mass))
    particle_density = density/particle_mass
    cross_section = 4.0_wp*M_PI*particle_radius**2
    mean_free_path = 1.0_wp/(sqrt(2.0_wp)*particle_density*cross_section)
    calc_diffusion_coeff = 1.0_wp/3.0_wp*thermal_velocity*mean_free_path
    
  end function calc_diffusion_coeff

end module mo_derived











