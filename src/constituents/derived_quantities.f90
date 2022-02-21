! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module derived_quantities

  ! In this module, more complex thermodynamic quantities are being calculated.
  
  use definitions, only: wp
  
  implicit none
  
  private
  
  public :: temperature_diagnostics
  public :: spec_heat_cap_diagnostics_v
  public :: spec_heat_cap_diagnostics_p
  public :: gas_constant_diagnostics
  public :: density_total
  public :: density_gas
  public :: calc_diffusion_coeff
  
  contains

  subroutine temperature_diagnostics(state,grid,diag)
    
    ! This function diagnoses the temperature of the gas phase.

    #pragma omp parallel for
    do (i = 0 i < NO_OF_SCALARS ++i)
      diagnostics -> temperature_gas[i] = (grid -> theta_bg[i] + state -> theta_pert[i])*(grid -> exner_bg[i] + state -> exner_pert[i])
    enddo
    
  end subroutine temperature_diagnostics

  function spec_heat_cap_diagnostics_v(state,grid_point_index,config)
  
    ! output argument
    real(wp)              :: spec_heat_cap_diagnostics_v
    
    rho_g = 0._wp
    no_of_relevant_constituents = 0
    if (.not. config%assume_lte) then
      no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS
      rho_g = density_gas(state,grid_point_index)
    endif
    
    if (config%lassume_lte) then
      no_of_relevant_constituents = 1
      rho_g = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + grid_point_index]
    endif
    
    spec_heat_cap_diagnostics_v = 0._wp
    do i=1,no_of_relevant_constituents
      spec_heat_cap_diagnostics_v += state -> rho[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*spec_heat_capacities_v_gas(i)
    enddo
    
  end function spec_heat_cap_diagnostics_v

  function spec_heat_cap_diagnostics_p(state,grid_point_index,config)
  
    ! output argument
    real(wp)              :: spec_heat_cap_diagnostics_p
    
    rho_g = 0._wp
    no_of_relevant_constituents = 0
    if (.not. config%assume_lte) then
      no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS
      rho_g = density_gas(state,grid_point_index)
    endif
    
    if (config%assume_lte) then
      no_of_relevant_constituents = 1
      rho_g = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + grid_point_index]
    endif
    
    spec_heat_cap_diagnostics_p = 0._wp
    do (i = 0 i < no_of_relevant_constituents ++i)
      spec_heat_cap_diagnostics_p += state -> rho[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*spec_heat_capacities_p_gas(i)
    enddo
    
  end function spec_heat_cap_diagnostics_p
    
  function gas_constant_diagnostics(state,grid_point_index,config)
  
    ! output argument
    real(wp)              :: gas_constant_diagnostics
    
    rho_g = 0._wp
    no_of_relevant_constituents = 0
    
    if (config%assume_lte == 0) then
      no_of_relevant_constituents = NO_OF_GASEOUS_CONSTITUENTS
      rho_g = density_gas(state,grid_point_index)
    endif
    
    if (config%assume_lte == 1)
      no_of_relevant_constituents = 1
      rho_g = state -> rho[NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + grid_point_index]
    endif
    
    gas_constant_diagnostics = 0
    
    do i=1,no_of_relevant_constituents
      gas_constant_diagnostics += state -> rho[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]/rho_g*specific_gas_constants(i)
    enddo
    
  end function gas_constant_diagnostics

  function density_total(state,grid_point_index)
  
    ! output argument
    real(wp)              :: density_total
    
    density_total = 0._wp
    
    do i=1,NO_OF_CONSTITUENTS
      density_total += state -> rho[i*NO_OF_SCALARS + grid_point_index]
    enddo
    
  end functiondensity_total

  function density_gas(state,grid_point_index)
    ! output argument
    real(wp)              :: density_total
    
    density_gas = 0._wp
    
    do i=1,no_of_gaseous_constituents
      density_gas = density_gas + state -> rho[(i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + grid_point_index]
    enddo
    
  end function density_gas

  function calc_diffusion_coeff(temperature,density)
  
    ! output argument
    real(wp)              :: calc_diffusion_coeff
    
    ! local variables
    real(wp)              :: particle_radius,particle_mass,thermal_velocity,particle_density,cross_section,mean_free_path
    
    ! This function calculates the molecular diffusion coefficient according to the kinetic gas theory.

    ! these things are hardly ever modified
    particle_radius = 130e-12
    particle_mass = mean_particle_masses_gas(0)
    
    ! actual calculation
    thermal_velocity = sqrt(8.0*K_B*temperature/(M_PI*particle_mass))
    particle_density = density/particle_mass
    cross_section = 4.0*M_PI*pow(particle_radius, 2.0)
    mean_free_path = 1.0/(sqrt(2.0)*particle_density*cross_section)
    calc_diffusion_coeff = 1.0/3.0*thermal_velocity*mean_free_path
    
  end function calc_diffusion_coeff

end module derived_quantities











