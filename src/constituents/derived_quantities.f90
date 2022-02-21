! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module derived_quantities

  ! In this module, more complex thermodynamic quantities are being calculated.
  
  use definitions,      only: wp,t_grid,t_state,t_diag,t_config
  use dictionary,       only: mean_particle_masses_gas
  use run_nml,          only: nlins,ncols,nlays
  use constituents_nml, only: no_of_condensed_constituents,no_of_gaseous_constituents,no_of_constituents
  
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
    
    ! input arguments
    type(t_state), intent(in)    :: state
    type(t_grid),  intent(in)    :: grid
    type(t_diag),  intent(inout) :: diag
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,comp_h,comp_v)
    do ji=1,nlins+2
      do jk=1,ncols+2
        do jl=1,nlays
          diag%temperature_gas(ji,jk,hl) = (grid%theta_bg(ji,jk,jl) + state%theta_pert(ji,jk,jl)) &
          *(grid%exner_bg(ji,jk,jl) + state%exner_pert(ji,jk,jl))
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
  end subroutine temperature_diagnostics

  function spec_heat_cap_diagnostics_v(state,ji,jk,jl,config)
  
    ! input arguments
    type(t_state), intent(in)    :: state
    integer, intent(in)          :: ji,jk,jl
    type(t_config), intent(in)   :: config
    
    ! output
    real(wp)              :: spec_heat_cap_diagnostics_v
    
    rho_g = 0._wp
    no_of_relevant_constituents = 0
    if (.not. config%lassume_lte) then
      no_of_relevant_constituents = no_of_gaseous_constituents
      rho_g = density_gas(state,ji,jk,jl)
    endif
    
    if (config%lassume_lte) then
      no_of_relevant_constituents = 1
      rho_g = state%rho(ji,jk,jl,no_of_condensed_constituents+1)
    endif
    
    spec_heat_cap_diagnostics_v = 0._wp
    do j_constituent=1,no_of_relevant_constituents
      spec_heat_cap_diagnostics_v = spec_heat_cap_diagnostics_v +  state%rho(ji,jk,jl,j_constituent) &
      /rho_g*spec_heat_capacities_v_gas(j_constituent)
    enddo
    
  end function spec_heat_cap_diagnostics_v

  function spec_heat_cap_diagnostics_p(state,ji,jk,jl,config)
  
    ! input arguments
    type(t_state), intent(in)  :: state
    integer, intent(in)        :: ji,jk,jl
    type(t_config), intent(in) :: config
    
    ! output
    real(wp)                   :: spec_heat_cap_diagnostics_p
    integer                    :: no_of_relevant_constituents
    
    rho_g = 0._wp
    no_of_relevant_constituents = 0
    
    if (.not. config%lassume_lte) then
      no_of_relevant_constituents = no_of_gaseous_constituents
      rho_g = density_gas(state,ji,jk,jl)
    endif
    
    if (config%lassume_lte) then
      no_of_relevant_constituents = 1
      rho_g = state%rho(NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + ji,jk,jl)
    endif
    
    spec_heat_cap_diagnostics_p = 0._wp
    do j_constituent=1,no_of_relevant_constituents
      spec_heat_cap_diagnostics_p = spec_heat_cap_diagnostics_p + state%rho(ji,jk,jl,j_constituents) &
      /rho_g*spec_heat_capacities_p_gas(i)
    enddo
    
  end function spec_heat_cap_diagnostics_p
    
  function gas_constant_diagnostics(state,ji,jk,jl,config)
  
    ! input arguments
    type(t_state), intent(in)  :: state
    integer, intent(in)        :: ji,jk,jl
    type(t_config), intent(in) :: config
    
    ! output
    real(wp)                   :: gas_constant_diagnostics
    integer                    :: no_of_relevant_constituents
    
    rho_g = 0._wp
    no_of_relevant_constituents = 0
    
    if (.not. config%lassume_lte) then
      no_of_relevant_constituents = no_of_gaseous_constituents
      rho_g = density_gas(state,ji,jk,jl)
    endif
    
    if (config%lassume_lte) then
      no_of_relevant_constituents = 1
      rho_g = state%rho(NO_OF_CONDENSED_CONSTITUENTS*NO_OF_SCALARS + ji,jk,jl)
    endif
    
    gas_constant_diagnostics = 0._wp
    
    do j_constituent=1,no_of_relevant_constituents
      gas_constant_diagnostics = gas_constant_diagnostics + state%rho((i + NO_OF_CONDENSED_CONSTITUENTS)*NO_OF_SCALARS + ji,jk,jl) &
      /rho_g*specific_gas_constants(j_constituent)
    enddo
    
  end function gas_constant_diagnostics

  function density_total(state,ji,jk,jl)
  
    ! This function calculates the total density of the air at a certain gridpoint.
    
    ! input arguments
    type(t_state), intent(in) :: state
    integer, intent(in)       :: ji,jk,jl
    
    ! output
    real(wp)                  :: density_total
    
    density_total = 0._wp
    
    do j_constituent=1,no_of_constituents
      density_gas = density_gas + state%rho(ji,jk,jl,j_constituent)
    enddo
    
  end function density_total

  function density_gas(state,ji,jk,jl)
    
    ! input arguments
    type(t_state), intent(in)    :: state
    integer, intent(in)          :: ji,jk,jl
    
    ! This function calculates the density of the gas phase at a certain gridpoint.
  
    ! input arguments
    integer, intent(in)   :: ji,jk,jl
    
    ! output
    real(wp)              :: density_gas
    
    density_gas = 0._wp
    
    do j_constituent=1,no_of_gaseous_constituents
      density_gas = density_gas + state%rho(ji,jk,jl,no_of_condensed_constituents+j_constituent)
    enddo
    
  end function density_gas

  function calc_diffusion_coeff(temperature,density)
  
    ! input arguments
    real(wp)              :: temperature,density
  
    ! output
    real(wp)              :: calc_diffusion_coeff
    
    ! local variables
    real(wp)              :: particle_radius,particle_mass,thermal_velocity,particle_density,cross_section,mean_free_path
    
    ! This function calculates the molecular diffusion coefficient according to the kinetic gas theory.

    ! these things are hardly ever modified
    particle_radius = 130e-12_wp
    particle_mass = mean_particle_masses_gas(0)
    
    ! actual calculation
    thermal_velocity = sqrt(8.0_wp*K_B*temperature/(M_PI*particle_mass))
    particle_density = density/particle_mass
    cross_section = 4.0_wp*M_PI*pow(particle_radius, 2.0_wp)
    mean_free_path = 1.0_wp/(sqrt(2.0_wp)*particle_density*cross_section)
    calc_diffusion_coeff = 1.0_wp/3.0_wp*thermal_velocity*mean_free_path
    
  end function calc_diffusion_coeff

end module derived_quantities











