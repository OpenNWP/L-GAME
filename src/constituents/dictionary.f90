! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module dictionary

  ! In this module, properties of the constituents are being stored.
  
  ! gaseous constituents IDs:
  ! 0: dry air
  ! 1: H2O

  use definitions, only: wp
  
  implicit none
  
  private
  
  public :: spec_heat_capacities_v_gas
  public :: spec_heat_capacities_p_gas
  public :: specific_gas_constants
  public :: mean_particle_masses_gas
  public :: phase_trans_heat
  public :: calc_o3_vmr
  public :: molar_fraction_in_dry_air
  
  contains

  function mean_particle_masses_gas(gas_constituent_id)
  
    ! mean particle masses of the constituents of the gas phase
    
    integer, intent(in) :: gas_constituent_id
    
    ! output
    real(wp)            :: mean_particle_masses_gas
  
    if (gas_constituent_id == 0) then
      mean_particle_masses_gas = 0.004810e-23;
    endif
  
    if (gas_constituent_id == 1) then
      mean_particle_masses_gas = 0.002991e-23
    endif
  
  end function mean_particle_masses_gas

  function spec_heat_capacities_v_gas(j_constituent)
    
    ! specific heat capacity at constant volume
    
    integer, intent(in) :: j_constituent
    
    ! output
    real(wp) :: spec_heat_capacities_v_gas
    
    if (j_constituent == 0) then
      spec_heat_capacities_v_gas = 717.942189_wp
    endif
    
    if (j_constituent == 1) then
      spec_heat_capacities_v_gas = 1396.475121_wp
    endif
    
  end function spec_heat_capacities_v_gas

  function spec_heat_capacities_p_gas(j_constituent)
    
    ! specific heat capacity at constant pressure
    
    integer, intent(in) :: j_constituent
    
    real(wp) :: spec_heat_capacities_p_gas
    
    if (j_constituent == 0) then
      spec_heat_capacities_p_gas = 1005.0_wp
    endif
    
    if (j_constituent == 1) then
      spec_heat_capacities_p_gas = 1858.0_wp
    endif
    
  end function spec_heat_capacities_p_gas
  
  function specific_gas_constants(j_constituent)
    
    ! specific gas constants
    
    integer, intent(in) :: j_constituent
    
    real(wp) :: specific_gas_constants
    
    if (j_constituent == 0) then
      specific_gas_constants = 287.057811_wp
    endif
    
    if (j_constituent == 1) then
      specific_gas_constants = 461.524879_wp
    endif
    
  end function specific_gas_constants

  function phase_trans_heat(direction,temperature)
    
    ! This function calculates the phase transition heat.

    ! input arguments
    integer  :: direction
    real(wp) :: temperature

    real(wp) :: phase_trans_heat
    
    ! directions:
    ! 0:  gas to liquid
    ! 1:  gas to solid
    ! 2: liquid to solid

    if (direction == 0) then
      phase_trans_heat = 2257000._wp
    endif
    if (direction == 1) then
      phase_trans_heat = 2257000._wp + 333500._wp
    endif
    if (direction == 2) then
      phase_trans_heat = 333500._wp
    endif
    
  end function phase_trans_heat
  
  real(wp) function calc_o3_vmr(height)

    ! This function calculates the ozone VMR as a function of height.
    ! assumes a Gaussian distribution
    
    real(wp)                          :: height            ! height above MSL
    
    ! local variables
    real(wp)                          :: fwhm = 20e3_wp      ! full width at half maximum
    real(wp)                          :: max = 34e3_wp     ! height of the maximum of the distribution
    real(wp)                          :: max_vmr = 8.5e-6_wp ! maximum volume mixing ratio
    real(wp)                          :: sigma               ! standard deviation
    real(wp)                          :: distance            ! distance from the maximum
    
    ! calculation of the result
    sigma = fwhm/(8*log(2._wp))**0.5_wp
    distance = height - max
    calc_o3_vmr = max_vmr*exp(-distance**2/(2*sigma**2))
    
  end function calc_o3_vmr

  real(wp) function molar_fraction_in_dry_air(gas_number)
    
    ! This function returns the molar fraction of certain gases in dry air.
    
    integer                           :: gas_number
  
    ! 0: dry air
    ! 1: H2O
    ! 2: N2
    ! 3: O2
    ! 4: Ar
    ! 5: CO2
    ! 6: Ne
    ! 7: He
    ! 8: CH4
    ! 9: CO
    ! 10: O3
    ! 11: N2O
  
    if (gas_number==2) then
      molar_fraction_in_dry_air = 0.7809_wp
    endif
    if (gas_number==3) then
      molar_fraction_in_dry_air = 0.2095_wp
    endif
    if (gas_number==4) then
      molar_fraction_in_dry_air = 0.0093_wp
    endif
    if (gas_number==5) then
      molar_fraction_in_dry_air = 0.0003_wp
    endif
    if (gas_number==6) then
      molar_fraction_in_dry_air=1.8e-5_wp
    endif
    if (gas_number==7) then
      molar_fraction_in_dry_air = 5.2e-6_wp
    endif
    if (gas_number==8) then
      molar_fraction_in_dry_air=1.5e-6_wp
    endif
    if (gas_number==9) then
      molar_fraction_in_dry_air=1.0e-7_wp
    endif
    if (gas_number==10) then
      molar_fraction_in_dry_air=1e-6_wp
    endif
    if (gas_number==11) then
	  ! https://www.epa.gov/climate-indicators/climate-change-indicators-atmospheric-concentrations-greenhouse-gases
      molar_fraction_in_dry_air = 0.3e-6_wp
    endif
  
  end function molar_fraction_in_dry_air

end module dictionary











