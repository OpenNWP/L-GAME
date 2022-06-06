! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module dictionary

  ! In this module, properties of the constituents are being stored.
  
  ! gaseous constituents IDs:
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

  use definitions, only: wp
  use constants,   only: T_0,N_A
  
  implicit none
  
  private
  
  public :: spec_heat_capacities_v_gas
  public :: spec_heat_capacities_p_gas
  public :: specific_gas_constants
  public :: mean_particle_masses_gas
  public :: phase_trans_heat
  public :: calc_o3_vmr
  public :: molar_fraction_in_dry_air
  public :: saturation_pressure_over_water
  public :: saturation_pressure_over_ice
  public :: enhancement_factor_over_water
  public :: enhancement_factor_over_ice
  public :: rel_humidity
  
  contains

  function mean_particle_masses_gas(gas_constituent_id)
  
    ! mean particle masses of the constituents of the gas phase
    
    integer, intent(in) :: gas_constituent_id ! number of the gas constituent
    ! the result
    real(wp)            :: mean_particle_masses_gas
  
    if (gas_constituent_id==0) then
      mean_particle_masses_gas = 0.004810e-23;
    endif
  
    if (gas_constituent_id==1) then
      mean_particle_masses_gas = 0.002991e-23
    endif
  
  end function mean_particle_masses_gas

  function spec_heat_capacities_v_gas(gas_constituent_id)
    
    ! specific heat capacity at constant volume
    
    integer, intent(in) :: gas_constituent_id ! number of the gas constituent
    ! the result
    real(wp) :: spec_heat_capacities_v_gas
    
    if (gas_constituent_id==0) then
      spec_heat_capacities_v_gas = 717.942189_wp
    endif
    
    if (gas_constituent_id==1) then
      spec_heat_capacities_v_gas = 1396.475121_wp
    endif
    
  end function spec_heat_capacities_v_gas

  function spec_heat_capacities_p_gas(gas_constituent_id)
    
    ! specific heat capacity at constant pressure
    
    integer, intent(in) :: gas_constituent_id ! number of the gas constituent
    ! the result
    real(wp) :: spec_heat_capacities_p_gas
    
    if (gas_constituent_id==0) then
      spec_heat_capacities_p_gas = 1005.0_wp
    endif
    
    if (gas_constituent_id==1) then
      spec_heat_capacities_p_gas = 1858.0_wp
    endif
    
  end function spec_heat_capacities_p_gas
  
  function specific_gas_constants(gas_constituent_id)
    
    ! specific gas constants
    
    integer, intent(in) :: gas_constituent_id ! number of the gas constituent
    ! the result
    real(wp) :: specific_gas_constants
    
    if (gas_constituent_id==0) then
      specific_gas_constants = 287.057811_wp
    endif
    
    if (gas_constituent_id==1) then
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
    ! 0: gas to liquid
    ! 1: gas to solid
    ! 2: liquid to solid

    phase_trans_heat = 0._wp
    if (direction==0) then
      phase_trans_heat = enthalpy_evaporation(temperature)
    endif
    if (direction==1) then
      phase_trans_heat = enthalpy_sublimation(temperature)
    endif
    if (direction==2) then
      phase_trans_heat = enthalpy_sublimation(temperature) - enthalpy_evaporation(temperature)
    endif
    
  end function phase_trans_heat
  
  function enthalpy_evaporation(temperature)
    
    ! This function returns the enthalpy of evaporation depending on the temperature.
    
    real(wp) :: temperature

    real(wp) :: enthalpy_evaporation
    ! local variable
    real(wp) :: temp_c
    
    ! temperature in degrees Celsius
    temp_c = temperature - T_0
    
    if (temperature<T_0) then
      ! This follows Eq. (9) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
      ! QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.
      ! clipping values that are too extreme for these approximations
      if (temperature<30._wp) then
        temperature = 30._wp
      endif
      enthalpy_evaporation = 56579._wp - 42.212_wp*temperature + exp(0.1149_wp*(281.6_wp - temperature))
      ! unit conversion from J/mol to J/kg
      enthalpy_evaporation = enthalpy_evaporation/(N_A*mean_particle_masses_gas(1))
    else
      ! This follows the formula (Eq. (8)) cited by Huang:
      ! A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
      ! clipping values that are too extreme for these approximations
      if (temperature>T_0+100._wp) then
        temperature = T_0 + 100._wp
      endif
      enthalpy_evaporation = 3151378._wp - 2386._wp*temperature
    endif

  end function enthalpy_evaporation
  
  function enthalpy_sublimation(temperature)
    
	! This function returns the enthalpy of sublimation depending on the temperature.
	! It follows Eq. (5) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
	! QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.
    
    ! input
    real(wp) :: temperature
    ! output
    real(wp) :: enthalpy_sublimation

    ! clipping values that are too extreme for this approximation
    if (temperature<30.0) then
      temperature = 30.0
    endif
    ! sublimation is not happening in thermodynamic equilibrium at temperatures > T_0
    if (temperature>T_0) then
      temperature = T_0
    endif

    enthalpy_sublimation = 46782.5 + 35.8925*temperature - 0.07414*temperature**2 + 541.5*exp(-(temperature/123.75)**2)

    ! unit conversion from J/mol to J/kg
    enthalpy_sublimation = enthalpy_sublimation/(N_A*mean_particle_masses_gas(1))
  
  end function enthalpy_sublimation
  
  real(wp) function calc_o3_vmr(height)

    ! This function calculates the ozone VMR as a function of height.
    ! assumes a Gaussian distribution
    
    real(wp) :: height ! height above MSL
    
    ! local variables
    real(wp) :: fwhm = 20e3_wp      ! full width at half maximum
    real(wp) :: max = 34e3_wp       ! height of the maximum of the distribution
    real(wp) :: max_vmr = 8.5e-6_wp ! maximum volume mixing ratio
    real(wp) :: sigma               ! standard deviation
    real(wp) :: distance            ! distance from the maximum
    
    ! calculation of the result
    sigma = fwhm/(8._wp*log(2._wp))**0.5_wp
    distance = height - max
    calc_o3_vmr = max_vmr*exp(-distance**2/(2._wp*sigma**2))
    
  end function calc_o3_vmr

  real(wp) function molar_fraction_in_dry_air(gas_number)
    
    ! This function returns the molar fraction of certain gases in dry air.
    
    integer :: gas_number
  
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
      molar_fraction_in_dry_air = 1.8e-5_wp
    endif
    if (gas_number==7) then
      molar_fraction_in_dry_air = 5.2e-6_wp
    endif
    if (gas_number==8) then
      molar_fraction_in_dry_air = 1.5e-6_wp
    endif
    if (gas_number==9) then
      molar_fraction_in_dry_air = 1.0e-7_wp
    endif
    if (gas_number==10) then
      molar_fraction_in_dry_air = 1e-6_wp
    endif
    if (gas_number==11) then
	  ! https://www.epa.gov/climate-indicators/climate-change-indicators-atmospheric-concentrations-greenhouse-gases
      molar_fraction_in_dry_air = 0.3e-6_wp
    endif
  
  end function molar_fraction_in_dry_air
	
  function saturation_pressure_over_water(temperature)

    ! This function returns the saturation pressure in Pa over liquid water as a function of the temperature in K.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
    
    ! input argument
    real(wp), intent(in) :: temperature
    ! output argument
    real(wp)             :: saturation_pressure_over_water
    
    ! local variables
    real(wp)  :: temp_c ! temperature in degrees Celsius

    temp_c = temperature - T_0
    ! clipping too extreme values for this approximation
    if (temp_c>100._wp) then
      temp_c = 100._wp
    endif
    
    if (temp_c>0._wp) then
      saturation_pressure_over_water = exp(34.494_wp-4924.99_wp/(temp_c+237.1_wp))/(temp_c+105._wp)**1.57_wp
    ! For super-cooled water we use the formula cited in Pruppacher and Klett (2010), p. 854, Eq. (A.4-1).
    else
      ! Clipping values that are too extreme for this approximation.
      if (temp_c<-50.0) then
        temp_c = -50.0
      endif
      saturation_pressure_over_water &
      = 6.107799961 &
      + 4.436518521e-1*temp_c &
      + 1.428945805e-2*temp_c**2 &
      + 2.650648471e-4*temp_c**3 &
      + 3.031240396e-6*temp_c**4 &
      + 2.034080948e-8*temp_c**5 &
      + 6.136820929e-11*temp_c**6
    endif

  end function saturation_pressure_over_water

  function saturation_pressure_over_ice(temperature)
    
    ! This function returns the saturation pressure in Pa over ice as a function of the temperature in K.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
    
    ! input argument
    real(wp), intent(in) :: temperature
    ! output argument
    real(wp)             :: saturation_pressure_over_ice
    
    ! local variables
    real(wp)             :: temp_c

    temp_c = temperature - T_0

    ! clipping too extreme values for this approximation
    if (temp_c<-70._wp) then
      temp_c = -70._wp
    endif
    if (temp_c>0._wp) then
      temp_c = 0._wp
    endif
    
    saturation_pressure_over_ice = exp(43.494_wp-6545.8_wp/(temp_c+278._wp))/(temp_c+868._wp)**2
    
  end function saturation_pressure_over_ice
  
  function enhancement_factor_over_water(air_pressure)

    ! This function calculates the enhancement factor over water, which accounts for the fact the the saturation vapour pressure is different in moist air compared to pure water vapour.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.

    ! input arguments
    real(wp), intent(in) :: air_pressure
    ! output argument
    real(wp)             :: enhancement_factor_over_water
    
    enhancement_factor_over_water = 0.99882*exp(0.00000008*air_pressure)

    end function enhancement_factor_over_water
  
  function enhancement_factor_over_ice(air_pressure)

    ! This function calculates the enhancement factor over ice, which accounts for the fact the the saturation vapour pressure is different in moist air compared to pure water vapour.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.

    ! input arguments
    real(wp), intent(in) :: air_pressure
    ! output argument
    real(wp)             :: enhancement_factor_over_ice
    
    enhancement_factor_over_ice = 0.99882*exp(0.00000008*air_pressure)

    end function enhancement_factor_over_ice

  function rel_humidity(abs_humidity,temperature)
    
    ! This function returns the relative humidity as a function of the absolute humidity in kg/m^3 and the temperature in K.
    
    ! input arguments
    real(wp), intent(in) :: abs_humidity,temperature
    ! output argument
    real(wp)             :: rel_humidity
    
    ! local variables
    real(wp)             :: vapour_pressure     ! actual water vapour pressure
    real(wp)             :: saturation_pressure ! saturation water vapour pressure
    
    ! calculation of the water vapour pressure according to the equation of state
    vapour_pressure = abs_humidity*specific_gas_constants(1)*temperature
    
    if (temperature>T_0) then
      saturation_pressure = saturation_pressure_over_water(temperature)
    endif
    if (temperature<=T_0) then
      saturation_pressure = saturation_pressure_over_ice(temperature)
    endif
    
    rel_humidity = vapour_pressure/saturation_pressure
    
  end function rel_humidity

end module dictionary











