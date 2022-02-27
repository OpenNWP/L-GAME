! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module humidity

  ! In this module, humidity quantities are being calculated.
  
  use definitions, only: wp
  use constants,   only: T_0
  use dictionary,  only: specific_gas_constants
  
  implicit none
  
  private
  
  public :: saturation_pressure_over_water
  public :: saturation_pressure_over_ice
  public :: rel_humidity
  
  contains
	
  function saturation_pressure_over_water(temperature)

    ! This function returns the saturation pressure in Pa over liquid water as a function of the temperature in K.
    
    ! input argument
    real(wp), intent(in) :: temperature                    ! gas temperature
    ! output argument
    real(wp)             :: saturation_pressure_over_water ! the result
    
    ! local variables
    real(wp)  :: temp_c ! temperature in degrees Celsius
    
    temp_c = temperature - T_0
    saturation_pressure_over_water = 100.0_wp*6.112_wp*exp(17.62_wp*temp_c/(243.12_wp+temp_c))

  end function saturation_pressure_over_water

  function saturation_pressure_over_ice(temperature)
    
    ! This function returns the saturation pressure in Pa over ice as a function of the temperature in K.
    
    ! input argument
    real(wp), intent(in) :: temperature
    ! output argument
    real(wp)             :: saturation_pressure_over_ice
    
    ! local variables
    real(wp)             :: temp_c

    temp_c = temperature - T_0
    saturation_pressure_over_ice = 100.0_wp*6.112_wp*exp(22.46_wp*temp_c/(272.62_wp+temp_c))

  end function saturation_pressure_over_ice

  function rel_humidity(abs_humidity,temperature)
    
    ! This function returns the relative humidity as a function of the absolute humidity in kg/m^3 and the temperature in K.
    
    ! input arguments
    real(wp), intent(in) :: abs_humidity,temperature
    ! output argument
    real(wp)             :: rel_humidity
    
    ! local variables
    real(wp)             :: vapour_pressure     ! actual water vapour pressure
    real(wp)             :: saturation_pressure ! saturation water vapour pressure
    
    vapour_pressure = abs_humidity*specific_gas_constants(1)*temperature
    
    if (temperature>T_0) then
      saturation_pressure = saturation_pressure_over_water(temperature)
    endif
    if (temperature<=T_0) then
      saturation_pressure = saturation_pressure_over_ice(temperature)
    endif
    
    rel_humidity = vapour_pressure/saturation_pressure
    
  end function rel_humidity

end module humidity











