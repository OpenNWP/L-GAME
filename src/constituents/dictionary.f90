! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module dictionary

  ! In this file, properties of the constituents are being stored.
  
  ! gaseous constituents IDs:
  ! 0: dry air
  ! 1: H2O

  use definitions, only: wp
  
  implicit none
  
  private
  
  public :: spec_heat_capacities_v_gas
  public :: spec_heat_capacities_p_gas
  public :: specific_gas_constants
  
  contains

  function spec_heat_capacities_v_gas(j_constituent)
    
    integer, intent(in) :: j_constituent
    
    ! specific heat capacity at constant volume
    real(wp) :: spec_heat_capacities_v_gas
    
    if (j_constituent == 0) then
      spec_heat_capacities_v_gas = 717.942189_wp
    endif
    
    if (j_constituent == 1) then
      spec_heat_capacities_v_gas = 1396.475121_wp
    endif
    
  end function spec_heat_capacities_v_gas

  function spec_heat_capacities_p_gas(j_constituent)
    
    integer, intent(in) :: j_constituent
    
    ! specific heat capacity at constant pressure
    real(wp) :: spec_heat_capacities_p_gas
    
    if (j_constituent == 0) then
      spec_heat_capacities_p_gas = 1005.0_wp
    endif
    
    if (j_constituent == 1) then
      spec_heat_capacities_p_gas = 1858.0_wp
    endif
    
  end function spec_heat_capacities_p_gas
  
  function specific_gas_constants(j_constituent)
    
    integer, intent(in) :: j_constituent
    
    ! specific gas constants
    real(wp) :: specific_gas_constants
    
    if (j_constituent == 0) then
      specific_gas_constants = 287.057811_wp
    endif
    
    if (j_constituent == 1) then
      specific_gas_constants = 461.524879_wp
    endif
    
  end function specific_gas_constants

end module dictionary











