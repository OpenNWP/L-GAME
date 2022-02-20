! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module thermodynamics

  ! In this file, thermodynamic relations are calculated.
  
  use definitions, only: wp
  
  implicit none
  
  private
  
  public :: spec_heat_cap_diagnostics_v
  public :: spec_heat_cap_diagnostics_p
  public :: gas_constant_diagnostics
  
  contains

  function spec_heat_cap_diagnostics_v(j_constituent)
    
    integer, intent(in) :: j_constituent
    
    ! specific heat capacity at constant volume
    real(wp) :: spec_heat_cap_diagnostics_v
    
    if (j_constituent == 0) then
      spec_heat_cap_diagnostics_v = 717.942189_wp
    endif
    
    if (j_constituent == 1) then
      spec_heat_cap_diagnostics_v = 1396.475121_wp
    endif
    
  end function spec_heat_cap_diagnostics_v

  function spec_heat_cap_diagnostics_p(j_constituent)
    
    integer, intent(in) :: j_constituent
    
    ! specific heat capacity at constant pressure
    real(wp) :: spec_heat_cap_diagnostics_p
    
    if (j_constituent == 0) then
      spec_heat_cap_diagnostics_p = 1005.0_wp
    endif
    
    if (j_constituent == 1) then
      spec_heat_cap_diagnostics_p = 1858.0_wp
    endif
    
  end function spec_heat_cap_diagnostics_p
  
  function gas_constant_diagnostics(j_constituent)
    
    integer, intent(in) :: j_constituent
    
    ! specific heat capacity at constant pressure
    real(wp) :: gas_constant_diagnostics
    
    if (j_constituent == 0) then
      gas_constant_diagnostics = 287.057811_wp
    endif
    
    if (j_constituent == 1) then
      gas_constant_diagnostics = 461.524879_wp
    endif
    
  end function gas_constant_diagnostics

end module thermodynamics











