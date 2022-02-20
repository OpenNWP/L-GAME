! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/L-GAME

module planetary_boundary_layer

  ! This module computes everything related to the planetary boundary layer.
  
  use definitions, only: wp
  use run_nml,     only: EPSILON_SECURITY,PRANDTL_HEIGHT,gravity
  
  implicit none
  
  private
  
  real(wp) :: KARMAN = 0.4_wp
  
  public :: momentum_flux_resistance
  
  contains
  
  function roughness_length_from_u10_sea(u10)
  
    ! This function returns the roughness length as a function of the mean wind speed at 10 m above a fully developed sea.

    ! input variable
    real(wp), intent(in) :: u10
    ! output variable
    real(wp)             :: roughness_length_from_u10_sea

    ! local variables
    real(wp)             :: swh, period, wavelength

    ! refer to Stensrud, Parameterization schemes (2007), p.130

    ! empirically determined formula for the SWH
    swh = 0.0248_wp*u10**2._wp

    ! empirically determined period of the waves
    period = 0.729_wp*u10

    ! deep-water gravity waves
    wavelength = gravity*period**2._wp/(2._wp*4._wp*atan(1.d0))

    ! final result
    roughness_length_from_u10_sea = 1200._wp*swh*swh/max(wavelength,EPSILON_SECURITY)**4.5_wp

    ! avoid too small values for stability
    roughness_length_from_u10_sea = max(0.0001_wp, roughness_length_from_u10_sea)
  
  end function roughness_length_from_u10_sea

  function scalar_flux_resistance(roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value)

    ! This function returns the surface flux resistance for scalar quantities.

    ! input variable
    real(wp), intent(in) :: roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value
    ! output variable
    real(wp)             :: scalar_flux_resistance

    ! local variables
    real(wp)             :: used_vertical_height

    ! height of the prandtl layer
    used_vertical_height = min(z_agl, PRANDTL_HEIGHT)

    scalar_flux_resistance = 1._wp/(KARMAN*roughness_velocity_value) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_h(used_vertical_height, monin_obukhov_length_value) &
    ! interfacial sublayer
    + log(7._wp))

    ! limitting the result for security
    if (scalar_flux_resistance < 1._wp) then
      scalar_flux_resistance = 1._wp
    endif 
    
  end function 

  function momentum_flux_resistance(wind_h_lowest_layer,z_agl,roughness_length_value,monin_obukhov_length_value)

    ! This function returns the surface flux resistance for momentum.

    ! input variable
    real(wp), intent(in) :: wind_h_lowest_layer,z_agl,roughness_length_value,monin_obukhov_length_value
    ! output variable
    real(wp)             :: momentum_flux_resistance

    ! local variables
    real(wp)             :: used_vertical_height

    ! height of the prandtl layer
    used_vertical_height = min(z_agl, PRANDTL_HEIGHT)

    momentum_flux_resistance = 1._wp/(KARMAN*roughness_velocity(wind_h_lowest_layer, z_agl, roughness_length_value)) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_m(used_vertical_height, monin_obukhov_length_value))

    ! limitting the result for security
    if (momentum_flux_resistance < 1._wp) then
      momentum_flux_resistance = 1._wp
    endif

  end function momentum_flux_resistance

  function roughness_velocity(wind_speed,z_agl,roughness_length_value)

    ! This function returns the roughness velocity.

    ! input variable
    real(wp), intent(in) :: wind_speed,z_agl,roughness_length_value
    ! output variable
    real(wp)             :: roughness_velocity

    ! local variables
    real(wp)             :: denominator

    denominator = log(z_agl/roughness_length_value)

    ! security
    if (abs(denominator) < EPSILON_SECURITY) then
      denominator = EPSILON_SECURITY
    endif

    roughness_velocity = wind_speed*KARMAN/denominator

    roughness_velocity = max(EPSILON_SECURITY, roughness_velocity)

  end function roughness_velocity

  function psi_h(z_eff,l)

    !This is a helper function for the correction to the surface scalar flux resistance for non-neutral conditions.

    ! input variable
    real(wp), intent(in) :: z_eff,l
    ! output variable
    real(wp)             :: psi_h

    ! local variables
    real(wp)             :: x,l_local

    ! z_eff: effective height above the surface
    ! l: Monin-Obukhov length

    ! avoiding l == 0
    l_local = l
    if (abs(l_local) < EPSILON_SECURITY) then
      l_local = EPSILON_SECURITY
    endif

    ! unstable conditions
    if (l_local < 0._wp) then
      ! helper variable
      x = (1._wp - 15._wp*z_eff/l_local)**0.25_wp
      psi_h = 2._wp*log((1._wp + x**2._wp)/2._wp)     
    ! neutral and stable conditions
    else
      psi_h = -4._wp*z_eff/l_local
    endif
    
  end function psi_h

  function psi_m(z_eff,l)

    ! This is a helper function for the correction to the surface momentum flux resistance for non-neutral conditions.

    ! input variable
    real(wp), intent(in) :: z_eff,l
    ! output variable
    real(wp)             :: psi_m

    ! local variables
    real(wp)             :: x,l_local

    ! z_eff: effective height above the surface
    ! l: Monin-Obukhov length

    ! avoiding l == 0
    l_local = l
    if (abs(l_local) < EPSILON_SECURITY) then
      l_local = EPSILON_SECURITY
    endif

    ! unstable conditions
    if (l_local < 0._wp) then
      ! helper variable
      x = (1._wp - 15._wp*z_eff/l_local)**0.25_wp

      psi_m = 2.0_wp*log((1._wp + x)/2._wp) + log((1._wp + x**2._wp)/2._wp) - 2._wp*atan(x) + 4._wp*atan(1.d0)/2._wp
    ! neutral and stable conditions
    else
      psi_m = -4.7_wp*z_eff/l_local
    endif
    
 end function psi_m

end module planetary_boundary_layer





