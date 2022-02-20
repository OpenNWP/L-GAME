! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/L-GAME

module planetary_boundary_layer

  ! This module computes everything related to the planetary boundary layer.
  
  use run_nml, only: EPSILON_SECURITY
  
  implicit none
  
  private
  
  public :: roughness_length_from_u10_sea
  public :: scalar_flux_resistance
  public :: momentum_flux_resistance
  public :: roughness_velocity
  public :: psi_h
  public :: psi_m
  
  contains
  
  function roughness_length_from_u10_sea(u10)
  
    ! This function returns the roughness length as a function of the mean wind speed at 10 m above a fully developed sea.

    ! refer to Stensrud, Parameterization schemes (2007), p.130

    ! empirically determined formula for the SWH
    swh = 0.0248_wp*u10**2._wp

    ! empirically determined period of the waves
    period = 0.729_wp*u10

    ! deep-water gravity waves
    wavelength = G_MEAN_SFC_ABS*period**2._wp/(2._wp*4._wp*atan(1.d0))

    ! final result
    roughness_length = 1200._wp*swh*swh/fmax(wavelength, EPSILON_SECURITY)**4.5_wp

    ! avoid too small values for stability
    return max(0.0001_wp, roughness_length)
  
  end function roughness_length_from_u10_sea

  function scalar_flux_resistance(roughness_velocity_value, z_agl, roughness_length_value, monin_obukhov_length_value)

    ! This function returns the surface flux resistance for scalar quantities.

    ! height of the prandtl layer
    used_vertical_height = fmin(z_agl, PRANDTL_HEIGHT)

    scalar_flux_resistance = 1._wp/(KARMAN*roughness_velocity_value)*
    ! neutral conditions
    (log(used_vertical_height/roughness_length_value)
    ! non-neutral conditions
    - psi_h(used_vertical_height, monin_obukhov_length_value)
    ! interfacial sublayer
    + log(7._wp))

    ! limitting the result for security
    if (scalar_flux_resistance < 1._wp) then
      scalar_flux_resistance = 1._wp
    endif 
    
  end function 

  function momentum_flux_resistance(wind_h_lowest_layer, z_agl, roughness_length_value, monin_obukhov_length_value)

    ! This function returns the surface flux resistance for momentum.

    ! height of the prandtl layer
    used_vertical_height = fmin(z_agl, PRANDTL_HEIGHT)

    momentum_flux_resistance = 1._wp/(KARMAN*roughness_velocity(wind_h_lowest_layer, z_agl, roughness_length_value))*
    ! neutral conditions
    (log(used_vertical_height/roughness_length_value)
    ! non-neutral conditions
    - psi_m(used_vertical_height, monin_obukhov_length_value))

    ! limitting the result for security
    if (momentum_flux_resistance < 1._wp) then
      momentum_flux_resistance = 1._wp
    endif

  end function momentum_flux_resistance

  function roughness_velocity(wind_speed, z_agl, roughness_length_value)

    ! This function returns the roughness velocity.

    denominator = log(z_agl/roughness_length_value)

    ! security
    if (abs(denominator) < EPSILON_SECURITY) then
      denominator = EPSILON_SECURITY
    endif

    roughness_velocity = wind_speed*KARMAN/denominator

    roughness_velocity = max(EPSILON_SECURITY, roughness_velocity)

  end function roughness_velocity

  function psi_h(z_eff, l)

    !This is a helper function for the correction to the surface scalar flux resistance for non-neutral conditions.


    ! z_eff: effective height above the surface
    ! l: Monin-Obukhov length

    ! avoiding l == 0
    if (abs(l) < EPSILON_SECURITY) then
      l = EPSILON_SECURITY
    endif

    result
    ! unstable conditions
    if (l < 0._wp) then
    ! helper variable
    x = (1._wp - 15._wp*z_eff/l)**0.25_wp

     = 2._wp*log((1._wp + x**2._wp)/2._wp)
     
    ! neutral and stable conditions
    else
      psi_h = -4._wp*z_eff/l
    endif
    
  end function psi_h

  function psi_m(z_eff, l)

    ! This is a helper function for the correction to the surface momentum flux resistance for non-neutral conditions.

    ! z_eff: effective height above the surface
    ! l: Monin-Obukhov length

    ! avoiding l == 0
    if (abs(l) < EPSILON_SECURITY) then
      l = EPSILON_SECURITY
    endif

    ! unstable conditions
    if (l < 0._wp) then
      ! helper variable
      x = (1._wp - 15._wp*z_eff/l)**0.25_wp

      psi_m = 2.0_wp*log((1._wp + x)/2._wp) + log((1._wp + x**2._wp)/2._wp) - 2._wp*atan(x) + M_PI/2._wp
    ! neutral and stable conditions
    else
      psi_m = -4.7_wp*z_eff/l
    endif
    
 end function psi_m

end module planetary_boundary_layer





