! This source file is part of the Limited-area GAME version (L-GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/L-GAME

module planetary_boundary_layer

  ! This module computes everything related to the planetary boundary layer.
  
  use definitions, only: wp,t_state,t_grid,t_diag
  use run_nml,     only: PRANDTL_HEIGHT,nlins,ncols,nlays,dtime
  use constants,   only: EPSILON_SECURITY,M_PI,gravity
  use surface_nml, only: lsoil
  
  implicit none
  
  private
  
  real(wp) :: KARMAN = 0.4_wp
  
  public :: momentum_flux_resistance
  public :: update_sfc_turb_quantities
  
  contains

  subroutine update_sfc_turb_quantities(state,diag,grid)

    ! This subroutine updates surface-related turbulence quantities.
    
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state with which to calculate the turbulence quantities
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(inout) :: grid  ! grid properties

    ! local variables
    real(wp) :: u_lowest_layer        ! wind speed in the lowest model layer
    real(wp) :: u10                   ! wind speed in 10 m height
    real(wp) :: z_agl                 ! height above ground f the lowest model layer
    real(wp) :: theta_lowest_layer    ! potential temperature in the lowest layer
    real(wp) :: theta_second_layer    ! potential temperature in the second-lowest layer
    real(wp) :: dz                    ! vertical grid point distance
    real(wp) :: dtheta_dz             ! vertical gradient of the potential temperature
    real(wp) :: w_pert                ! vertical velocity peturbation near the surface
    real(wp) :: theta_pert            ! potential temperature peturbation near the surface
    real(wp) :: w_pert_theta_pert_avg ! correlation between vertical velocity and potential temperature peturbations
    real(wp) :: prop_coeff            ! semi-empirical coefficient for computing w_pert_theta_pert_avg
    integer  :: ji,jk                 ! loop variables

    prop_coeff = 0.2_wp
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,u_lowest_layer,u10,z_agl,theta_lowest_layer,theta_second_layer,dz,dtheta_dz, &
    !$OMP w_pert,theta_pert,w_pert_theta_pert_avg)
    do ji=1,nlins
      do jk=1,ncols
        z_agl = grid%z_geo_scal(ji,jk,nlays) - grid%z_geo_w(ji,jk,nlays+1)

        ! wind speed in the lowest layer
        u_lowest_layer = diag%v_squared(ji,jk,nlays)**0.5_wp

        ! calculating the 10 m wind velocity from the logarithmic wind profile
        u10 = u_lowest_layer*log(10._wp/grid%roughness_length(ji,jk))/log(z_agl/grid%roughness_length(ji,jk))

        ! only over the sea the roughness length is time-dependant (because of the waves)
        if (grid%is_land(ji,jk)==0) then
          ! calculating the roughness length fom the wind velocity
          grid%roughness_length(ji,jk) = roughness_length_from_u10_sea(u10)
        endif

        ! updating the roughness velocity
        diag%roughness_velocity(ji,jk) = roughness_velocity(u_lowest_layer,z_agl,grid%roughness_length(ji,jk))

        ! theta in the lowest layer
        theta_lowest_layer = grid%theta_bg(ji,jk,nlays) + state%theta_pert(ji,jk,nlays)
        ! theta in the second-lowest layer
        theta_second_layer = grid%theta_bg(ji,jk,nlays-1) + state%theta_pert(ji,jk,nlays-1)

        ! delta z
        dz = grid%z_geo_scal(ji,jk,nlays-1) - grid%z_geo_scal(ji,jk,nlays)

        ! vertical gradient of theta
        dtheta_dz = (theta_second_layer - theta_lowest_layer)/dz

        ! the perturbation of the vertical velocity is assumed to be proportional to the 10 m wind speed
        ! times a stability-dependant factor
        w_pert = u10*max(0.001_wp,0.02_wp*(1._wp - dtheta_dz/0.01_wp))
        theta_pert = -0.2_wp*dtime*w_pert*dtheta_dz
        w_pert_theta_pert_avg = prop_coeff*w_pert*theta_pert

        ! security
        if (abs(w_pert_theta_pert_avg)<EPSILON_SECURITY) then
          w_pert_theta_pert_avg = EPSILON_SECURITY
        endif

        ! computing the Monin-Obukhov length
        diag%monin_obukhov_length(ji,jk) = -theta_lowest_layer*diag%roughness_velocity(ji,jk)**3 &
        /(KARMAN*gravity*w_pert_theta_pert_avg)
      
      enddo
    enddo 
    !$OMP END DO
    !$OMP END PARALLEL

    ! updating the surface flux resistance acting on scalar quantities (moisture and sensible heat)
    if (lsoil) then
      !$OMP PARALLEL
      !$OMP DO PRIVATE(ji,jk)
      do ji=1,nlays
        do jk=1,ncols
          diag%scalar_flux_resistance(ji,jk) = scalar_flux_resistance(diag%roughness_velocity(ji,jk), &
          grid%z_geo_scal(ji+1,jk+1,nlays) - grid%z_geo_w(ji+1,jk+1,nlays+1), &
          grid%roughness_length(ji,jk),diag%monin_obukhov_length(ji,jk))
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    endif

  end subroutine update_sfc_turb_quantities
  
  function roughness_length_from_u10_sea(u10)
  
    ! This function returns the roughness length as a function of the mean wind speed at 10 m above a fully developed sea.

    ! input variable
    real(wp),intent(in) :: u10
    ! output variable
    real(wp)             :: roughness_length_from_u10_sea

    ! local variables
    real(wp) :: swh,period,wavelength

    ! refer to Stensrud,Parameterization schemes (2007),p.130

    ! empirically determined formula for the SWH
    swh = 0.0248_wp*u10**2

    ! empirically determined period of the waves
    period = 0.729_wp*u10

    ! deep-water gravity waves
    wavelength = gravity*period**2/(2._wp*M_PI)

    ! final result
    roughness_length_from_u10_sea = 1200._wp*swh*swh/max(wavelength,EPSILON_SECURITY)**4.5_wp

    ! avoid too small values for stability
    roughness_length_from_u10_sea = max(0.0001_wp,roughness_length_from_u10_sea)
  
  end function roughness_length_from_u10_sea

  function scalar_flux_resistance(roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value)

    ! This function returns the surface flux resistance for scalar quantities.

    ! input variable
    real(wp),intent(in) :: roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value
    ! output variable
    real(wp)             :: scalar_flux_resistance

    ! local variables
    real(wp)             :: used_vertical_height

    ! height of the prandtl layer
    used_vertical_height = min(z_agl,PRANDTL_HEIGHT)

    scalar_flux_resistance = 1._wp/(KARMAN*roughness_velocity_value) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_h(used_vertical_height,monin_obukhov_length_value) &
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
    used_vertical_height = min(z_agl,PRANDTL_HEIGHT)

    momentum_flux_resistance = 1._wp/(KARMAN*roughness_velocity(wind_h_lowest_layer,z_agl,roughness_length_value)) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_m(used_vertical_height,monin_obukhov_length_value))

    ! limitting the result for security
    if (momentum_flux_resistance < 1._wp) then
      momentum_flux_resistance = 1._wp
    endif

  end function momentum_flux_resistance

  function roughness_velocity(wind_speed,z_agl,roughness_length_value)

    ! This function returns the roughness velocity.

    ! input variables
    real(wp), intent(in) :: wind_speed
    real(wp), intent(in) :: z_agl
    real(wp), intent(in) :: roughness_length_value
    ! output variable
    real(wp)             :: roughness_velocity

    ! local variables
    real(wp) :: denominator

    denominator = log(z_agl/roughness_length_value)

    ! security
    if (abs(denominator) < EPSILON_SECURITY) then
      denominator = EPSILON_SECURITY
    endif

    roughness_velocity = wind_speed*KARMAN/denominator

    roughness_velocity = max(EPSILON_SECURITY,roughness_velocity)

  end function roughness_velocity

  function psi_h(z_eff,l)

    ! This is a helper function for the correction to the surface scalar flux resistance for non-neutral conditions.

    ! input variables
    real(wp), intent(in) :: z_eff ! effective height above the surface
    real(wp), intent(in) :: l     ! Monin-Obukhov length
    ! output variable
    real(wp)             :: psi_h

    ! local variables
    real(wp) :: x       ! helper variable
    real(wp) :: l_local ! used to avoid l==0
    
    ! avoiding l==0
    l_local = l
    if (abs(l_local) < EPSILON_SECURITY) then
      l_local = EPSILON_SECURITY
    endif

    ! unstable conditions
    if (l_local<0._wp) then
      x = (1._wp - 15._wp*z_eff/l_local)**0.25_wp
      psi_h = 2._wp*log((1._wp + x**2)/2._wp)     
    ! neutral and stable conditions
    else
      psi_h = -4._wp*z_eff/l_local
    endif
    
  end function psi_h

  function psi_m(z_eff,l)

    ! This is a helper function for the correction to the surface momentum flux resistance for non-neutral conditions.

    ! input variables
    real(wp), intent(in) :: z_eff ! effective height above the surface
    real(wp), intent(in) :: l     ! Monin-Obukhov length
    ! output variable
    real(wp)             :: psi_m

    ! local variables
    real(wp) :: x       ! helper variable
    real(wp) :: l_local ! used to avoid l==0

    ! avoiding l == 0
    l_local = l
    if (abs(l_local) < EPSILON_SECURITY) then
      l_local = EPSILON_SECURITY
    endif

    ! unstable conditions
    if (l_local<0._wp) then
      x = (1._wp - 15._wp*z_eff/l_local)**0.25_wp

      psi_m = 2.0_wp*log((1._wp + x)/2._wp) + log((1._wp + x**2)/2._wp) - 2._wp*atan(x) + M_PI/2._wp
    ! neutral and stable conditions
    else
      psi_m = -4.7_wp*z_eff/l_local
    endif
    
 end function psi_m

end module planetary_boundary_layer





