! This source file is part of the Limited-area GAME version (L-GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_pbl
  
  ! This module computes everything related to the planetary boundary layer.
  
  use mo_definitions, only: wp,t_state,t_grid,t_diag
  use mo_run_nml,     only: ny,nx,n_layers,n_levels,dtime
  use mo_constants,   only: EPSILON_SECURITY,M_PI,gravity
  use mo_surface_nml, only: lprog_soil_temp
  use mo_diff_nml,    only: h_prandtl,karman
  use mo_bc_nml,      only: lperiodic
  
  implicit none
  
  contains

  subroutine update_sfc_turb_quantities(state,diag,grid)

    ! This subroutine updates surface-related turbulence quantities.
    
    type(t_state), intent(in)    :: state ! state with which to calculate the turbulence quantities
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(inout) :: grid  ! grid properties

    ! local variables
    real(wp) :: u_lowest_layer          ! wind speed in the lowest model layer
    real(wp) :: u10                     ! wind speed in 10 m height
    real(wp) :: agl                     ! height above ground f the lowest model layer
    real(wp) :: theta_v_lowest_layer    ! virtual potential temperature in the lowest layer
    real(wp) :: theta_v_second_layer    ! virtual potential temperature in the second-lowest layer
    real(wp) :: dz                      ! vertical gridpoint distance
    real(wp) :: dtheta_v_dz             ! vertical gradient of the virtual potential temperature
    real(wp) :: w_pert                  ! vertical velocity perturbation near the surface
    real(wp) :: theta_v_pert            ! virtual potential temperature perturbation near the surface
    real(wp) :: w_pert_theta_v_pert_avg ! correlation between vertical velocity and virtual potential temperature perturbations
    real(wp) :: w_theta_v_corr          ! semi-empirical coefficient for computing w_pert_theta_v_pert_avg
    integer  :: ji                      ! horizontal index
    integer  :: jk                      ! horizontal index

    w_theta_v_corr = 0.2_wp
    
    !$omp parallel do private(ji,jk,u_lowest_layer,u10,agl,theta_v_lowest_layer,theta_v_second_layer,dz,dtheta_v_dz, &
    !$omp w_pert,theta_v_pert,w_pert_theta_v_pert_avg)
    do jk=1,nx
      do ji=1,ny
        agl = grid%z_scalar(ji,jk,n_layers) - grid%z_w(ji,jk,n_levels)

        ! wind speed in the lowest layer
        u_lowest_layer = diag%v_squared(ji,jk,n_layers)**0.5_wp
        
        ! calculating the 10 m wind velocity from the logarithmic wind profile
        u10 = u_lowest_layer*log(10._wp/grid%roughness_length(ji,jk))/log(agl/grid%roughness_length(ji,jk))

        ! only over the sea the roughness length is time-dependant (because of the waves)
        if (grid%land_fraction(ji,jk)==0._wp) then
          ! calculating the roughness length fom the wind velocity
          grid%roughness_length(ji,jk) = roughness_length_from_u10_sea(u10)
        endif

        ! updating the roughness velocity
        diag%roughness_velocity(ji,jk) = roughness_velocity(u_lowest_layer,agl,grid%roughness_length(ji,jk))

        ! theta_v in the lowest layer
        theta_v_lowest_layer = grid%theta_v_bg(ji,jk,n_layers) + state%theta_v_pert(ji,jk,n_layers)
        ! theta_v in the second-lowest layer
        theta_v_second_layer = grid%theta_v_bg(ji,jk,n_layers-1) + state%theta_v_pert(ji,jk,n_layers-1)

        ! delta z
        dz = grid%z_scalar(ji,jk,n_layers-1) - grid%z_scalar(ji,jk,n_layers)

        ! vertical gradient of theta_v
        dtheta_v_dz = (theta_v_second_layer - theta_v_lowest_layer)/dz

        ! the perturbation of the vertical velocity is assumed to be proportional to the 10 m wind speed
        ! times a stability-dependant factor
        w_pert = u10*max(0.001_wp,0.02_wp*(1._wp - dtheta_v_dz/0.01_wp))
        theta_v_pert = -0.2_wp*dtime*w_pert*dtheta_v_dz
        w_pert_theta_v_pert_avg = w_theta_v_corr*w_pert*theta_v_pert

        ! security
        if (abs(w_pert_theta_v_pert_avg)<EPSILON_SECURITY) then
          w_pert_theta_v_pert_avg = EPSILON_SECURITY
        endif

        ! computing the Monin-Obukhov length
        diag%monin_obukhov_length(ji,jk) = -theta_v_lowest_layer*diag%roughness_velocity(ji,jk)**3 &
        /(karman*gravity*w_pert_theta_v_pert_avg)
        
      enddo
    enddo
    !$omp end parallel do

    ! updating the surface flux resistance acting on scalar quantities (moisture and sensible heat)
    if (lprog_soil_temp) then
      !$omp parallel do private(ji,jk)
      do jk=1,nx
        do ji=1,ny
          diag%scalar_flux_resistance(ji,jk) = scalar_flux_resistance(diag%roughness_velocity(ji,jk), &
          grid%z_scalar(ji,jk,n_layers) - grid%z_w(ji,jk,n_levels), &
          grid%roughness_length(ji,jk),diag%monin_obukhov_length(ji,jk))
        enddo
      enddo
      !$omp end parallel do
    endif

  end subroutine update_sfc_turb_quantities
  
  subroutine pbl_wind_tendency(state,diag,grid)
    
    ! This subroutine computes the interaction of the horizontal wind with the surface.
    
    type(t_state), intent(in)    :: state ! state with which to calculate the horizontal diffusion
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    real(wp) :: flux_resistance            ! momentum flux resistance value
    real(wp) :: wind_speed_lowest_layer    ! wind speed in the lowest layer
    real(wp) :: z_agl                      ! height above ground level of the lowest layer
    real(wp) :: roughness_length           ! roughness length
    real(wp) :: layer_thickness            ! thickness of the lowest layer
    real(wp) :: monin_obukhov_length_value ! Monin-Obukhov length
    real(wp) :: wind_rescale_factor        ! factor taking into account the distance of the wind vector from the surface
    integer  :: ji                         ! horizontal index
    integer  :: jk                         ! horizontal index

    !$omp parallel do private(ji,jk,wind_speed_lowest_layer,z_agl,layer_thickness,roughness_length,monin_obukhov_length_value, &
    !$omp flux_resistance,wind_rescale_factor)
    do ji=1,ny
      do jk=2,nx
        
        ! averaging some quantities to the vector point
        wind_speed_lowest_layer = 0.5_wp*(diag%v_squared(ji,jk-1,n_layers)**0.5_wp + diag%v_squared(ji,jk,n_layers)**0.5_wp)
        z_agl = grid%z_u(ji,jk,n_layers) - 0.5_wp*(grid%z_w(ji,jk-1,n_levels) + grid%z_w(ji,jk,n_levels))
        layer_thickness = 0.5_wp*(grid%z_w(ji,jk-1,n_layers) + grid%z_w(ji,jk,n_layers)) &
        - 0.5_wp*(grid%z_w(ji,jk-1,n_levels) + grid%z_w(ji,jk,n_levels))
        roughness_length = 0.5_wp*(grid%roughness_length(ji,jk-1) + grid%roughness_length(ji,jk))
        monin_obukhov_length_value = 0.5_wp*(diag%monin_obukhov_length(ji,jk-1) &
        + diag%monin_obukhov_length(ji,jk))
        
        ! calculating the flux resistance at the vector point
        flux_resistance = momentum_flux_resistance(wind_speed_lowest_layer,z_agl,roughness_length,monin_obukhov_length_value)
        
        ! rescaling the wind if the lowest wind vector is above the height of the Prandtl layer
        wind_rescale_factor = 1.0_wp
        if (z_agl>h_prandtl) then
          wind_rescale_factor = log(h_prandtl/roughness_length)/log(z_agl/roughness_length)
        endif
        
        ! adding the momentum flux into the surface as an acceleration
        diag%mom_diff_tend_x(ji,jk,n_layers) = diag%mom_diff_tend_x(ji,jk,n_layers) - &
        wind_rescale_factor*state%wind_u(ji,jk,n_layers)/flux_resistance/layer_thickness
          
      enddo
      
      ! periodic boundary conditions
      if (lperiodic) then
        
        ! averaging some quantities to the vector point
        wind_speed_lowest_layer = 0.5_wp*(diag%v_squared(ji,nx,n_layers)**0.5_wp + diag%v_squared(ji,1,n_layers)**0.5_wp)
        z_agl = grid%z_u(ji,1,n_layers) - 0.5_wp*(grid%z_w(ji,nx,n_levels) + grid%z_w(ji,1,n_levels))
        layer_thickness = 0.5_wp*(grid%z_w(ji,nx,n_layers) + grid%z_w(ji,1,n_layers)) &
        - 0.5_wp*(grid%z_w(ji,nx,n_levels) + grid%z_w(ji,1,n_levels))
        roughness_length = 0.5_wp*(grid%roughness_length(ji,nx) + grid%roughness_length(ji,1))
        monin_obukhov_length_value = 0.5_wp*(diag%monin_obukhov_length(ji,nx) &
        + diag%monin_obukhov_length(ji,1))
        
        ! calculating the flux resistance at the vector point
        flux_resistance = momentum_flux_resistance(wind_speed_lowest_layer,z_agl,roughness_length,monin_obukhov_length_value)
        
        ! rescaling the wind if the lowest wind vector is above the height of the Prandtl layer
        wind_rescale_factor = 1.0_wp
        if (z_agl>h_prandtl) then
          wind_rescale_factor = log(h_prandtl/roughness_length)/log(z_agl/roughness_length)
        endif
        
        ! adding the momentum flux into the surface as an acceleration
        diag%mom_diff_tend_x(ji,1,n_layers) = diag%mom_diff_tend_x(ji,1,n_layers) - &
        wind_rescale_factor*state%wind_u(ji,1,n_layers)/flux_resistance/layer_thickness
        
        diag%mom_diff_tend_x(ji,nx+1,n_layers) = diag%mom_diff_tend_x(ji,1,n_layers)
        
      endif
    enddo
    !$omp end parallel do

    !$omp parallel do private(ji,jk,wind_speed_lowest_layer,z_agl,layer_thickness,roughness_length,monin_obukhov_length_value, &
    !$omp flux_resistance,wind_rescale_factor)
    do jk=1,nx
      do ji=2,ny

          ! averaging some quantities to the vector point
          wind_speed_lowest_layer = 0.5_wp*(diag%v_squared(ji-1,jk,n_layers)**0.5_wp + diag%v_squared(ji,jk,n_layers)**0.5_wp)
          z_agl = grid%z_v(ji,jk,n_layers) - 0.5_wp*(grid%z_w(ji-1,jk,n_levels) + grid%z_w(ji,jk,n_levels))
          layer_thickness = 0.5_wp*(grid%z_w(ji-1,jk,n_layers) + grid%z_w(ji,jk,n_layers)) &
          - 0.5_wp*(grid%z_w(ji-1,jk,n_levels) + grid%z_w(ji,jk,n_levels))
          roughness_length = 0.5_wp*(grid%roughness_length(ji-1,jk) + grid%roughness_length(ji,jk))
          monin_obukhov_length_value = 0.5_wp*(diag%monin_obukhov_length(ji-1,jk) &
          + diag%monin_obukhov_length(ji,jk))

          ! calculating the flux resistance at the vector point
          flux_resistance = momentum_flux_resistance(wind_speed_lowest_layer,z_agl,roughness_length,monin_obukhov_length_value)

          ! rescaling the wind if the lowest wind vector is above the height of the Prandtl layer
          wind_rescale_factor = 1.0_wp
          if (z_agl>h_prandtl) then
            wind_rescale_factor = log(h_prandtl/roughness_length)/log(z_agl/roughness_length)
          endif

          ! adding the momentum flux into the surface as an acceleration
          diag%mom_diff_tend_y(ji,jk,n_layers) = diag%mom_diff_tend_y(ji,jk,n_layers) - &
          wind_rescale_factor*state%wind_v(ji,jk,n_layers)/flux_resistance/layer_thickness
          
        enddo

        ! periodic boundary conditions
        if (lperiodic) then
          ! averaging some quantities to the vector point
          wind_speed_lowest_layer = 0.5_wp*(diag%v_squared(ny,jk,n_layers)**0.5_wp + diag%v_squared(1,jk,n_layers)**0.5_wp)
          z_agl = grid%z_v(1,jk,n_layers) - 0.5_wp*(grid%z_w(ny,jk,n_levels) + grid%z_w(1,jk,n_levels))
          layer_thickness = 0.5_wp*(grid%z_w(ny,jk,n_layers) + grid%z_w(1,jk,n_layers)) &
          - 0.5_wp*(grid%z_w(ny,jk,n_levels) + grid%z_w(1,jk,n_levels))
          roughness_length = 0.5_wp*(grid%roughness_length(ny,jk) + grid%roughness_length(1,jk))
          monin_obukhov_length_value = 0.5_wp*(diag%monin_obukhov_length(ny,jk) &
          + diag%monin_obukhov_length(1,jk))

          ! calculating the flux resistance at the vector point
          flux_resistance = momentum_flux_resistance(wind_speed_lowest_layer,z_agl,roughness_length,monin_obukhov_length_value)

          ! rescaling the wind if the lowest wind vector is above the height of the Prandtl layer
          wind_rescale_factor = 1.0_wp
          if (z_agl>h_prandtl) then
            wind_rescale_factor = log(h_prandtl/roughness_length)/log(z_agl/roughness_length)
          endif

          ! adding the momentum flux into the surface as an acceleration
          diag%mom_diff_tend_y(1,jk,n_layers) = diag%mom_diff_tend_y(1,jk,n_layers) - &
          wind_rescale_factor*state%wind_v(1,jk,n_layers)/flux_resistance/layer_thickness
          
          diag%mom_diff_tend_y(ny+1,jk,n_layers) = diag%mom_diff_tend_y(1,jk,n_layers)
          
        endif
        
    enddo
    !$omp end parallel do
    
  end subroutine pbl_wind_tendency
  
  function roughness_length_from_u10_sea(u10)
    
    ! This function returns the roughness length as a function of the mean wind speed at 10 m above a fully developed sea.
    
    real(wp), intent(in) :: u10                           ! wind velocity 10 m above the surface
    real(wp)             :: roughness_length_from_u10_sea ! result

    ! local variables
    real(wp) :: swh        ! significant wave height
    real(wp) :: period     ! mean period of the sea surface waves
    real(wp) :: wavelength ! mean wavelength of the sea surface waves

    ! refer to Stensrud,Parameterization schemes (2007),p.130

    ! empirically determined formula for the SWH
    swh = 0.0248_wp*u10**2

    ! empirically determined period of the waves
    period = 0.729_wp*u10

    ! deep-water gravity waves
    wavelength = gravity*period**2/(2._wp*M_PI)

    ! final result
    roughness_length_from_u10_sea = 1200._wp*swh*(swh/max(wavelength,EPSILON_SECURITY))**4.5_wp

    ! avoid too small values for stability
    roughness_length_from_u10_sea = max(0.0001_wp,roughness_length_from_u10_sea)
    
  end function roughness_length_from_u10_sea

  function scalar_flux_resistance(roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value)
    
    ! This function returns the surface flux resistance for scalar quantities.
    
    real(wp), intent(in) :: roughness_velocity_value   ! roughness velocity
    real(wp), intent(in) :: z_agl                      ! height above ground level
    real(wp), intent(in) :: roughness_length_value     ! roughness length
    real(wp), intent(in) :: monin_obukhov_length_value ! Monin-Obukhov length
    real(wp)             :: scalar_flux_resistance     ! result
    
    ! local variables
    real(wp) :: used_vertical_height ! vertical height above ground level actually used in the calculation
    
    ! height of the prandtl layer
    used_vertical_height = min(z_agl,h_prandtl)
    
    scalar_flux_resistance = 1._wp/(karman*roughness_velocity_value) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_h(used_vertical_height,monin_obukhov_length_value) &
    ! interfacial sublayer
    + log(7._wp))
    
    ! limitting the result for security
    if (scalar_flux_resistance<dtime/z_agl) then
      scalar_flux_resistance = dtime/z_agl
    endif
    
  end function scalar_flux_resistance
  
  function momentum_flux_resistance(wind_h_lowest_layer,z_agl,roughness_length_value,monin_obukhov_length_value)
    
    ! This function returns the surface flux resistance for momentum.

    real(wp), intent(in) :: wind_h_lowest_layer        ! wind speed in the lowest model layer
    real(wp), intent(in) :: z_agl                      ! height of the target gridpoint above the surface (m)
    real(wp), intent(in) :: roughness_length_value     ! roughness length  of the surface (m)
    real(wp), intent(in) :: monin_obukhov_length_value ! Monin-Obukhiv-length (m)
    real(wp)             :: momentum_flux_resistance   ! result

    ! local variables
    real(wp) :: used_vertical_height ! vertical height above ground level actually used in the calculation

    ! height of the prandtl layer
    used_vertical_height = min(z_agl,h_prandtl)

    momentum_flux_resistance = 1._wp/(karman*roughness_velocity(wind_h_lowest_layer,z_agl,roughness_length_value)) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_m(used_vertical_height,monin_obukhov_length_value))

    ! limitting the result for security
    if (momentum_flux_resistance<dtime/z_agl) then
      momentum_flux_resistance = dtime/z_agl
    endif
    
  end function momentum_flux_resistance

  function roughness_velocity(wind_speed,z_agl,roughness_length_value)
    
    ! This function returns the roughness velocity.
    
    real(wp), intent(in) :: wind_speed             ! wind speed at a certain height
    real(wp), intent(in) :: z_agl                  ! height at which the wind speed is valid
    real(wp), intent(in) :: roughness_length_value ! roughness length at this point
    real(wp)             :: roughness_velocity     ! the result

    ! local variables
    real(wp) :: denominator ! helper variable

    denominator = log(z_agl/roughness_length_value)

    ! security
    if (abs(denominator)<EPSILON_SECURITY) then
      denominator = EPSILON_SECURITY
    endif

    roughness_velocity = wind_speed*karman/denominator

    roughness_velocity = max(EPSILON_SECURITY,roughness_velocity)

  end function roughness_velocity

  function psi_h(eff,l)

    ! This is a helper function for the correction to the surface scalar flux resistance for non-neutral conditions.
    
    real(wp), intent(in) :: eff   ! effective height above the surface
    real(wp), intent(in) :: l     ! Monin-Obukhov length
    real(wp)             :: psi_h ! the value of the helper function

    ! local variables
    real(wp) :: x       ! helper variable
    real(wp) :: l_local ! used to avoid l==0
    
    ! avoiding l==0
    l_local = l
    if (abs(l_local)<EPSILON_SECURITY) then
      l_local = EPSILON_SECURITY
    endif

    ! unstable conditions
    if (l_local<0._wp) then
      x = (1._wp - 15._wp*eff/l_local)**0.25_wp
      psi_h = 2._wp*log((1._wp + x**2)/2._wp)
    ! neutral and stable conditions
    else
      psi_h = -4._wp*eff/l_local
    endif
    
  end function psi_h

  function psi_m(eff,l)

    ! This is a helper function for the correction to the surface momentum flux resistance for non-neutral conditions.
    
    real(wp), intent(in) :: eff   ! effective height above the surface
    real(wp), intent(in) :: l     ! Monin-Obukhov length
    real(wp)             :: psi_m ! the value of the helper function

    ! local variables
    real(wp) :: x       ! helper variable
    real(wp) :: l_local ! used to avoid l==0

    ! avoiding l==0
    l_local = l
    if (abs(l_local)<EPSILON_SECURITY) then
      l_local = EPSILON_SECURITY
    endif

    ! unstable conditions
    if (l_local<0._wp) then
      x = (1._wp - 15._wp*eff/l_local)**0.25_wp

      psi_m = 2.0_wp*log((1._wp + x)/2._wp) + log((1._wp + x**2)/2._wp) - 2._wp*atan(x) + M_PI/2._wp
    ! neutral and stable conditions
    else
      psi_m = -4.7_wp*eff/l_local
    endif
    
 end function psi_m

end module mo_pbl





