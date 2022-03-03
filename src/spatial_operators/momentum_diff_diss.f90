! This ! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/L-GAME

module momentum_diff_diss

  ! This module handles momentum diffusion and dissipation.
  
  use definitions,              only: t_grid,t_diag,t_irrev,t_state
  use divergence_operators,     only: div_h,add_vertical_div
  use gradient_operators,       only: grad_hor,grad_vert_cov
  use run_nml,                  only: nlins,ncols,nlays,wp
  use diff_nml,                 only: h_prandtl
  use inner_product,            only: inner
  use derived_quantities,       only: density_gas
  use effective_diff_coeffs,    only: hori_div_viscosity,vert_vert_mom_viscosity,hori_curl_viscosity
  use multiplications,          only: scalar_times_scalar
  use planetary_boundary_layer, only: momentum_flux_resistance
  use bc_nml,                   only: lperiodic
  use vorticities,              only: rel_vort
  
  implicit none
  
  private
  
  public :: mom_diff_h
  public :: mom_diff_v
  public :: simple_dissipation_rate
  
  contains
  
  subroutine mom_diff_h(state,diag,irrev,grid)
  
    ! This subroutine handles horizontal momentum diffusion.
    
    type(t_state), intent(in)    :: state ! state with which to calculate the horizontal diffusion
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: upper_index,lower_index ! vertical interpolation indices
    integer :: ji,jk,jl                ! loop indices
    
    ! Preparation of kinematic properties of the wind field
    ! -----------------------------------------------------
    ! calculating the divergence of the horizontal wind field
    call div_h(state%wind_u,state%wind_v,diag%scalar_placeholder,grid)
    ! calculating the relative vorticity of the wind field
    call rel_vort(state,diag,grid)
    
    ! Computing the necessary diffusion coefficients
    ! ----------------------------------------------
    ! computing the relevant diffusion coefficient
    call hori_div_viscosity(state,diag,diag%scalar_placeholder,irrev,grid)
    ! calculating the diffusion coefficient acting on rotational movements
    call hori_curl_viscosity(state,diag,irrev,grid)
    
    ! Computing the gradient of divergence component
    ! ----------------------------------------------
    ! multiplying the divergence by the diffusion coefficient acting on divergent movements
    call scalar_times_scalar(irrev%viscosity_coeff_div,diag%scalar_placeholder,diag%scalar_placeholder)
    ! calculating the horizontal gradient of the divergence
    call grad_hor(diag%scalar_placeholder,irrev%mom_diff_tend_x,irrev%mom_diff_tend_y,irrev%mom_diff_tend_z,grid)
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,upper_index,lower_index)
    do ji=1,nlins
      do jk=1,ncols+1
        do jl=1,nlays
          diag%u_placeholder(ji,jk,jl) = (irrev%viscosity_coeff_curl_dual(ji+1,jk,jl)*diag%zeta_z(ji+1,jk,jl) - &
          irrev%viscosity_coeff_curl_dual(ji,jk,jl)*diag%zeta_z(ji,jk,jl))/grid%dy_dual(ji,jk,jl)
          
          ! terrain-following correction
          diag%u_placeholder(ji,jk,jl) = diag%u_placeholder(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,upper_index,lower_index)
    do ji=1,nlins+1
      do jk=1,ncols
        do jl=1,nlays
          diag%v_placeholder(ji,jk,jl) = (irrev%viscosity_coeff_curl_dual(ji,jk+1,jl)*diag%zeta_z(ji,jk+1,jl) - &
          irrev%viscosity_coeff_curl_dual(ji,jk,jl)*diag%zeta_z(ji,jk,jl))/grid%dx_dual(ji,jk,jl)
          
          ! terrain-following correction
          diag%v_placeholder(ji,jk,jl) = diag%v_placeholder(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! adding up the two components and dividing by the averaged density
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jl=1,nlays
        do jk=2,ncols
          irrev%mom_diff_tend_x(ji,jk,jl) = (irrev%mom_diff_tend_x(ji,jk,jl) + diag%u_placeholder(ji,jk,jl)) &
          /(0.5_wp*(density_gas(state,ji,jk-1,jl) + density_gas(state,ji,jk,jl)))
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          irrev%mom_diff_tend_x(ji,1,jl) = (irrev%mom_diff_tend_x(ji,1,jl) + diag%u_placeholder(ji,1,jl)) &
          /(0.5_wp*(density_gas(state,ji,ncols,jl) + density_gas(state,ji,1,jl)))
        endif
        irrev%mom_diff_tend_x(ji,ncols+1,jl) = irrev%mom_diff_tend_x(ji,1,jl)
        
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do jk=1,ncols
      do jl=1,nlays
        do ji=2,nlins
          irrev%mom_diff_tend_y(ji,jk,jl) = (irrev%mom_diff_tend_y(ji,jk,jl) + diag%v_placeholder(ji,jk,jl)) &
          /(0.5_wp*(density_gas(state,ji-1,jk,jl) + density_gas(state,ji,jk,jl)))
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          irrev%mom_diff_tend_y(1,jk,jl) = (irrev%mom_diff_tend_y(1,jk,jl) + diag%v_placeholder(1,jk,jl)) &
          /(0.5_wp*(density_gas(state,nlins,jk,jl) + density_gas(state,1,jk,jl)))
        endif
        irrev%mom_diff_tend_y(nlins+1,jk,jl) = irrev%mom_diff_tend_y(1,jk,jl)
        
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine mom_diff_h

  subroutine mom_diff_v(state,diag,irrev,grid)
  
    ! This subroutine handles vertical momentum diffusion.
    
    type(t_state), intent(in)    :: state ! state with which to calculate the horizontal diffusion
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    real(wp) :: flux_resistance,wind_speed_lowest_layer,z_agl,roughness_length ! variables needed for the surface friction
    real(wp) :: layer_thickness,monin_obukhov_length_value,wind_rescale_factor ! variables needed for the surface friction
    integer  :: ji,jk,jl                                                       ! loop indices
    
    ! 1.) vertical diffusion of horizontal velocity
    ! ---------------------------------------------
    ! calculating the vertical gradients of the velocity components
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl)
    do jl=2,nlays
      diag%du_dz(:,:,jl) = (state%wind_u(:,:,jl-1) - state%wind_u(:,:,jl))/(grid%z_u(:,:,jl-1) - grid%z_u(:,:,jl))
      diag%dv_dz(:,:,jl) = (state%wind_v(:,:,jl-1) - state%wind_v(:,:,jl))/(grid%z_v(:,:,jl-1) - grid%z_v(:,:,jl))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! extrapolation to the TOA
    !$OMP PARALLEL
    !$OMP WORKSHARE
    diag%du_dz(:,:,1) = diag%du_dz(:,:,2)
    diag%dv_dz(:,:,1) = diag%dv_dz(:,:,2)
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! calculation at the surface
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=2,ncols
      diag%du_dz(:,jk,nlays+1) = state%wind_u(:,jk,nlays)/(grid%z_u(:,jk,nlays) &
      - 0.5_wp*(grid%z_w(:,jk-1,nlays+1)+grid%z_w(:,jk,nlays+1)))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! periodic boundary conditions
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      diag%du_dz(:,1,nlays+1) = state%wind_u(:,1,nlays)/(grid%z_u(:,1,nlays) &
      - 0.5_wp*(grid%z_w(:,ncols,nlays+1)+grid%z_w(:,1,nlays+1)))
      diag%du_dz(:,ncols+1,nlays+1) = diag%du_dz(:,1,nlays+1)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji)
    do ji=2,nlins
      diag%dv_dz(ji,:,nlays+1) = state%wind_v(ji,:,nlays)/(grid%z_v(ji,:,nlays) &
      - 0.5_wp*(grid%z_w(ji-1,:,nlays+1)+grid%z_w(ji,:,nlays+1)))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! periodic boundary conditions
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      diag%dv_dz(1,:,nlays+1) = state%wind_v(1,:,nlays)/(grid%z_v(1,:,nlays) &
      - 0.5_wp*(grid%z_w(nlins,:,nlays+1)+grid%z_w(1,:,nlays+1)))
      diag%dv_dz(nlins+1,:,nlays+1) = diag%dv_dz(1,:,nlays+1)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! calculating the acceleration
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jl=1,nlays
        do jk=2,ncols
          irrev%mom_diff_tend_x(ji,jk,jl) = irrev%mom_diff_tend_x(ji,jk,jl) &
          + (irrev%vert_hor_viscosity_u(ji,jk,jl)*diag%du_dz(ji,jk,jl) &
          - irrev%vert_hor_viscosity_u(ji,jk,jl+1)*diag%du_dz(ji,jk,jl+1)) &
          /(0.5_wp*(grid%z_w(ji,jk-1,jl) + grid%z_w(ji,jk,jl) - grid%z_w(ji,jk-1,jl+1) - grid%z_w(ji,jk,jl+1))) &
          /(0.5_wp*(density_gas(state,ji,jk-1,jl) + density_gas(state,ji,jk,jl)))
        enddo
      
        ! periodic boundary conditions
        if (lperiodic) then
          irrev%mom_diff_tend_x(ji,1,jl) = irrev%mom_diff_tend_x(ji,1,jl) &
          + (irrev%vert_hor_viscosity_u(ji,1,jl)*diag%du_dz(ji,1,jl) &
          - irrev%vert_hor_viscosity_u(ji,1,jl+1)*diag%du_dz(ji,1,jl+1)) &
          /(0.5_wp*(grid%z_w(ji,ncols,jl) + grid%z_w(ji,1,jl) - grid%z_w(ji,ncols,jl+1) - grid%z_w(ji,1,jl+1))) &
          /(0.5_wp*(density_gas(state,ji,ncols,jl) + density_gas(state,ji,1,jl)))
          irrev%mom_diff_tend_x(ji,ncols+1,jl) = irrev%mom_diff_tend_x(ji,1,jl)
        endif
      
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jl)
    do jk=1,ncols
      do jl=1,nlays
        do ji=2,nlins
          irrev%mom_diff_tend_y(ji,jk,jl) = irrev%mom_diff_tend_y(ji,jk,jl) &
          + (irrev%vert_hor_viscosity_v(ji,jk,jl)*diag%dv_dz(ji,jk,jl) &
          - irrev%vert_hor_viscosity_v(ji,jk,jl+1)*diag%dv_dz(ji,jk,jl+1)) &
          /(0.5_wp*(grid%z_w(ji-1,jk,jl) + grid%z_w(ji,jk,jl) - grid%z_w(ji-1,jk,jl+1) - grid%z_w(ji,jk,jl+1))) &
          /(0.5_wp*(density_gas(state,ji-1,jk,jl) + density_gas(state,ji,jk,jl)))
        enddo
      
        ! periodic boundary conditions
        if (lperiodic) then
          irrev%mom_diff_tend_y(1,jk,jl) = irrev%mom_diff_tend_y(1,jk,jl) &
          + (irrev%vert_hor_viscosity_v(1,jk,jl)*diag%dv_dz(1,jk,jl) &
          - irrev%vert_hor_viscosity_v(1,jk,jl+1)*diag%dv_dz(1,jk,jl+1)) &
          /(0.5_wp*(grid%z_w(nlins,jk,jl) + grid%z_w(1,jk,jl) - grid%z_w(nlins,jk,jl+1) - grid%z_w(1,jk,jl+1))) &
          /(0.5_wp*(density_gas(state,nlins,jk,jl) + density_gas(state,1,jk,jl)))
          irrev%mom_diff_tend_y(nlins+1,jk,jl) = irrev%mom_diff_tend_y(1,jk,jl)
        endif
      
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! 2.) vertical diffusion of vertical velocity
    ! -------------------------------------------
    ! resetting the placeholder field
    !$OMP PARALLEL
    !$OMP WORKSHARE
    diag%scalar_placeholder = 0._wp
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! computing something like dw/dz
    call add_vertical_div(state%wind_w,diag%scalar_placeholder,grid)
    ! computing and multiplying by the respective diffusion coefficient
    call vert_vert_mom_viscosity(state,diag,irrev,grid)
    ! taking the second derivative to compute the diffusive tendency
    call grad_vert_cov(diag%scalar_placeholder,irrev%mom_diff_tend_z,grid)


    ! 3.) horizontal diffusion of vertical velocity
    ! ---------------------------------------------
    ! the diffusion coefficient is the same as the one for vertical diffusion of horizontal velocity
    ! averaging the vertical velocity vertically to cell centers, using the inner product weights
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl)
    do jl=1,nlays
      diag%scalar_placeholder(:,:,jl) = &
      grid%inner_product_weights(:,:,jl,5)*state%wind_w(:,:,jl) &
      + grid%inner_product_weights(:,:,jl,6)*state%wind_w(:,:,jl+1)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! computing the horizontal gradient of the vertical velocity field
    call grad_hor(diag%scalar_placeholder,diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
    ! multiplying by the already computed diffusion coefficient
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl)
    do jl=1,nlays
      diag%u_placeholder(:,:,jl) = 0.5_wp*(irrev%vert_hor_viscosity_u(:,:,jl) + irrev%vert_hor_viscosity_u(:,:,jl+1)) &
      *diag%u_placeholder(:,:,jl)
      diag%v_placeholder(:,:,jl) = 0.5_wp*(irrev%vert_hor_viscosity_v(:,:,jl) + irrev%vert_hor_viscosity_v(:,:,jl+1)) &
      *diag%v_placeholder(:,:,jl)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! the divergence of the diffusive flux density results in the diffusive acceleration
    call div_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder,grid)
    ! vertically averaging the divergence to half levels and dividing by the density
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=2,nlays
          irrev%mom_diff_tend_z(ji,jk,jl) = irrev%mom_diff_tend_z(ji,jk,jl) + &
          0.5_wp*(diag%scalar_placeholder(ji,jk,jl-1) + diag%scalar_placeholder(ji,jk,jl))
          ! dividing by the density
          irrev%mom_diff_tend_z(ji,jk,jl) = irrev%mom_diff_tend_z(ji,jk,jl)/ &
          (0.5_wp*(density_gas(state,ji,jk,jl-1) + density_gas(state,ji,jk,jl)))
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! 4.) interaction of the horizontal wind with the surface
    ! -------------------------------------------------------

    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,wind_speed_lowest_layer,layer_thickness,roughness_length,monin_obukhov_length_value, &
    !$OMP flux_resistance,wind_rescale_factor)
    do ji=1,nlins
        do jk=2,ncols

          ! averaging some quantities to the vector point
          wind_speed_lowest_layer = 0.5_wp*(diag%v_squared(ji,jk-1,nlays)**0.5_wp + diag%v_squared(ji,jk,nlays)**0.5_wp)
          z_agl = grid%z_u(ji,jk,nlays) - 0.5_wp*(grid%z_w(ji,jk-1,nlays+1) + grid%z_w(ji,jk,nlays+1))
          layer_thickness = 0.5_wp*(grid%z_w(ji,jk-1,nlays) + grid%z_w(ji,jk,nlays)) &
          - 0.5_wp*(grid%z_w(ji,jk-1,nlays+1) + grid%z_w(ji,jk,nlays+1))
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
          irrev%mom_diff_tend_x(ji,jk,nlays) = irrev%mom_diff_tend_x(ji,jk,nlays) - &
          wind_rescale_factor*state%wind_u(ji,jk,nlays)/flux_resistance/layer_thickness
          
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
        
          ! averaging some quantities to the vector point
          wind_speed_lowest_layer = 0.5_wp*(diag%v_squared(ji,ncols,nlays)**0.5_wp + diag%v_squared(ji,1,nlays)**0.5_wp)
          z_agl = grid%z_u(ji,1,nlays) - 0.5_wp*(grid%z_w(ji,ncols,nlays+1) + grid%z_w(ji,1,nlays+1))
          layer_thickness = 0.5_wp*(grid%z_w(ji,ncols,nlays) + grid%z_w(ji,1,nlays)) &
          - 0.5_wp*(grid%z_w(ji,ncols,nlays+1) + grid%z_w(ji,1,nlays+1))
          roughness_length = 0.5_wp*(grid%roughness_length(ji,ncols) + grid%roughness_length(ji,1))
          monin_obukhov_length_value = 0.5_wp*(diag%monin_obukhov_length(ji,ncols) &
          + diag%monin_obukhov_length(ji,1))

          ! calculating the flux resistance at the vector point
          flux_resistance = momentum_flux_resistance(wind_speed_lowest_layer,z_agl,roughness_length,monin_obukhov_length_value)

          ! rescaling the wind if the lowest wind vector is above the height of the Prandtl layer
          wind_rescale_factor = 1.0_wp
          if (z_agl>h_prandtl) then
            wind_rescale_factor = log(h_prandtl/roughness_length)/log(z_agl/roughness_length)
          endif

          ! adding the momentum flux into the surface as an acceleration
          irrev%mom_diff_tend_x(ji,1,nlays) = irrev%mom_diff_tend_x(ji,1,nlays) - &
          wind_rescale_factor*state%wind_u(ji,1,nlays)/flux_resistance/layer_thickness
          
          irrev%mom_diff_tend_x(ji,ncols+1,nlays) = irrev%mom_diff_tend_x(ji,1,nlays)
          
        endif
        
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,wind_speed_lowest_layer,layer_thickness,roughness_length,monin_obukhov_length_value, &
    !$OMP flux_resistance,wind_rescale_factor)
    do jk=1,ncols
      do ji=2,nlins

          ! averaging some quantities to the vector point
          wind_speed_lowest_layer = 0.5_wp*(diag%v_squared(ji-1,jk,nlays)**0.5_wp + diag%v_squared(ji,jk,nlays)**0.5_wp)
          z_agl = grid%z_v(ji,jk,nlays) - 0.5_wp*(grid%z_w(ji-1,jk,nlays+1) + grid%z_w(ji,jk,nlays+1))
          layer_thickness = 0.5_wp*(grid%z_w(ji-1,jk,nlays) + grid%z_w(ji,jk,nlays)) &
          - 0.5_wp*(grid%z_w(ji-1,jk,nlays+1) + grid%z_w(ji,jk,nlays+1))
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
          irrev%mom_diff_tend_y(ji,jk,nlays) = irrev%mom_diff_tend_y(ji,jk,nlays) - &
          wind_rescale_factor*state%wind_v(ji,jk,nlays)/flux_resistance/layer_thickness
          
        enddo

        ! periodic boundary conditions
        if (lperiodic) then
          ! averaging some quantities to the vector point
          wind_speed_lowest_layer = 0.5_wp*(diag%v_squared(nlins,jk,nlays)**0.5_wp + diag%v_squared(1,jk,nlays)**0.5_wp)
          z_agl = grid%z_v(1,jk,nlays) - 0.5_wp*(grid%z_w(nlins,jk,nlays+1) + grid%z_w(1,jk,nlays+1))
          layer_thickness = 0.5_wp*(grid%z_w(nlins,jk,nlays) + grid%z_w(1,jk,nlays)) &
          - 0.5_wp*(grid%z_w(nlins,jk,nlays+1) + grid%z_w(1,jk,nlays+1))
          roughness_length = 0.5_wp*(grid%roughness_length(nlins,jk) + grid%roughness_length(1,jk))
          monin_obukhov_length_value = 0.5_wp*(diag%monin_obukhov_length(nlins,jk) &
          + diag%monin_obukhov_length(1,jk))

          ! calculating the flux resistance at the vector point
          flux_resistance = momentum_flux_resistance(wind_speed_lowest_layer,z_agl,roughness_length,monin_obukhov_length_value)

          ! rescaling the wind if the lowest wind vector is above the height of the Prandtl layer
          wind_rescale_factor = 1.0_wp
          if (z_agl>h_prandtl) then
            wind_rescale_factor = log(h_prandtl/roughness_length)/log(z_agl/roughness_length)
          endif

          ! adding the momentum flux into the surface as an acceleration
          irrev%mom_diff_tend_y(1,jk,nlays) = irrev%mom_diff_tend_y(1,jk,nlays) - &
          wind_rescale_factor*state%wind_v(1,jk,nlays)/flux_resistance/layer_thickness
          
          irrev%mom_diff_tend_y(nlins+1,jk,nlays) = irrev%mom_diff_tend_y(1,jk,nlays)
          
        endif
        
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine mom_diff_v
  
  subroutine simple_dissipation_rate(state,irrev,grid)
  
    ! This subroutine calculates a simplified dissipation rate.
    
    ! input arguments and output
    type(t_state), intent(in)    :: state ! the state with which to calculate the dissipation rates
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! calculating the inner product of the momentum diffusion acceleration and the wind
    call inner(state%wind_u,state%wind_v,state%wind_w, &
    irrev%mom_diff_tend_x,irrev%mom_diff_tend_y,irrev%mom_diff_tend_z,irrev%heating_diss,grid)
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          irrev%heating_diss(ji,jk,jl) = -density_gas(state,ji,jk,jl)*irrev%heating_diss(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine simple_dissipation_rate

end module momentum_diff_diss
















