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
    
    ! calculating the divergence of the horizontal wind field
    call div_h(state%wind_u,state%wind_v,diag%scalar_placeholder,grid)
    
    ! computing the relevant diffusion coefficient
    call hori_div_viscosity(state,diag,diag%scalar_placeholder,irrev,grid)
    
    ! multiplying the divergence by the diffusion coefficient acting on divergent movements
    call scalar_times_scalar(irrev%viscosity_coeff_div,diag%scalar_placeholder,diag%scalar_placeholder)
    
    ! calculating the horizontal gradient of the divergence
    call grad_hor(diag%scalar_placeholder,irrev%mom_diff_tend_x,irrev%mom_diff_tend_y,irrev%mom_diff_tend_z,grid)
    
    ! calculating the diffusion coefficient acting on rotational movements
    call hori_curl_viscosity(diag,irrev,grid)
  
  end subroutine mom_diff_h

  subroutine mom_diff_v(state,diag,irrev,grid)
  
    ! This subroutine handles vertical momentum diffusion.
    
    type(t_state), intent(in)    :: state ! state with which to calculate the horizontal diffusion
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji,jk,jl ! loop indices
    real(wp) :: flux_resistance,wind_speed_lowest_layer,z_agl,roughness_length ! variables needed for the surface friction
    real(wp) :: layer_thickness,monin_obukhov_length_value,wind_rescale_factor ! variables needed for the surface friction
    	
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
















