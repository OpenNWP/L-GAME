! This ! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_momentum_diff_diss
  
  ! In this module momentum diffusion and dissipation operators are calculated.
  
  use mo_constants,            only: M_PI
  use mo_definitions,          only: t_grid,t_diag,t_state,wp
  use mo_divergence_operators, only: div_h,add_vertical_div
  use mo_gradient_operators,   only: grad_hor,grad_vert
  use mo_run_nml,              only: ny,nx,n_layers,n_levels,toa,dtime
  use mo_diff_nml,             only: h_prandtl,lklemp,klemp_begin_rel,klemp_damp_max
  use mo_constituents_nml,     only: n_constituents,n_condensed_constituents
  use mo_grid_generator,       only: n_damping_levels
  use mo_inner_product,        only: inner_product
  use mo_eff_diff_coeffs,      only: hor_viscosity,vert_vert_mom_viscosity
  use mo_bc_nml,               only: lperiodic
  use mo_vorticities,          only: rel_vort
  
  implicit none
  
  contains
  
  subroutine mom_diff_h(state,diag,grid)
    
    ! This subroutine handles horizontal momentum diffusion.
    
    type(t_state), intent(in)    :: state ! state with which to calculate the horizontal diffusion
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: lower_index       ! lower vertical interpolation index
    integer  :: upper_index       ! upper vertical interpolation index
    real(wp) :: slope             ! tangential slope at an edge
    real(wp) :: vertical_gradient ! vertical gradient needed for computing the diffusion operator
    integer  :: ji                ! horizontal index
    integer  :: jk                ! horizontal index
    integer  :: jl                ! layer index
    
    ! Preparation of kinematic properties of the wind field
    ! -----------------------------------------------------
    ! calculating the divergence of the horizontal wind field
    call div_h(state%wind_u,state%wind_v,diag%scalar_placeholder,grid)
    ! calculating the relative vorticity of the wind field
    call rel_vort(state,diag,grid)
    
    ! Computing the necessary diffusion coefficients
    ! ----------------------------------------------
    ! calculating the effective horizontal kinematic viscosity
    call hor_viscosity(state,diag)
    
    ! Computing the gradient of divergence component
    ! ----------------------------------------------
    ! multiplying the divergence by the diffusion coefficient acting on divergent movements
    !$omp parallel workshare
    diag%scalar_placeholder = diag%viscosity*diag%scalar_placeholder
    !$omp end parallel workshare
    ! calculating the horizontal gradient of the divergence
    call grad_hor(diag%scalar_placeholder,diag%mom_diff_tend_x,diag%mom_diff_tend_y,diag%mom_diff_tend_z,grid)
    
    ! Computing the curl of vorticity component
    ! -----------------------------------------
    !$omp parallel do private(ji,jk,jl,upper_index,lower_index,slope,vertical_gradient)
    do jl=1,n_layers
      do jk=1,nx+1
        do ji=1,ny
          diag%u_placeholder(ji,jk,jl) = (diag%viscosity_coeff_curl_dual(ji+1,jk,jl)*diag%zeta_z(ji+1,jk,jl) - &
          diag%viscosity_coeff_curl_dual(ji,jk,jl)*diag%zeta_z(ji,jk,jl))/grid%dy_dual(ji,jk,jl)
          
          ! terrain-following correction
          slope = (grid%z_area_dual_z(ji+1,jk,jl) - grid%z_area_dual_z(ji,jk,jl))/grid%dy_dual(ji,jk,jl)
          upper_index = max(1,jl-1)
          lower_index = min(n_layers,jl+1)
          ! computing the vertical gradient of the vertical vorticity (times the respective diffusion coefficient)
          ! and averaging it to the vector point
          vertical_gradient = 0.5_wp*(diag%viscosity_coeff_curl_dual(ji+1,jk,upper_index)*diag%zeta_z(ji+1,jk,upper_index) &
          - diag%viscosity_coeff_curl_dual(ji+1,jk,lower_index)*diag%zeta_z(ji+1,jk,lower_index)) &
          /(grid%z_area_dual_z(ji+1,jk,upper_index) - grid%z_area_dual_z(ji+1,jk,lower_index)) &
          + 0.5_wp*(diag%viscosity_coeff_curl_dual(ji,jk,upper_index)*diag%zeta_z(ji,jk,upper_index) &
          - diag%viscosity_coeff_curl_dual(ji,jk,lower_index)*diag%zeta_z(ji,jk,lower_index)) &
          /(grid%z_area_dual_z(ji,jk,upper_index) - grid%z_area_dual_z(ji,jk,lower_index))
          ! adding the terrain-following correction
          diag%u_placeholder(ji,jk,jl) = diag%u_placeholder(ji,jk,jl) - slope*vertical_gradient
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl,upper_index,lower_index,slope,vertical_gradient)
    do jl=1,n_layers
      do jk=1,nx
        do ji=1,ny+1
          diag%v_placeholder(ji,jk,jl) = (diag%viscosity_coeff_curl_dual(ji,jk+1,jl)*diag%zeta_z(ji,jk+1,jl) - &
          diag%viscosity_coeff_curl_dual(ji,jk,jl)*diag%zeta_z(ji,jk,jl))/grid%dx_dual(ji,jk,jl)
          
          ! terrain-following correction
          slope = (grid%z_area_dual_z(ji,jk+1,jl) - grid%z_area_dual_z(ji,jk,jl))/grid%dx_dual(ji,jk,jl)
          upper_index = max(1,jl-1)
          lower_index = min(n_layers,jl+1)
          ! computing the vertical gradient of the vertical vorticity (times the respective diffusion coefficient)
          ! and averaging it to the vector point
          vertical_gradient = 0.5_wp*(diag%viscosity_coeff_curl_dual(ji,jk+1,upper_index)*diag%zeta_z(ji,jk+1,upper_index) &
          - diag%viscosity_coeff_curl_dual(ji,jk+1,lower_index)*diag%zeta_z(ji,jk+1,lower_index)) &
          /(grid%z_area_dual_z(ji,jk+1,upper_index) - grid%z_area_dual_z(ji,jk+1,lower_index)) &
          + 0.5_wp*(diag%viscosity_coeff_curl_dual(ji,jk,upper_index)*diag%zeta_z(ji,jk,upper_index) &
          - diag%viscosity_coeff_curl_dual(ji,jk,lower_index)*diag%zeta_z(ji,jk,lower_index)) &
          /(grid%z_area_dual_z(ji,jk,upper_index) - grid%z_area_dual_z(ji,jk,lower_index))
          ! adding the terrain-following correction
          diag%v_placeholder(ji,jk,jl) = diag%v_placeholder(ji,jk,jl) - slope*vertical_gradient
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! adding up the two components and dividing by the averaged density
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do ji=1,ny
        do jk=2,nx
          diag%mom_diff_tend_x(ji,jk,jl) = (diag%mom_diff_tend_x(ji,jk,jl) + diag%u_placeholder(ji,jk,jl)) &
          /(0.5_wp*(sum(state%rho(ji,jk-1,jl,1:n_condensed_constituents+1)) &
                  + sum(state%rho(ji,jk,jl,1:n_condensed_constituents+1))))
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          diag%mom_diff_tend_x(ji,1,jl) = (diag%mom_diff_tend_x(ji,1,jl) + diag%u_placeholder(ji,1,jl)) &
          /(0.5_wp*(sum(state%rho(ji,nx,jl,1:n_condensed_constituents+1)) &
                  + sum(state%rho(ji,1,jl,1:n_condensed_constituents+1))))
        endif
        diag%mom_diff_tend_x(ji,nx+1,jl) = diag%mom_diff_tend_x(ji,1,jl)
        
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do jk=1,nx
        do ji=2,ny
          diag%mom_diff_tend_y(ji,jk,jl) = (diag%mom_diff_tend_y(ji,jk,jl) + diag%v_placeholder(ji,jk,jl)) &
          /(0.5_wp*(sum(state%rho(ji-1,jk,jl,1:n_condensed_constituents+1)) &
                  + sum(state%rho(ji,jk,jl,1:n_condensed_constituents+1))))
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          diag%mom_diff_tend_y(1,jk,jl) = (diag%mom_diff_tend_y(1,jk,jl) + diag%v_placeholder(1,jk,jl)) &
          /(0.5_wp*(sum(state%rho(ny,jk,jl,1:n_condensed_constituents+1)) &
                  + sum(state%rho(1,jk,jl,1:n_condensed_constituents+1))))
        endif
        diag%mom_diff_tend_y(ny+1,jk,jl) = diag%mom_diff_tend_y(1,jk,jl)
        
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine mom_diff_h
  
  subroutine mom_diff_v(state,diag,grid)
    
    ! This subroutine handles vertical momentum diffusion.
    
    type(t_state), intent(in)    :: state ! state with which to calculate the horizontal diffusion
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji ! horizontal index
    integer :: jk ! horizontal index
    integer :: jl ! layer index
    
    ! 1.) vertical diffusion of horizontal velocity
    ! ---------------------------------------------
    ! calculating the vertical gradients of the velocity components
    !$omp parallel do private(jl)
    do jl=2,n_layers
      diag%du_dz(:,:,jl) = (state%wind_u(:,:,jl-1) - state%wind_u(:,:,jl))/(grid%z_u(:,:,jl-1) - grid%z_u(:,:,jl))
      diag%dv_dz(:,:,jl) = (state%wind_v(:,:,jl-1) - state%wind_v(:,:,jl))/(grid%z_v(:,:,jl-1) - grid%z_v(:,:,jl))
    enddo
    !$omp end parallel do
    ! extrapolation to the TOA
    !$omp parallel workshare
    diag%du_dz(:,:,1) = diag%du_dz(:,:,2)
    diag%dv_dz(:,:,1) = diag%dv_dz(:,:,2)
    !$omp end parallel workshare
    ! calculation at the surface
    !$omp parallel do private(jk)
    do jk=2,nx
      diag%du_dz(:,jk,n_levels) = state%wind_u(:,jk,n_layers)/(grid%z_u(:,jk,n_layers) &
      - 0.5_wp*(grid%z_w(:,jk-1,n_levels)+grid%z_w(:,jk,n_levels)))
    enddo
    !$omp end parallel do
    
    ! periodic boundary conditions
    if (lperiodic) then
      !$omp parallel workshare
      diag%du_dz(:,1,n_levels) = state%wind_u(:,1,n_layers)/(grid%z_u(:,1,n_layers) &
      - 0.5_wp*(grid%z_w(:,nx,n_levels)+grid%z_w(:,1,n_levels)))
      diag%du_dz(:,nx+1,n_levels) = diag%du_dz(:,1,n_levels)
      !$omp end parallel workshare
    endif
    
    !$omp parallel do private(ji)
    do ji=2,ny
      diag%dv_dz(ji,:,n_levels) = state%wind_v(ji,:,n_layers)/(grid%z_v(ji,:,n_layers) &
      - 0.5_wp*(grid%z_w(ji-1,:,n_levels)+grid%z_w(ji,:,n_levels)))
    enddo
    !$omp end parallel do
    
    ! periodic boundary conditions
    if (lperiodic) then
      !$omp parallel workshare
      diag%dv_dz(1,:,n_levels) = state%wind_v(1,:,n_layers)/(grid%z_v(1,:,n_layers) &
      - 0.5_wp*(grid%z_w(ny,:,n_levels)+grid%z_w(1,:,n_levels)))
      diag%dv_dz(ny+1,:,n_levels) = diag%dv_dz(1,:,n_levels)
      !$omp end parallel workshare
    endif
    
    ! calculating the acceleration
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do ji=1,ny
        do jk=2,nx
          diag%mom_diff_tend_x(ji,jk,jl) = diag%mom_diff_tend_x(ji,jk,jl) &
          + (diag%vert_hor_viscosity_u(ji,jk,jl)*diag%du_dz(ji,jk,jl) &
          - diag%vert_hor_viscosity_u(ji,jk,jl+1)*diag%du_dz(ji,jk,jl+1)) &
          /(0.5_wp*(grid%z_w(ji,jk-1,jl) + grid%z_w(ji,jk,jl) - grid%z_w(ji,jk-1,jl+1) - grid%z_w(ji,jk,jl+1))) &
          /(0.5_wp*(sum(state%rho(ji,jk-1,jl,1:n_condensed_constituents+1)) &
                  + sum(state%rho(ji,jk,jl,1:n_condensed_constituents+1))))
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          diag%mom_diff_tend_x(ji,1,jl) = diag%mom_diff_tend_x(ji,1,jl) &
          + (diag%vert_hor_viscosity_u(ji,1,jl)*diag%du_dz(ji,1,jl) &
          - diag%vert_hor_viscosity_u(ji,1,jl+1)*diag%du_dz(ji,1,jl+1)) &
          /(0.5_wp*(grid%z_w(ji,nx,jl) + grid%z_w(ji,1,jl) - grid%z_w(ji,nx,jl+1) - grid%z_w(ji,1,jl+1))) &
          /(0.5_wp*(sum(state%rho(ji,nx,jl,1:n_condensed_constituents+1)) &
                  + sum(state%rho(ji,1,jl,1:n_condensed_constituents+1))))
          diag%mom_diff_tend_x(ji,nx+1,jl) = diag%mom_diff_tend_x(ji,1,jl)
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do jk=1,nx
        do ji=2,ny
          diag%mom_diff_tend_y(ji,jk,jl) = diag%mom_diff_tend_y(ji,jk,jl) &
          + (diag%vert_hor_viscosity_v(ji,jk,jl)*diag%dv_dz(ji,jk,jl) &
          - diag%vert_hor_viscosity_v(ji,jk,jl+1)*diag%dv_dz(ji,jk,jl+1)) &
          /(0.5_wp*(grid%z_w(ji-1,jk,jl) + grid%z_w(ji,jk,jl) - grid%z_w(ji-1,jk,jl+1) - grid%z_w(ji,jk,jl+1))) &
          /(0.5_wp*(sum(state%rho(ji-1,jk,jl,1:n_condensed_constituents+1)) &
                  + sum(state%rho(ji,jk,jl,1:n_condensed_constituents+1))))
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          diag%mom_diff_tend_y(1,jk,jl) = diag%mom_diff_tend_y(1,jk,jl) &
          + (diag%vert_hor_viscosity_v(1,jk,jl)*diag%dv_dz(1,jk,jl) &
          - diag%vert_hor_viscosity_v(1,jk,jl+1)*diag%dv_dz(1,jk,jl+1)) &
          /(0.5_wp*(grid%z_w(ny,jk,jl) + grid%z_w(1,jk,jl) - grid%z_w(ny,jk,jl+1) - grid%z_w(1,jk,jl+1))) &
          /(0.5_wp*(sum(state%rho(ny,jk,jl,1:n_condensed_constituents+1)) &
                  + sum(state%rho(1,jk,jl,1:n_condensed_constituents+1))))
          diag%mom_diff_tend_y(ny+1,jk,jl) = diag%mom_diff_tend_y(1,jk,jl)
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    ! 2.) vertical diffusion of vertical velocity
    ! -------------------------------------------
    ! resetting the placeholder field
    !$omp parallel workshare
    diag%scalar_placeholder = 0._wp
    !$omp end parallel workshare
    ! computing something like dw/dz
    call add_vertical_div(state%wind_w,diag%scalar_placeholder,grid)
    ! computing and multiplying by the respective diffusion coefficient
    call vert_vert_mom_viscosity(state,diag,grid)
    ! taking the second derivative to compute the diffusive tendency
    call grad_vert(diag%scalar_placeholder,diag%mom_diff_tend_z,grid)
    
    ! 3.) horizontal diffusion of vertical velocity
    ! ---------------------------------------------
    ! the diffusion coefficient is the same as the one for vertical diffusion of horizontal velocity
    ! averaging the vertical velocity vertically to cell centers, using the inner product weights
    !$omp parallel do private(jl)
    do jl=1,n_layers
      diag%scalar_placeholder(:,:,jl) = &
      grid%inner_product_weights(5,:,:,jl)*state%wind_w(:,:,jl) &
      + grid%inner_product_weights(6,:,:,jl)*state%wind_w(:,:,jl+1)
    enddo
    !$omp end parallel do
    ! computing the horizontal gradient of the vertical velocity field
    call grad_hor(diag%scalar_placeholder,diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
    ! multiplying by the already computed diffusion coefficient
    !$omp parallel do private(jl)
    do jl=1,n_layers
      diag%u_placeholder(:,:,jl) = 0.5_wp*(diag%vert_hor_viscosity_u(:,:,jl) + diag%vert_hor_viscosity_u(:,:,jl+1)) &
      *diag%u_placeholder(:,:,jl)
      diag%v_placeholder(:,:,jl) = 0.5_wp*(diag%vert_hor_viscosity_v(:,:,jl) + diag%vert_hor_viscosity_v(:,:,jl+1)) &
      *diag%v_placeholder(:,:,jl)
    enddo
    !$omp end parallel do
    
    ! the divergence of the diffusive flux density results in the diffusive acceleration
    call div_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder,grid)
    ! vertically averaging the divergence to half levels and dividing by the density
    !$omp parallel do private(ji,jk,jl)
    do jl=2,n_layers
      do jk=1,nx
        do ji=1,ny
          diag%mom_diff_tend_z(ji,jk,jl) = diag%mom_diff_tend_z(ji,jk,jl) + &
          0.5_wp*(diag%scalar_placeholder(ji,jk,jl-1) + diag%scalar_placeholder(ji,jk,jl))
          ! dividing by the density
          diag%mom_diff_tend_z(ji,jk,jl) = diag%mom_diff_tend_z(ji,jk,jl)/ &
          (0.5_wp*(sum(state%rho(ji,jk,jl-1,1:n_condensed_constituents+1)) &
                 + sum(state%rho(ji,jk,jl,1:n_condensed_constituents+1))))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine mom_diff_v
  
  subroutine simple_dissipation_rate(state,diag,grid)
    
    ! This subroutine calculates a simplified dissipation rate.
    
    type(t_state), intent(in)    :: state ! the state with which to calculate the dissipation rates
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji                   ! horizontal index
    integer  :: jk                   ! horizontal index
    integer  :: jl                   ! layer index
    real(wp) :: damping_start_height ! the height in which the Klemp layer begins
    real(wp) :: z_above_damping      ! height of a given gridpoint above damping_start_height
    real(wp) :: damping_coeff        ! coefficient needed for the Klemp damping layer
    real(wp) :: damping_prefactor    ! coefficient needed for the Klemp damping layer
    
    ! copying the vertical diffusion acceleration into w_placeholder
    !$omp parallel workshare
    diag%w_placeholder = diag%mom_diff_tend_z
    !$omp end parallel workshare
    
    ! heating rate in the swamp layer
    if (lklemp) then
      
      damping_start_height = klemp_begin_rel*toa
      
      !$omp parallel do private(ji,jk,jl,z_above_damping,damping_coeff,damping_prefactor)
      do jl=2,n_damping_levels
        do jk=1,nx
          do ji=1,ny
            z_above_damping = grid%z_w(ji,jk,jl) - damping_start_height
            if (z_above_damping<0._wp) then
              damping_coeff = 0._wp
            else
              damping_coeff = klemp_damp_max*sin(0.5_wp*M_PI*z_above_damping/(toa - damping_start_height))**2
            endif
            damping_prefactor = 1._wp + dtime*damping_coeff
            
            ! adding the acceleration due to the swamp layer
            diag%w_placeholder(ji,jk,jl) = diag%w_placeholder(ji,jk,jl) &
                                           + (state%wind_w(ji,jk,jl)/damping_prefactor - state%wind_w(ji,jk,jl))/dtime
          enddo
        enddo
      enddo
      !$omp end parallel do
    endif
    
    ! calculating the inner product of the momentum diffusion acceleration and the wind
    call inner_product(state%wind_u,state%wind_v,state%wind_w,diag%mom_diff_tend_x,diag%mom_diff_tend_y, &
                       diag%w_placeholder,diag%heating_diss,grid)
    
    !$omp parallel workshare
    diag%heating_diss = -sum(state%rho(:,:,:,1:n_condensed_constituents+1),4)*diag%heating_diss
    !$omp end parallel workshare
    
  end subroutine simple_dissipation_rate
  
end module mo_momentum_diff_diss
















