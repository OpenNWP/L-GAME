! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_vector_tend_expl
  
  ! This module manages the calculation of the explicit part of the wind tendencies.
  
  use mo_definitions,        only: t_grid,t_state,t_diag,t_tend,wp
  use mo_constants,          only: impl_thermo_weight
  use mo_inner_product,      only: inner_product
  use mo_gradient_operators, only: grad_hor,grad_vert
  use mo_constituents_nml,   only: n_condensed_constituents
  use mo_vorticities,        only: calc_pot_vort
  use mo_multiplications,    only: scalar_times_vector_h,scalar_times_vector_v
  use mo_vorticity_flux,     only: calc_vorticity_flux_term
  use mo_diff_nml,           only: lmom_diff_h,lmom_diff_v,lmass_diff_h,ltemp_diff_h,lklemp
  use mo_momentum_diff_diss, only: mom_diff_h,mom_diff_v,simple_dissipation_rate
  use mo_eff_diff_coeffs,    only: update_n_squared
  use mo_tke,                only: tke_update
  use mo_pbl,                only: pbl_wind_tendency
  use mo_run_nml,            only: ny,nx,n_layers,llinear,lcorio,n_levels
  use mo_surface_nml,        only: lpbl
  use mo_bc_nml,             only: lfreeslip
  
  implicit none
  
  contains
  
  subroutine vector_tend_expl(state,tend,diag,grid,rk_step,total_step_counter)
    
    ! This subroutine manages the calculation of the explicit part of the wind tendencies.
    
    type(t_state), intent(in)    :: state              ! state to use for calculating the tendencies
    type(t_tend),  intent(inout) :: tend               ! the tendency
    type(t_diag),  intent(inout) :: diag               ! diagnostic properties
    type(t_grid),  intent(in)    :: grid               ! grid properties
    integer,       intent(in)    :: rk_step            ! Runge-Kutta step
    integer,       intent(in)    :: total_step_counter ! time step counter of the model integration
    
    ! local variables
    integer  :: ji                       ! horizontal index
    integer  :: jk                       ! horizontal index
    real(wp) :: density_value            ! individual density value
    real(wp) :: old_hor_pgrad_weight     ! old time step pressure gradient weight
    real(wp) :: current_hor_pgrad_weight ! current time step horizontal pressure gradient weight
    real(wp) :: current_ver_pgrad_weight ! current time step vertical pressure gradient weight
    real(wp) :: old_weight               ! weight of the old predictor-corrector substep
    real(wp) :: new_weight               ! weight of the new predictor-corrector substep
    
    ! momentum advection
    if ((rk_step==2 .or. total_step_counter==0) .and. ((.not. llinear) .or. lcorio)) then
    
      ! calculating the potential vorticity
      call calc_pot_vort(state,diag,grid)
      
      ! calculating the mass flux density
      call scalar_times_vector_h(state%rho(:,:,:,n_condensed_constituents+1),state%wind_u,state%wind_v, &
                                 diag%u_placeholder,diag%v_placeholder)
      call scalar_times_vector_v(state%rho(:,:,:,n_condensed_constituents+1),state%wind_w,diag%w_placeholder)
      ! calculating the vertical mass flux density at the surface for free slip boundary conditions
      if (lfreeslip) then
        !$omp parallel do private(ji,jk,density_value)
        do jk=1,nx
          do ji=1,ny
            density_value = state%rho(ji,jk,n_layers,n_condensed_constituents+1) &
            + (state%rho(ji,jk,n_layers-1,n_condensed_constituents+1)-state%rho(ji,jk,n_layers,n_condensed_constituents+1)) &
            /(grid%z_scalar(ji,jk,n_layers-1)-grid%z_scalar(ji,jk,n_layers)) &
            *(grid%z_w(ji,jk,n_levels)-grid%z_scalar(ji,jk,n_layers))
            diag%w_placeholder(ji,jk,n_levels) = density_value*state%wind_w(ji,jk,n_levels)
          enddo
        enddo
        !$omp end parallel do
      endif
      
      ! calculating the potential voritcity flux term
      call calc_vorticity_flux_term(diag,grid)
      
      ! resetting the vertical mass flux density at the surface
      if (lfreeslip) then
        !$omp parallel workshare
        diag%w_placeholder(:,:,n_levels) = 0._wp
        !$omp end parallel workshare
      endif
      
      if (.not. llinear) then
        ! Kinetic energy is prepared for the gradient term of the Lamb transformation.
        call inner_product(state%wind_u,state%wind_v,state%wind_w,state%wind_u,state%wind_v,state%wind_w,diag%v_squared,grid)
        ! taking the gradient of the kinetic energy
        call grad_vert(diag%v_squared,diag%v_squared_grad_z,grid)
        call grad_hor(diag%v_squared,diag%v_squared_grad_x,diag%v_squared_grad_y,diag%v_squared_grad_z,grid)
      endif
    endif
    
    ! momentum diffusion and dissipation (only updated at the first RK step)
    if (rk_step==1) then
      ! updating the Brunt-Väisälä frequency and the TKE if any diffusion is switched on because it is required for computing the diffusion coefficients
      if (lmom_diff_h .or. lmass_diff_h .or. ltemp_diff_h) then
        call update_n_squared(state,diag,grid)
        call tke_update(state,diag,grid)
      endif
      ! horizontal momentum diffusion
      if (lmom_diff_h) then
        call mom_diff_h(state,diag,grid)
      endif
      ! vertical momentum diffusion
      if (lmom_diff_v) then
        call mom_diff_v(state,diag,grid)
      endif
      ! planetary boundary layer
      if (lpbl) then
        call pbl_wind_tendency(state,diag,grid)
      endif
      ! calculation of the dissipative heating rate
      if (lmom_diff_h .or. lpbl .or. lklemp) then
        call simple_dissipation_rate(state,diag,grid)
      endif
    endif
    
    ! Now the explicit forces are added up.
    new_weight = 1._wp
    if (rk_step==2) then
      new_weight = 0.5_wp
    endif
    old_weight = 1._wp - new_weight
    ! the weights for the pressure gradient
    current_hor_pgrad_weight = 0.5_wp + impl_thermo_weight
    old_hor_pgrad_weight = 1._wp - current_hor_pgrad_weight
    current_ver_pgrad_weight = 1._wp - impl_thermo_weight
    
    ! x-direction
    !$omp parallel workshare
    tend%wind_u = old_weight*tend%wind_u &
    + new_weight*( &
    ! old time step pressure gradient component
    old_hor_pgrad_weight*diag%p_grad_acc_old_u &
    ! new time step pressure gradient component
    - current_hor_pgrad_weight*(diag%p_grad_acc_neg_nl_u + diag%p_grad_acc_neg_l_u) &
    ! momentum advection
    - 0.5_wp*diag%v_squared_grad_x + diag%pot_vort_tend_x &
    ! momentum diffusion
    + diag%mom_diff_tend_x)
    !$omp end parallel workshare
    
    ! y-direction
    !$omp parallel workshare
    tend%wind_v = old_weight*tend%wind_v &
    + new_weight*( &
    ! old time step pressure gradient component
    old_hor_pgrad_weight*diag%p_grad_acc_old_v &
    ! new time step pressure gradient component
    - current_hor_pgrad_weight*(diag%p_grad_acc_neg_nl_v + diag%p_grad_acc_neg_l_v) &
    ! momentum advection
    - 0.5_wp*diag%v_squared_grad_y + diag%pot_vort_tend_y &
    ! momentum diffusion
    + diag%mom_diff_tend_y)
    !$omp end parallel workshare
    
    ! z-direction
    !$omp parallel workshare
    tend%wind_w = old_weight*tend%wind_w &
    + new_weight*( &
    ! old time step component of the pressure gradient acceleration
    -current_ver_pgrad_weight*(diag%p_grad_acc_neg_nl_w + diag%p_grad_acc_neg_l_w) &
    ! momentum advection
    - 0.5_wp*diag%v_squared_grad_z + diag%pot_vort_tend_z &
    ! effect of condensates on the pressure gradient acceleration
    + diag%pressure_grad_condensates_w &
    ! momentum diffusion
    + diag%mom_diff_tend_z)
    !$omp end parallel workshare
    
  end subroutine vector_tend_expl
  
end module mo_vector_tend_expl










