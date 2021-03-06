! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module explicit_vector_tendencies

  ! This module manages the calculation of the explicit part of the wind tendencies.

  use definitions,              only: t_grid,t_state,t_diag,t_irrev,t_tend,wp
  use inner_product,            only: inner
  use gradient_operators,       only: grad
  use run_nml,                  only: nlins,ncols,nlays,impl_weight,llinear,lcorio
  use constituents_nml,         only: no_of_condensed_constituents
  use vorticities,              only: calc_pot_vort
  use multiplications,          only: scalar_times_vector
  use vorticity_flux,           only: calc_vorticity_flux_term
  use diff_nml,                 only: lmom_diff_h,lmom_diff_v
  use momentum_diff_diss,       only: mom_diff_h,mom_diff_v,simple_dissipation_rate
  use planetary_boundary_layer, only: pbl_wind_tendency
  use surface_nml,              only: lpbl

  implicit none
  
  private
  
  public :: expl_vector_tend
  
  contains

  subroutine expl_vector_tend(state,tend,diag,irrev,grid,rk_step,total_step_counter)
  
    ! This subroutine manages the calculation of the explicit part of the wind tendencies.
  
    type(t_state), intent(in)    :: state              ! state to use for calculating the tendencies
    type(t_tend),  intent(inout) :: tend               ! the tendency
    type(t_diag),  intent(inout) :: diag               ! diagnostic properties
    type(t_irrev), intent(inout) :: irrev              ! irreversible quantities
    type(t_grid),  intent(in)    :: grid               ! grid properties
    integer,       intent(in)    :: rk_step            ! Runge-Kutta step
    integer,       intent(in)    :: total_step_counter ! time step counter of the model integration
    
    ! local variables
    real(wp) :: old_hor_pgrad_weight  ! old time step pressure gradient weight
    real(wp) :: new_hor_pgrad_weight  ! new time step pressure gradient weight
    real(wp) :: old_weight,new_weight ! Runge-Kutta weights
    
    new_hor_pgrad_weight = 0.5_wp + impl_weight
    old_hor_pgrad_weight = 1._wp - new_hor_pgrad_weight
     
    ! momentum advection
    if ((rk_step==2 .or. total_step_counter==0) .and. ((.not. llinear) .or. lcorio)) then
      ! calculating the mass flux density
      call scalar_times_vector(state%rho(:,:,:,no_of_condensed_constituents+1),state%wind_u,state%wind_v,state%wind_w, &
      diag%u_placeholder,diag%v_placeholder,diag%w_placeholder)
      ! calculating the potential vorticity
      call calc_pot_vort(state,diag,grid)
      ! calculating the potential voritcity flux term
      call calc_vorticity_flux_term(diag,grid)
      
      if (.not. llinear) then
        ! Kinetic energy is prepared for the gradient term of the Lamb transformation.
        call inner(state%wind_u,state%wind_v,state%wind_w,state%wind_u,state%wind_v,state%wind_w,diag%v_squared,grid)
        ! taking the gradient of the kinetic energy
        call grad(diag%v_squared,diag%v_squared_grad_x,diag%v_squared_grad_y,diag%v_squared_grad_z,grid)
      endif
    endif
    
    ! momentum diffusion and dissipation (only updated at the first RK step and if advection is updated as well)
    if (rk_step==1) then
      ! horizontal momentum diffusion
      if (lmom_diff_h) then
        call mom_diff_h(state,diag,irrev,grid)
      endif
      ! vertical momentum diffusion
      if (lmom_diff_v) then
        call mom_diff_v(state,diag,irrev,grid)
      endif
      
      ! planetary boundary layer
      if (lpbl) then
        call pbl_wind_tendency(state,diag,irrev,grid)
      endif
      
      ! dissipation
      if (lmom_diff_h .or. lpbl) then
        call simple_dissipation_rate(state,irrev,grid)
      endif
    endif
    
    new_weight = 1._wp
    if (rk_step==2) then
      new_weight = 0.5_wp
    endif
    old_weight = 1._wp - new_weight
    
    ! adding up the explicit wind tendencies
    
    ! x-direction
    !$omp parallel workshare
    tend%wind_u = old_weight*tend%wind_u &
    + new_weight*( &
    ! old time step pressure gradient component
    old_hor_pgrad_weight*diag%p_grad_acc_old_u &
    ! new time step pressure gradient component
    - new_hor_pgrad_weight*(diag%p_grad_acc_neg_nl_u + diag%p_grad_acc_neg_l_u) &
    ! momentum advection
    - 0.5_wp*diag%v_squared_grad_x + diag%pot_vort_tend_x & 
    ! momentum diffusion
    + irrev%mom_diff_tend_x)
    !$omp end parallel workshare
    
    ! y-direction
    !$omp parallel workshare
    tend%wind_v = old_weight*tend%wind_v &
    + new_weight*( &
    ! old time step pressure gradient component
    old_hor_pgrad_weight*diag%p_grad_acc_old_v &
    ! new time step pressure gradient component
    - new_hor_pgrad_weight*(diag%p_grad_acc_neg_nl_v + diag%p_grad_acc_neg_l_v) &
    ! momentum advection
    - 0.5_wp*diag%v_squared_grad_y + diag%pot_vort_tend_y &
    ! momentum diffusion
    + irrev%mom_diff_tend_y)
    !$omp end parallel workshare
    
    ! z-direction
    !$omp parallel workshare
    tend%wind_w = old_weight*tend%wind_w &
    + new_weight*( &
    ! old time step component of the pressure gradient acceleration
    -(1._wp - impl_weight) &
    *(diag%p_grad_acc_neg_nl_w + diag%p_grad_acc_neg_l_w) &
    ! momentum advection
    - 0.5_wp*diag%v_squared_grad_z + diag%pot_vort_tend_z &
    ! momentum diffusion
    + irrev%mom_diff_tend_z)
    !$omp end parallel workshare
    
  end subroutine expl_vector_tend

end module explicit_vector_tendencies










