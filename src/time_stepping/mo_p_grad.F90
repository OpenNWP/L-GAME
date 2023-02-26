! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_p_grad
  
  ! This module manages the handling of the explicit component of the pressure gradient.
  
  use mo_constants,          only: c_d_p
  use mo_run_nml,            only: theta_adv_order,luse_bg_state
  use mo_gradient_operators, only: grad_vert,grad_hor
  use mo_definitions,        only: t_state,t_diag,t_grid,wp
  use mo_multiplications,    only: scalar_times_vector_h,scalar_times_vector_h2,scalar_times_vector_v,scalar_times_vector_v2, &
                                   theta_v_adv_3rd_order
  use mo_constituents_nml,   only: n_condensed_constituents
  
  implicit none
  
  contains
  
  subroutine manage_pressure_gradient(state,diag,grid,lfirst)
    
    ! This subroutine manages the calculation of the pressure gradient.
    
    type(t_state), intent(in)    :: state  ! state to work with
    type(t_diag),  intent(inout) :: diag   ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid   ! model grid
    logical,       intent(in)    :: lfirst ! true for the lfirst model step
    
    ! local variables
    integer :: use_bg_switch ! switch set to one when using the hydrostatic background state, zero otherwise
    
    use_bg_switch = 0
    if (luse_bg_state) then
      use_bg_switch = 1
    endif
    
    ! saving the old pressure gradient acceleration before it is overwritten with the new one
    if (.not. lfirst) then
      !$omp parallel workshare
      diag%p_grad_acc_old_u = -diag%p_grad_acc_neg_nl_u - use_bg_switch*diag%p_grad_acc_neg_l_u
      diag%p_grad_acc_old_v = -diag%p_grad_acc_neg_nl_v - use_bg_switch*diag%p_grad_acc_neg_l_v
      !$omp end parallel workshare
    endif
    
    ! the nonlinear pressure gradient term
    ! calculating the gradient of the perturbed Exner pressure
    call grad_vert(state%exner_pert,diag%p_grad_acc_neg_nl_w,grid)
    call grad_hor(state%exner_pert,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v,diag%p_grad_acc_neg_nl_w,grid)
    
    ! calculating the full virtual potential temperature
    !$omp parallel workshare
    diag%scalar_placeholder = grid%theta_v_bg + state%theta_v_pert
    !$omp end parallel workshare
    ! multiplying the perturbed Exner pressure gradient by the full virtual potential temperature
    if (theta_adv_order==2) then
      call scalar_times_vector_h2(diag%scalar_placeholder,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v)
    elseif (theta_adv_order==3) then
      call theta_v_adv_3rd_order(state,diag,grid)
      !$omp parallel workshare
      diag%p_grad_acc_neg_nl_u = diag%theta_v_u*diag%p_grad_acc_neg_nl_u
      diag%p_grad_acc_neg_nl_v = diag%theta_v_v*diag%p_grad_acc_neg_nl_v
      !$omp end parallel workshare
    endif
    call scalar_times_vector_v2(diag%scalar_placeholder,diag%p_grad_acc_neg_nl_w)
    
    ! the linear pressure gradient term
    if (luse_bg_state) then
      !$omp parallel workshare
      diag%scalar_placeholder = state%theta_v_pert
      !$omp end parallel workshare
      ! multiplying the background Exner pressure gradient by the perturbed virtual potential temperature
      if (theta_adv_order==2) then
        call scalar_times_vector_h(diag%scalar_placeholder,grid%exner_bg_grad_u,grid%exner_bg_grad_v, &
                                   diag%p_grad_acc_neg_l_u,diag%p_grad_acc_neg_l_v)
      elseif (theta_adv_order==3) then
        call theta_v_adv_3rd_order(state,diag,grid)
        !$omp parallel workshare
        diag%p_grad_acc_neg_l_u = diag%theta_v_u*grid%exner_bg_grad_u
        diag%p_grad_acc_neg_l_v = diag%theta_v_v*grid%exner_bg_grad_v
        !$omp end parallel workshare
      endif
      call scalar_times_vector_v(diag%scalar_placeholder,grid%exner_bg_grad_w,diag%p_grad_acc_neg_l_w)
    endif
    
    !$omp parallel workshare
    diag%p_grad_decel_factor = c_d_p*state%rho(:,:,:,n_condensed_constituents+1)/ &
                                          sum(state%rho(:,:,:,1:n_condensed_constituents+1),4)
    !$omp end parallel workshare
    call scalar_times_vector_h2(diag%p_grad_decel_factor,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v)
    call scalar_times_vector_v2(diag%p_grad_decel_factor,diag%p_grad_acc_neg_nl_w)
    if (luse_bg_state) then
      call scalar_times_vector_h2(diag%p_grad_decel_factor,diag%p_grad_acc_neg_l_u,diag%p_grad_acc_neg_l_v)
      call scalar_times_vector_v2(diag%p_grad_decel_factor,diag%p_grad_acc_neg_l_w)
    endif
    
    ! At the very first step of the model integration, the "old" pressure gradient acceleration is saved for the first time.
    if (lfirst) then
      !$omp parallel workshare
      diag%p_grad_acc_old_u = -diag%p_grad_acc_neg_nl_u - use_bg_switch*diag%p_grad_acc_neg_l_u
      diag%p_grad_acc_old_v = -diag%p_grad_acc_neg_nl_v - use_bg_switch*diag%p_grad_acc_neg_l_v
      !$omp end parallel workshare
    endif
    
  end subroutine manage_pressure_gradient
  
  subroutine calc_pressure_grad_condensates_v(state,diag,grid)
    
    ! This subroutine computes the correction to the vertical pressure gradient acceleration due to condensates.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    !$omp parallel workshare
    diag%p_grad_decel_factor = state%rho(:,:,:,n_condensed_constituents+1) &
                                          /sum(state%rho(:,:,:,1:n_condensed_constituents+1),4) - 1._wp
    !$omp end parallel workshare
    call scalar_times_vector_v(diag%p_grad_decel_factor,grid%gravity_m_v,diag%p_grad_condensates_w)
    
  end subroutine calc_pressure_grad_condensates_v

end module mo_p_grad








