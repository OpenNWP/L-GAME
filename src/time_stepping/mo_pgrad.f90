! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_pgrad

  ! This module manages the handling of the explicit component of the pressure gradient.

  use mo_constants,          only: c_d_p
  use mo_gradient_operators, only: grad_vert,grad_hor
  use mo_definitions,        only: t_state,t_diag,t_grid,wp
  use mo_multiplications,    only: scalar_times_vector_h,scalar_times_vector_v
  use mo_constituents_nml,   only: n_condensed_constituents

  implicit none
  
  contains

  subroutine manage_pressure_gradient(state,diag,grid,lfirst)
  
    ! This subroutine manages the calculation of the pressure gradient.
  
    type(t_state), intent(in)    :: state  ! state to work with
    type(t_diag),  intent(inout) :: diag   ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid   ! model grid
    logical,       intent(in)    :: lfirst ! true for the lfirst model step
  
    ! saving the old pressure gradient acceleration before it is overwritten with the new one
    if (.not. lfirst) then
      !$omp parallel workshare
      diag%p_grad_acc_old_u = -diag%p_grad_acc_neg_nl_u - diag%p_grad_acc_neg_l_u
      diag%p_grad_acc_old_v = -diag%p_grad_acc_neg_nl_v - diag%p_grad_acc_neg_l_v
      diag%p_grad_acc_old_w = -diag%p_grad_acc_neg_nl_w - diag%p_grad_acc_neg_l_w
      !$omp end parallel workshare
    endif
    
    ! the nonlinear pressure gradient term
    ! calculating the gradient of the perturbed Exner pressure
    call grad_vert(state%exner_pert,diag%p_grad_acc_neg_nl_w,grid)
    call grad_hor(state%exner_pert,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v,diag%p_grad_acc_neg_nl_w,grid)
    ! calculating the full virtual potential temperature
    !$omp parallel workshare
    diag%scalar_placeholder = c_d_p*(grid%theta_v_bg + state%theta_v_pert)
    !$omp end parallel workshare
    ! multiplying the perturbed Exner pressure gradient by the full virtual potential temperature
    call scalar_times_vector_h(diag%scalar_placeholder,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v, &
                               diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v)
    call scalar_times_vector_v(diag%scalar_placeholder,diag%p_grad_acc_neg_nl_w,diag%p_grad_acc_neg_nl_w)
    
    ! the linear pressure gradient term
    !$omp parallel workshare
    diag%scalar_placeholder = c_d_p*state%theta_v_pert
    !$omp end parallel workshare
    ! multiplying the background Exner pressure gradient by the perturbed virtual potential temperature
    call scalar_times_vector_h(diag%scalar_placeholder,grid%exner_bg_grad_u,grid%exner_bg_grad_v, &
                               diag%p_grad_acc_neg_l_u,diag%p_grad_acc_neg_l_v)
    call scalar_times_vector_v(diag%scalar_placeholder,grid%exner_bg_grad_w,diag%p_grad_acc_neg_l_w)
    
    !$omp parallel workshare
    diag%pressure_gradient_decel_factor = state%rho(:,:,:,n_condensed_constituents+1)/ &
                                          sum(state%rho(:,:,:,1:n_condensed_constituents+1),4)
    !$omp end parallel workshare
    call scalar_times_vector_h(diag%pressure_gradient_decel_factor,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v, &
                               diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v)
    call scalar_times_vector_v(diag%pressure_gradient_decel_factor,diag%p_grad_acc_neg_nl_w,diag%p_grad_acc_neg_nl_w)
    call scalar_times_vector_h(diag%pressure_gradient_decel_factor,diag%p_grad_acc_neg_l_u,diag%p_grad_acc_neg_l_v, &
                               diag%p_grad_acc_neg_l_u,diag%p_grad_acc_neg_l_v)
    call scalar_times_vector_v(diag%pressure_gradient_decel_factor,diag%p_grad_acc_neg_l_w,diag%p_grad_acc_neg_l_w)
    
    ! At the very first step of the model integration, the "old" pressure gradient acceleration is saved for the first time.
    if (lfirst) then
      !$omp parallel workshare
      diag%p_grad_acc_old_u = -diag%p_grad_acc_neg_nl_u - diag%p_grad_acc_neg_l_u
      diag%p_grad_acc_old_v = -diag%p_grad_acc_neg_nl_v - diag%p_grad_acc_neg_l_v
      diag%p_grad_acc_old_w = -diag%p_grad_acc_neg_nl_w - diag%p_grad_acc_neg_l_w
      !$omp end parallel workshare
    endif
    
  end subroutine manage_pressure_gradient
  
  subroutine calc_pressure_grad_condensates_v(state,diag,grid)
    
    ! This subroutine computes the correction to the vertical pressure gradient acceleration due to condensates.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(inout) :: grid  ! grid quantities
    
    !$omp parallel workshare
    diag%pressure_gradient_decel_factor = state%rho(:,:,:,n_condensed_constituents+1) &
                                          /sum(state%rho(:,:,:,1:n_condensed_constituents+1),4) - 1._wp
    !$omp end parallel workshare
    call scalar_times_vector_v(diag%pressure_gradient_decel_factor,grid%gravity_m_v,diag%pressure_grad_condensates_w)
  
  end subroutine calc_pressure_grad_condensates_v

end module mo_pgrad








