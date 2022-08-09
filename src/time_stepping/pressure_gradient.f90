! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module pressure_gradient

  ! This module manages the handling of the explicit component of the pressure gradient.

  use constants,          only: c_d_p
  use gradient_operators, only: grad
  use run_nml,            only: wp
  use definitions,        only: t_state,t_diag,t_grid
  use multiplications,    only: scalar_times_vector

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
    call grad(state%exner_pert,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v,diag%p_grad_acc_neg_nl_w,grid)
    ! calculating the full virtual potential temperature
    !$omp parallel workshare
    diag%scalar_placeholder = c_d_p*(grid%theta_v_bg + state%theta_v_pert)
    !$omp end parallel workshare
    ! multiplying the perturbed Exner pressure gradient by the full virtual potential temperature
    call scalar_times_vector(diag%scalar_placeholder,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v, &
    diag%p_grad_acc_neg_nl_w,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v,diag%p_grad_acc_neg_nl_w)
    
    ! the linear pressure gradient term
    !$omp parallel workshare
    diag%scalar_placeholder = c_d_p*state%theta_v_pert
    !$omp end parallel workshare
    ! multiplying the background Exner pressure gradient by the perturbed virtual potential temperature
    call scalar_times_vector(diag%scalar_placeholder,grid%exner_bg_grad_u,grid%exner_bg_grad_v, &
    grid%exner_bg_grad_w,diag%p_grad_acc_neg_l_u,diag%p_grad_acc_neg_l_v,diag%p_grad_acc_neg_l_w)
    
    ! At the first step, the "old" pressure gradient acceleration is saved for the first time.
    if (lfirst) then
      !$omp parallel workshare
      diag%p_grad_acc_old_u = -diag%p_grad_acc_neg_nl_u - diag%p_grad_acc_neg_l_u
      diag%p_grad_acc_old_v = -diag%p_grad_acc_neg_nl_v - diag%p_grad_acc_neg_l_v
      diag%p_grad_acc_old_w = -diag%p_grad_acc_neg_nl_w - diag%p_grad_acc_neg_l_w
      !$omp end parallel workshare
    endif
    
  end subroutine manage_pressure_gradient

end module pressure_gradient








