! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module pressure_gradient

  ! This module manages the handling of the explicit component of the pressure gradient.

  use gradient_operators, only: grad
  use run_nml,            only: wp
  use definitions,        only: t_state,t_diag,t_grid
  use multiplications,    only: scalar_times_vector
  use dictionary,         only: spec_heat_capacities_p_gas

  implicit none

  private
  
  public :: manage_pressure_gradient
  
  contains

  subroutine manage_pressure_gradient(state,diag,grid,lfirst)
  
    ! This subroutine manages the calculation of the pressure gradient.
  
    type(t_state), intent(in)    :: state  ! state to work with
    type(t_diag),  intent(inout) :: diag   ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid   ! model grid
    logical,       intent(in)    :: lfirst ! true for the lfirst model step
    
    ! local variables
    real(wp) :: c_p ! as usual
    
    c_p = spec_heat_capacities_p_gas(0)
  
    ! saving the old pressure gradient acceleration before it is overwritten with the new one
    if (.not. lfirst) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      diag%p_grad_acc_old_u = -diag%p_grad_acc_neg_nl_u - diag%p_grad_acc_neg_l_u
      diag%p_grad_acc_old_v = -diag%p_grad_acc_neg_nl_v - diag%p_grad_acc_neg_l_v
      diag%p_grad_acc_old_w = -diag%p_grad_acc_neg_nl_w - diag%p_grad_acc_neg_l_w
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! the nonlinear pressure gradient term
    ! calculating the gradient of the perturbed Exner pressure
    call grad(state%exner_pert,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v,diag%p_grad_acc_neg_nl_w,grid)
    ! calculating the full potential temperature
    !$OMP PARALLEL
    !$OMP WORKSHARE
    diag%scalar_placeholder = c_p*(grid%theta_bg + state%theta_pert)
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! multiplying the perturbed Exner pressure gradient by the full potential temperature
    call scalar_times_vector(diag%scalar_placeholder,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v, &
    diag%p_grad_acc_neg_nl_w,diag%p_grad_acc_neg_nl_u,diag%p_grad_acc_neg_nl_v,diag%p_grad_acc_neg_nl_w)
    
    ! the linear pressure gradient term
    !$OMP PARALLEL
    !$OMP WORKSHARE
    diag%scalar_placeholder = c_p*state%theta_pert
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! multiplying the background Exner pressure gradient by the perturbed potential temperature
    call scalar_times_vector(diag%scalar_placeholder,grid%exner_bg_grad_u,grid%exner_bg_grad_v, &
    grid%exner_bg_grad_w,diag%p_grad_acc_neg_l_u,diag%p_grad_acc_neg_l_v,diag%p_grad_acc_neg_l_w)
    
    ! At the first step, the "old" pressure gradient acceleration is saved for the first time.
    if (lfirst) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      diag%p_grad_acc_old_u = -diag%p_grad_acc_neg_nl_u - diag%p_grad_acc_neg_l_u
      diag%p_grad_acc_old_v = -diag%p_grad_acc_neg_nl_v - diag%p_grad_acc_neg_l_v
      diag%p_grad_acc_old_w = -diag%p_grad_acc_neg_nl_w - diag%p_grad_acc_neg_l_w
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
  end subroutine manage_pressure_gradient

end module pressure_gradient








