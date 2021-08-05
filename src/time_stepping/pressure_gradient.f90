! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module pressure_gradient

	! This module manages the handling of the explicit component of the pressure gradient.

	use gradient_operators, only: grad
	use run_nml,            only: nlins,ncols
	use definitions,        only: t_state,t_diag,t_grid,t_bg
	use multiplications,    only: scalar_times_vector_for_gradient

	implicit none

	private
	
	public :: manage_pressure_gradient
	
	contains

	subroutine manage_pressure_gradient(state,diag,bg,grid,first)
	
		type(t_state), intent(in)    :: state ! state to work with
		type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
		type(t_bg),    intent(in)    :: bg    ! the background state
		type(t_grid),  intent(in)    :: grid  ! model grid
		logical,       intent(in)    :: first ! true for the first model step
	
		! saving the old pressure gradient acceleration before it is overwritten with the new one
		if (.not. first) then
			diag%p_grad_acc_old_u(:,:,:) = diag%p_grad_acc_l_u(:,:,:) + diag%p_grad_acc_nl_u(:,:,:)
			diag%p_grad_acc_old_v(:,:,:) = diag%p_grad_acc_l_v(:,:,:) + diag%p_grad_acc_nl_v(:,:,:)
			diag%p_grad_acc_old_w(:,:,:) = diag%p_grad_acc_l_w(:,:,:) + diag%p_grad_acc_nl_w(:,:,:)
		endif
		
		! the linear pressure gradient term
		! calculating the full Exner pressure
		diag%scalar_placeholder(:,:,:) = bg%exner(:,:,:) + state%exner_pert(:,:,:)
		! calculating the gradient of the Exner pressure
		call grad(diag%scalar_placeholder(2:nlins+1,2:ncols+1,:),diag%p_grad_acc_l_u,diag%p_grad_acc_l_v,diag%p_grad_acc_l_w,grid)
		! multiplying the Exner pressure gradient by the potential temperature
		call scalar_times_vector_for_gradient(bg%theta,diag%p_grad_acc_l_u,diag%p_grad_acc_l_v, &
		diag%p_grad_acc_l_w,diag%p_grad_acc_l_u,diag%p_grad_acc_l_v,diag%p_grad_acc_l_w,grid)
		
		! the nonlinear pressure gradient term
		! calculating the gradient of the Exner pressure
		call grad(state%exner_pert(2:nlins+1,2:ncols+1,:),diag%p_grad_acc_nl_u,diag%p_grad_acc_nl_v,diag%p_grad_acc_nl_w,grid)
		! multiplying the perturbed Exner pressure gradient by the perturbed potential temperature
		call scalar_times_vector_for_gradient(state%theta_pert,diag%p_grad_acc_nl_u,diag%p_grad_acc_nl_v, &
		diag%p_grad_acc_nl_w,diag%p_grad_acc_nl_u,diag%p_grad_acc_nl_v,diag%p_grad_acc_nl_w,grid)
		
		! At the first step, the "old" pressure gradient acceleration is saved for the first time.
		if (first) then
			diag%p_grad_acc_old_u(:,:,:) = diag%p_grad_acc_l_u(:,:,:) + diag%p_grad_acc_nl_u(:,:,:)
			diag%p_grad_acc_old_v(:,:,:) = diag%p_grad_acc_l_v(:,:,:) + diag%p_grad_acc_nl_v(:,:,:)
			diag%p_grad_acc_old_w(:,:,:) = diag%p_grad_acc_l_w(:,:,:) + diag%p_grad_acc_nl_w(:,:,:)
		endif
	
	end subroutine manage_pressure_gradient

end module pressure_gradient








