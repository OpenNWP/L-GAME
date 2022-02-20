! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module manage_rkhevi

  ! In this file, the RKHEVI time stepping is managed.

  use definitions,                only: t_grid,t_state,t_diag,wp,t_tend
  use linear_combine_two_states,  only: lin_combination
  use run_nml,                    only: dtime,nlins,ncols
  use pressure_gradient,          only: manage_pressure_gradient
  use explicit_vector_tendencies, only: expl_vector_tend
  use explicit_scalar_tendencies, only: expl_scalar_tend
  use column_solvers,             only: three_band_solver_ver
  use boundaries,                 only: bc

  implicit none
  
  private
  
  public :: rkhevi

  contains
  
  subroutine rkhevi(state_old,state_new,tend,grid,diag,total_step_counter)
    
    type(t_state),  intent(inout) :: state_old          ! the state at the old timestep
    type(t_state),  intent(inout) :: state_new          ! the state at the new timestep
    type(t_tend),   intent(inout) :: tend               ! the tendency
    type(t_grid),   intent(in)    :: grid               ! the grid of the model
    type(t_diag),   intent(inout) :: diag               ! diagnostic quantities
    integer,        intent(in)    :: total_step_counter ! time step counter
    
    ! local variables
    integer  :: rk_step          ! index of the Runge-Kutta step
    
    do rk_step = 1,2
    
      ! state_old remains unchanged the whole time.
      ! At rk_step == 1, it is state_old == state_new.
      
      ! 1.) Explicit component of the momentum equation.
      ! ------------------------------------------------
      ! Update of the pressure gradient.
      if (rk_step == 1) then
        call manage_pressure_gradient(state_new,diag,grid,total_step_counter==0)
      endif
      ! Only the horizontal momentum is a forward tendency.
      call expl_vector_tend(state_new,tend,diag,grid,rk_step,total_step_counter)
      ! time stepping for the horizontal momentum can be directly executed
      state_new%wind_u(2:nlins+1,2:ncols  ,:) = state_old%wind_u(2:nlins+1,2:ncols  ,:) + dtime*tend%wind_u(:,:,:)
      state_new%wind_v(2:nlins  ,2:ncols+1,:) = state_old%wind_v(2:nlins  ,2:ncols+1,:) + dtime*tend%wind_v(:,:,:)
      ! Horizontal velocity can be considered to be updated from now on.

      ! 2.) Explicit component of the generalized density equations.
      ! ------------------------------------------------------------
      call expl_scalar_tend(grid,state_new,tend,diag)
      
      ! 3.) implicit dynamic vertical solver
      ! ------------------------------------
      call three_band_solver_ver(state_old,state_new,tend,grid,rk_step)
  
    enddo
    
    ! calling the boundary conditions subroutine
    call bc()
    
  end subroutine rkhevi

end module manage_rkhevi






