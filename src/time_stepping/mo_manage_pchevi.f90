! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_manage_pchevi

  ! In this module, the RKHEVI time stepping is managed.

  use mo_definitions,             only: t_grid,t_state,t_diag,t_irrev,t_tend,t_bc,wp
  use linear_combine_two_states,  only: lin_combination
  use run_nml,                    only: dtime,ny,nx
  use pressure_gradient,          only: manage_pressure_gradient
  use explicit_vector_tendencies, only: expl_vector_tend
  use explicit_scalar_tendencies, only: expl_scalar_tend,moisturizer
  use mo_column_solvers,          only: three_band_solver_ver,three_band_solver_gen_densities
  use boundaries,                 only: update_boundaries
  use mo_derived,                 only: temperature_diagnostics
  use constituents_nml,           only: n_constituents
  use diff_nml,                   only: lmom_diff_v
  use surface_nml,                only: lsfc_sensible_heat_flux,lsfc_phase_trans,lpbl
  use mo_pbl,                     only: update_sfc_turb_quantities
  use bc_nml,                     only: lperiodic
  use mo_manage_radiation_calls,  only: call_radiation

  implicit none

  contains
  
  subroutine pchevi(state_old,state_new,tend,bc,grid,diag,irrev,total_step_counter,lrad_update,t_0)
  
    ! This subroutine manages the predictor-corrector HEVI time stepping.
    
    type(t_state), intent(inout) :: state_old          ! the state at the old timestep
    type(t_state), intent(inout) :: state_new          ! the state at the new timestep
    type(t_tend),  intent(inout) :: tend               ! the tendency
    type(t_bc),    intent(inout) :: bc                 ! boundary conditions
    type(t_grid),  intent(inout) :: grid               ! the grid of the model
    type(t_diag),  intent(inout) :: diag               ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev              ! irreversible quantities
    integer,       intent(in)    :: total_step_counter ! time step counter
    logical,       intent(in)    :: lrad_update        ! radiation update switch
    real(wp),      intent(in)    :: t_0                ! Unix time
    
    ! local variables
    integer :: rk_step ! index of the Runge-Kutta step
    
    ! diagnosing the temperature
    call temperature_diagnostics(state_old,diag,grid)
    
    ! upating radiation if necessary
    if (lrad_update) then
      call call_radiation(state_old,grid,diag,t_0)
    endif
    
    ! updating surface-related turbulence quantities if it is necessary
    if (lsfc_sensible_heat_flux .or. lsfc_phase_trans .or. lpbl) then
      call update_sfc_turb_quantities(state_old,diag,grid)
    endif
    
    ! loop over the two Runge-Kutta steps
    do rk_step=1,2
    
      ! state_old remains unchanged the whole time.
      ! At rk_step==1, state_new contains garbage.
      
      ! 1.) Explicit component of the momentum equation.
      ! ------------------------------------------------
      ! Update of the pressure gradient.
      if (rk_step==1) then
        call manage_pressure_gradient(state_old,diag,grid,total_step_counter==0)
      endif
      ! Only the horizontal momentum is a forward tendency.
      if (rk_step==1) then
        call expl_vector_tend(state_old,tend,diag,irrev,grid,rk_step,total_step_counter)
      endif
      if (rk_step==2) then
        call expl_vector_tend(state_new,tend,diag,irrev,grid,rk_step,total_step_counter)
      endif
      ! time stepping for the horizontal momentum can be directly executed
      !$omp parallel workshare
      state_new%wind_u = state_old%wind_u + dtime*tend%wind_u
      state_new%wind_v = state_old%wind_v + dtime*tend%wind_v
      !$omp end parallel workshare
      ! Horizontal velocity can be considered to be updated from now on.

      ! 2.) Explicit component of the generalized density equations.
      ! ------------------------------------------------------------
      if (rk_step==1) then
        call expl_scalar_tend(grid,state_old,state_new,tend,diag,irrev,rk_step)
      endif
      if (rk_step==2) then
        call expl_scalar_tend(grid,state_new,state_new,tend,diag,irrev,rk_step)
      endif
      
      ! 3.) implicit dynamic vertical solver
      ! ------------------------------------
      if (rk_step==1) then
        call three_band_solver_ver(state_old,state_old,state_new,diag,tend,grid,rk_step)
      endif
      if (rk_step==2) then
        call three_band_solver_ver(state_old,state_new,state_new,diag,tend,grid,rk_step)
      endif
      
      ! 4.) vertical tracer advection
      ! -----------------------------
      if (n_constituents>1) then
        call three_band_solver_gen_densities(state_old,state_new,tend,grid)
      endif
      
    enddo
    
    ! saturation adjustment, calculation of latent heating rates, evaporation at the surface
    if (n_constituents>1) then
      call moisturizer(state_new,diag,irrev,grid)
    endif
    
    ! calling the boundary conditions subroutine in real-data simulations
    if (.not. lperiodic) then
      call update_boundaries(state_new,bc,(total_step_counter+1)*dtime,grid)
    endif
    
  end subroutine pchevi

end module mo_manage_pchevi






