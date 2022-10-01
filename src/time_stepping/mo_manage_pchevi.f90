! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_manage_pchevi

  ! In this module, the RKHEVI time stepping is managed.

  use mo_definitions,            only: t_grid,t_state,t_diag,t_tend,t_bc,wp
  use mo_run_nml,                only: dtime,ny,nx
  use mo_pgrad,                  only: manage_pressure_gradient,calc_pressure_grad_condensates_v
  use mo_scalar_tend_expl,       only: scalar_tend_expl
  use mo_vector_tend_expl,       only: vector_tend_expl
  use mo_column_solvers,         only: three_band_solver_ver,three_band_solver_gen_densities
  use mo_boundaries,             only: update_boundaries
  use mo_derived,                only: temperature_diagnostics
  use mo_constituents_nml,       only: n_constituents
  use mo_surface_nml,            only: lsfc_sensible_heat_flux,lsfc_phase_trans,lpbl
  use mo_pbl,                    only: update_sfc_turb_quantities
  use mo_bc_nml,                 only: lperiodic
  use mo_manage_radiation_calls, only: update_rad_fluxes

  implicit none

  contains
  
  subroutine pchevi(state_old,state_new,tend,bc,grid,diag,total_step_counter,lrad_update,t_0)
  
    ! This subroutine manages the predictor-corrector HEVI time stepping.
    
    type(t_state), intent(inout) :: state_old          ! the state at the old timestep
    type(t_state), intent(inout) :: state_new          ! the state at the new timestep
    type(t_tend),  intent(inout) :: tend               ! the tendency
    type(t_bc),    intent(inout) :: bc                 ! boundary conditions
    type(t_grid),  intent(inout) :: grid               ! the grid of the model
    type(t_diag),  intent(inout) :: diag               ! diagnostic quantities
    integer,       intent(in)    :: total_step_counter ! time step counter
    logical,       intent(in)    :: lrad_update        ! radiation update switch
    real(wp),      intent(in)    :: t_0                ! Unix time
    
    ! local variables
    integer :: rk_step ! index of the Runge-Kutta step
    
    ! Preparations
    ! ------------
    
    ! diagnosing the temperature
    call temperature_diagnostics(state_old,diag,grid)
    
    ! updating surface-related turbulence quantities if it is necessary
    if (lsfc_sensible_heat_flux .or. lsfc_phase_trans .or. lpbl) then
      call update_sfc_turb_quantities(state_old,diag,grid)
    endif
    
    ! upating radiation if necessary
    if (lrad_update) then
      call update_rad_fluxes(state_old,grid,diag,t_0)
    endif
      
    ! Loop over the RK substeps
    ! -------------------------
    
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
        call calc_pressure_grad_condensates_v(state_old,diag,grid)
        call vector_tend_expl(state_old,tend,diag,grid,rk_step,total_step_counter)
      endif
      if (rk_step==2) then
        call calc_pressure_grad_condensates_v(state_old,diag,grid)
        call vector_tend_expl(state_new,tend,diag,grid,rk_step,total_step_counter)
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
        call scalar_tend_expl(grid,state_old,state_new,tend,diag,rk_step)
      endif
      if (rk_step==2) then
        call scalar_tend_expl(grid,state_new,state_new,tend,diag,rk_step)
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
        call three_band_solver_gen_densities(state_old,state_new,tend,diag,grid,rk_step)
      endif
      
    enddo
    
    ! calling the boundary conditions subroutine in real-data simulations
    if (.not. lperiodic) then
      call update_boundaries(state_new,bc,(total_step_counter+1)*dtime,grid)
    endif
    
  end subroutine pchevi

end module mo_manage_pchevi






