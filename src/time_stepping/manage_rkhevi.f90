! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module manage_rkhevi

  ! In this module, the RKHEVI time stepping is managed.

  use definitions,                only: t_grid,t_state,t_diag,t_irrev,t_tend,t_bc,wp
  use linear_combine_two_states,  only: lin_combination
  use run_nml,                    only: dtime,nlins,ncols
  use pressure_gradient,          only: manage_pressure_gradient
  use explicit_vector_tendencies, only: expl_vector_tend
  use explicit_scalar_tendencies, only: expl_scalar_tend,moisturizer
  use column_solvers,             only: three_band_solver_ver,three_band_solver_gen_densities
  use boundaries,                 only: update_boundaries
  use derived_quantities,         only: temperature_diagnostics
  use constituents_nml,           only: no_of_constituents
  use diff_nml,                   only: lmom_diff_v
  use surface_nml,                only: lsoil
  use planetary_boundary_layer,   only: update_sfc_turb_quantities
  use bc_nml,                     only: lperiodic
  use manage_radiation_calls,     only: call_radiation

  implicit none
  
  private
  
  public :: rkhevi

  contains
  
  subroutine rkhevi(state_old,state_new,tend,bc,grid,diag,irrev,total_step_counter,lrad_update)
    
    type(t_state), intent(inout) :: state_old          ! the state at the old timestep
    type(t_state), intent(inout) :: state_new          ! the state at the new timestep
    type(t_tend),  intent(inout) :: tend               ! the tendency
    type(t_bc),    intent(inout) :: bc                 ! boundary conditions
    type(t_grid),  intent(inout) :: grid               ! the grid of the model
    type(t_diag),  intent(inout) :: diag               ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev              ! irreversible quantities
    integer,       intent(in)    :: total_step_counter ! time step counter
    logical,       intent(in)    :: lrad_update        ! radiation update switch
    
    ! local variables
    integer :: rk_step ! index of the Runge-Kutta step
    
    ! diagnosing the temperature
    call temperature_diagnostics(state_old,diag,grid)
    
    ! upating radiation if necessary
    if (lrad_update) then
      call call_radiation(state_old,grid,diag,irrev)
    endif
    
    ! updating surface-related turbulence quantities if it is necessary
    if (lsoil .or. lmom_diff_v) then
      call update_sfc_turb_quantities(state_old,diag,grid)
    endif
    
    ! loop over the two Runge-Kutta steps
    do rk_step=1,2
    
      ! state_old remains unchanged the whole time.
      ! At rk_step==1, it is state_old==state_new.
      
      ! 1.) Explicit component of the momentum equation.
      ! ------------------------------------------------
      ! Update of the pressure gradient.
      if (rk_step == 1) then
        call manage_pressure_gradient(state_new,diag,grid,total_step_counter==0)
      endif
      ! Only the horizontal momentum is a forward tendency.
      call expl_vector_tend(state_new,tend,diag,irrev,grid,rk_step,total_step_counter)
      ! time stepping for the horizontal momentum can be directly executed
      !$OMP PARALLEL
      !$OMP WORKSHARE
      state_new%wind_u = state_old%wind_u + dtime*tend%wind_u
      state_new%wind_v = state_old%wind_v + dtime*tend%wind_v
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      ! Horizontal velocity can be considered to be updated from now on.

      ! 2.) Explicit component of the generalized density equations.
      ! ------------------------------------------------------------
      call expl_scalar_tend(grid,state_new,tend,diag,irrev,rk_step)
      
      ! 3.) implicit dynamic vertical solver
      ! ------------------------------------
      call three_band_solver_ver(state_old,state_new,diag,tend,grid,rk_step)
      
      ! 3.) vertical tracer advection
      ! -----------------------------
      if (no_of_constituents>1) then
        call three_band_solver_gen_densities(state_old,state_new,tend,grid)
      endif
  
    enddo
    
    ! saturation adjustment, calculation of latent heating rates, evaporation at the surface
    call moisturizer(state_new,diag,irrev,grid)
    
    ! calling the boundary conditions subroutine for real-data simulation
    if (.not. lperiodic) then
      call update_boundaries(bc)
    endif
    
  end subroutine rkhevi

end module manage_rkhevi






