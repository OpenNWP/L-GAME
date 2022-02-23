! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module explicit_scalar_tendencies

  ! This module manages the calculation of the explicit component of the scalar tendencies.

  use definitions,          only: wp,t_grid,t_state,t_diag,t_irrev,t_tend
  use multiplications,      only: scalar_times_vector_h
  use divergence_operators, only: divv_h
  use run_nml,              only: nlins,ncols,dtime
  use phase_trans,          only: calc_h2otracers_source_rates
  use constituents_nml,     only: no_of_condensed_constituents,no_of_constituents,lassume_lte
  use dictionary,           only: spec_heat_capacities_p_gas

  implicit none
  
  private
  
  public :: expl_scalar_tend
  public :: moisturizer
  
  contains
  
  subroutine expl_scalar_tend(grid,state,tend,diag,irrev,rk_step)
  
    type(t_grid),   intent(in)    :: grid       ! model grid
    type(t_state),  intent(in)    :: state      ! state with which to calculate the divergence
    type(t_tend),   intent(inout) :: tend       ! state which will contain the tendencies
    type(t_diag),   intent(inout) :: diag       ! diagnostic quantities
    type(t_irrev),  intent(inout) :: irrev      ! irreversible quantities
    integer,        intent(in)    :: rk_step    ! RK substep index
    
    ! local variables
    integer                       :: j_constituent                   ! loop variable
    real(wp)                      :: old_weight(no_of_constituents)  ! time stepping weight
    real(wp)                      :: new_weight(no_of_constituents)  ! time stepping weight
    real(wp)                      :: c_p                             ! as usual
    
    c_p = spec_heat_capacities_p_gas(0)
    
    ! setting the time stepping weights
    do j_constituent=1,no_of_constituents
      new_weight(j_constituent) = 1._wp
        if (rk_step == 2 .and. j_constituent /= no_of_condensed_constituents+1) then
          new_weight(j_constituent) = 0.5_wp
        endif
      old_weight(j_constituent) = 1._wp - new_weight(j_constituent)
    enddo
    
    ! loop over all constituents
    do j_constituent=1,no_of_constituents
    
      ! explicit mass densities integration
      ! -----------------------------------
      ! calculating the mass flux density
      call scalar_times_vector_h(state%rho(:,:,:,j_constituent),state%wind_u,state%wind_v,diag%u_placeholder,diag%v_placeholder)
      ! calculating the divergence of the mass flux density
      call divv_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder,grid)
      !$OMP PARALLEL
      !$OMP WORKSHARE
      tend%rho(:,:,:,j_constituent) = old_weight(j_constituent)*tend%rho(:,:,:,j_constituent) &
      + new_weight(j_constituent)*(-diag%scalar_placeholder)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    
      ! explicit potential temperature density integration
      ! --------------------------------------------------
      ! calculating the potential temperature density flux
      if (j_constituent == no_of_condensed_constituents+1) then
        diag%scalar_placeholder = grid%theta_bg + state%theta_pert
        call scalar_times_vector_h(diag%scalar_placeholder,diag%u_placeholder,diag%v_placeholder, &
        diag%u_placeholder,diag%v_placeholder)
        ! calculating the divergence of the potential temperature flux density
        call divv_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder,grid)
        !$OMP PARALLEL
        !$OMP WORKSHARE
        tend%rhotheta = -diag%scalar_placeholder &
        ! dissipative heating
        + (irrev%heating_diss + diag%radiation_tendency)/(c_p*(grid%exner_bg+state%exner_pert))
        !$OMP END WORKSHARE
        !$OMP END PARALLEL
      endif
      
      ! explicit temperature density integration for condensates
      ! --------------------------------------------------------
      if (j_constituent <= no_of_condensed_constituents .and. (.not. lassume_lte)) then
        call scalar_times_vector_h(state%condensed_rho_t(:,:,:,j_constituent),state%wind_u,state%wind_v, &
        diag%u_placeholder,diag%v_placeholder)
        call divv_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder,grid)
        !$OMP PARALLEL
        !$OMP WORKSHARE
        tend%condensed_rho_t(:,:,:,j_constituent) = old_weight(j_constituent)*tend%condensed_rho_t(:,:,:,j_constituent) &
        + new_weight(j_constituent)*(-diag%scalar_placeholder)
        !$OMP END WORKSHARE
        !$OMP END PARALLEL
      endif
      
    enddo
        
  end subroutine
  
  subroutine moisturizer(state,diag,irrev,grid)
  
    ! This subroutine manages the calculation of the phase transition rates.
    
    type(t_state), intent(inout) :: state
    type(t_diag), intent(inout)  :: diag
    type(t_irrev), intent(inout) :: irrev
    type(t_grid), intent(in)     :: grid
      
    if (no_of_constituents > 1) then
    
      ! calculating the source rates
      call calc_h2otracers_source_rates(state,diag,irrev,grid)
      
      ! condensates
      !$OMP PARALLEL
      !$OMP WORKSHARE
      state%rho(:,:,:,1:no_of_condensed_constituents) = state%rho(:,:,:,1:no_of_condensed_constituents) &
      + dtime*irrev%mass_source_rates(:,:,:,1:no_of_condensed_constituents)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      ! water vapour
      !$OMP PARALLEL
      !$OMP WORKSHARE
      state%rho(:,:,:,no_of_constituents) = state%rho(:,:,:,no_of_constituents) &
      + dtime*irrev%mass_source_rates(:,:,:,no_of_constituents-1)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      
    endif
  
  end subroutine moisturizer

end module explicit_scalar_tendencies

















