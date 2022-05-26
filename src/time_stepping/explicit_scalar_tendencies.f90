! This source file is part of the Limited-area GAME version (L-GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module explicit_scalar_tendencies

  ! This module manages the calculation of the explicit component of the scalar tendencies.

  use definitions,           only: wp,t_grid,t_state,t_diag,t_irrev,t_tend
  use multiplications,       only: scalar_times_vector_h,scalar_times_vector_h_upstream,scalar_times_vector_v
  use divergence_operators,  only: div_h,div_h_tracers,add_vertical_div
  use run_nml,               only: dtime
  use phase_trans,           only: calc_h2otracers_source_rates
  use constituents_nml,      only: no_of_condensed_constituents,no_of_constituents,lassume_lte
  use dictionary,            only: spec_heat_capacities_p_gas
  use diff_nml,              only: ltemp_diff_h,ltemp_diff_v,ltracer_diff_h,ltracer_diff_v
  use effective_diff_coeffs, only: temp_diffusion_coeffs,mass_diffusion_coeffs
  use gradient_operators,    only: grad

  implicit none
  
  private
  
  public :: expl_scalar_tend
  public :: moisturizer
  
  contains
  
  subroutine expl_scalar_tend(grid,state,tend,diag,irrev,rk_step)
  
    ! This subroutine manages the calculation of the explicit part of the scalar tendencies.
  
    type(t_grid),  intent(in)    :: grid      ! model grid
    type(t_state), intent(in)    :: state     ! state with which to calculate the divergence
    type(t_tend),  intent(inout) :: tend      ! state which will contain the tendencies
    type(t_diag),  intent(inout) :: diag      ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev     ! irreversible quantities
    integer,       intent(in)    :: rk_step   ! RK substep index
    
    ! local variables
    integer  :: j_constituent                  ! loop variable
    real(wp) :: old_weight(no_of_constituents) ! time stepping weight
    real(wp) :: new_weight(no_of_constituents) ! time stepping weight
    real(wp) :: c_p                            ! as usual
    integer  :: diff_switch                    ! diffusion switch
    
    c_p = spec_heat_capacities_p_gas(0)
    
    ! setting the time stepping weights
    do j_constituent=1,no_of_constituents
      new_weight(j_constituent) = 1._wp
        if (rk_step==2 .and. j_constituent/=no_of_condensed_constituents+1) then
          new_weight(j_constituent) = 0.5_wp
        endif
      old_weight(j_constituent) = 1._wp - new_weight(j_constituent)
    enddo
    
    ! Temperature diffusion gets updated here,but only at the first RK step and if heat conduction is switched on.
    if (ltemp_diff_h) then
      ! Now we need to calculate the temperature diffusion coefficients.
      call temp_diffusion_coeffs(state,diag,irrev,grid)
      ! The diffusion of the temperature depends on its gradient.
      call grad(diag%temperature_gas,diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
      ! Now the diffusive temperature flux density can be obtained.
      call scalar_times_vector_h(irrev%scalar_diff_coeff_h,diag%u_placeholder,diag%v_placeholder, &
      diag%flux_density_u,diag%flux_density_v)
      ! The divergence of the diffusive temperature flux density is the diffusive temperature heating.
      call div_h(diag%flux_density_u,diag%flux_density_v,irrev%temp_diff_heating,grid)
      ! vertical temperature diffusion
      if (ltemp_diff_v) then
        call scalar_times_vector_v(irrev%scalar_diff_coeff_v,diag%w_placeholder,diag%flux_density_w)
        call add_vertical_div(diag%flux_density_w,irrev%temp_diff_heating,grid)
      endif
    endif
    
    ! loop over all constituents
    do j_constituent=1,no_of_constituents
    
      ! explicit mass densities integration
      ! -----------------------------------
      ! calculating the divergence of the mass flux density
      ! main gaseous constituent
      if (j_constituent==no_of_condensed_constituents+1) then
        call scalar_times_vector_h(state%rho(:,:,:,j_constituent),state%wind_u,state%wind_v,diag%flux_density_u,diag%flux_density_v)
        call div_h(diag%flux_density_u,diag%flux_density_v,diag%flux_density_div,grid)
      ! all other constituents
      else
        call scalar_times_vector_h_upstream( &
        state%rho(:,:,:,j_constituent),state%wind_u,state%wind_v,diag%flux_density_u,diag%flux_density_v)
        call div_h_tracers(diag%flux_density_u,diag%flux_density_v,state%rho(:,:,:,j_constituent), &
        state%wind_u,state%wind_v,diag%flux_density_div,grid)
      endif

      ! mass diffusion, only for gaseous tracers
      diff_switch = 0
      if (j_constituent>no_of_condensed_constituents+1 .and. ltracer_diff_h) then
        diff_switch = 1
        ! firstly, we need to calculate the mass diffusion coeffcients
        call mass_diffusion_coeffs(state,diag,irrev,grid)
        ! The diffusion of the tracer density depends on its gradient.
        call grad(state%rho(:,:,:,j_constituent),diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
        ! Now the diffusive mass flux density can be obtained.
        call scalar_times_vector_h(irrev%scalar_diff_coeff_h,diag%u_placeholder,diag%v_placeholder, &
        diag%u_placeholder,diag%v_placeholder)
        ! The divergence of the diffusive mass flux density is the diffusive mass source rate.
        call div_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder,grid)
        ! vertical mass diffusion
        if (ltracer_diff_v) then
          call scalar_times_vector_v(irrev%scalar_diff_coeff_v,diag%w_placeholder,diag%w_placeholder)
          call add_vertical_div(diag%w_placeholder,diag%scalar_placeholder,grid)
        endif
      endif
      
      !$OMP PARALLEL
      !$OMP WORKSHARE
      tend%rho(:,:,:,j_constituent) = old_weight(j_constituent)*tend%rho(:,:,:,j_constituent) &
      + new_weight(j_constituent)*(-diag%flux_density_div)+diff_switch*diag%scalar_placeholder
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    
      ! explicit potential temperature density integration
      ! --------------------------------------------------
      ! calculating the potential temperature density flux
      if (j_constituent==no_of_condensed_constituents+1) then
        !$OMP PARALLEL
        !$OMP WORKSHARE
        diag%scalar_placeholder = grid%theta_bg + state%theta_pert
        !$OMP END WORKSHARE
        !$OMP END PARALLEL
        call scalar_times_vector_h(diag%scalar_placeholder,diag%flux_density_u,diag%flux_density_v,&
        diag%flux_density_u,diag%flux_density_v)
        ! calculating the divergence of the potential temperature flux density
        call div_h(diag%flux_density_u,diag%flux_density_v,diag%flux_density_div,grid)
        
        !$OMP PARALLEL
        !$OMP WORKSHARE
        tend%rhotheta = -diag%flux_density_div &
        ! diabatic heating rates
        + (irrev%heating_diss + diag%radiation_tendency + sum(irrev%heat_source_rates,dim=4) &
        + irrev%temp_diff_heating) &
        /(c_p*(grid%exner_bg+state%exner_pert))
        !$OMP END WORKSHARE
        !$OMP END PARALLEL
        
      endif
      
    enddo
        
  end subroutine
  
  subroutine moisturizer(state,diag,irrev,grid)
  
    ! This subroutine manages the calculation of the phase transition rates.
    
    type(t_state),intent(inout) :: state ! the state with which to calculate the phase transition rates
    type(t_diag), intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev),intent(inout) :: irrev ! irreversible quantities (phase transitions are irreversible)
    type(t_grid), intent(in)    :: grid  ! grid properties
      
    if (no_of_constituents>1) then
    
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

















