! This source file is part of the Limited-area GAME version (L-GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module explicit_scalar_tendencies

  ! This module manages the calculation of the explicit component of the scalar tendencies.

  use constants,             only: c_d_p
  use mo_definitions,        only: wp,t_grid,t_state,t_diag,t_tend
  use mo_multiplications,    only: scalar_times_vector_h,scalar_times_vector_h_upstream,scalar_times_vector_v
  use divergence_operators,  only: div_h,div_h_tracers,add_vertical_div
  use run_nml,               only: dtime
  use mo_phase_trans,        only: calc_h2otracers_source_rates
  use constituents_nml,      only: n_condensed_constituents,n_constituents
  use diff_nml,              only: ltemp_diff_h,ltemp_diff_v,lmass_diff_h,lmass_diff_v
  use effective_diff_coeffs, only: temp_diffusion_coeffs,mass_diffusion_coeffs
  use gradient_operators,    only: grad

  implicit none
  
  contains
  
  subroutine expl_scalar_tend(grid,state_scalar,state_vector,tend,diag,rk_step)
  
    ! This subroutine manages the calculation of the explicit part of the scalar tendencies.
  
    type(t_grid),  intent(in)    :: grid         ! model grid
    type(t_state), intent(in)    :: state_scalar ! state from which to use the scalar quantities
    type(t_state), intent(in)    :: state_vector ! state from which to use the wind
    type(t_tend),  intent(inout) :: tend         ! state which will contain the tendencies
    type(t_diag),  intent(inout) :: diag         ! diagnostic quantities
    integer,       intent(in)    :: rk_step      ! RK substep index
    
    ! local variables
    integer  :: jc                         ! loop variable
    real(wp) :: old_weight(n_constituents) ! time stepping weight
    real(wp) :: new_weight(n_constituents) ! time stepping weight
    integer  :: diff_switch                ! diffusion switch
    
    ! setting the time stepping weights
    do jc=1,n_constituents
      new_weight(jc) = 1._wp
        if (rk_step==2 .and. jc/=n_condensed_constituents+1) then
          new_weight(jc) = 0.5_wp
        endif
      old_weight(jc) = 1._wp - new_weight(jc)
    enddo
    
    ! Temperature diffusion gets updated here,but only at the first RK step and if heat conduction is switched on.
    if (ltemp_diff_h) then
      ! Now we need to calculate the temperature diffusion coefficients.
      call temp_diffusion_coeffs(state_scalar,diag,grid)
      ! The diffusion of the temperature depends on its gradient.
      call grad(diag%temperature,diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
      ! Now the diffusive temperature flux density can be obtained.
      call scalar_times_vector_h(diag%scalar_diff_coeff_h,diag%u_placeholder,diag%v_placeholder, &
      diag%flux_density_u,diag%flux_density_v)
      ! The divergence of the diffusive temperature flux density is the diffusive temperature heating.
      call div_h(diag%flux_density_u,diag%flux_density_v,diag%temp_diff_heating,grid)
      ! vertical temperature diffusion
      if (ltemp_diff_v) then
        call scalar_times_vector_v(diag%scalar_diff_coeff_v,diag%w_placeholder,diag%flux_density_w)
        call add_vertical_div(diag%flux_density_w,diag%temp_diff_heating,grid)
      endif
    endif
    
    ! loop over all constituents
    do jc=1,n_constituents
    
      ! explicit mass densities integration
      ! -----------------------------------
      ! calculating the divergence of the mass flux density
      ! main gaseous constituent
      if (jc==n_condensed_constituents+1) then
        call scalar_times_vector_h(state_scalar%rho(:,:,:,jc),state_scalar%wind_u,state_scalar%wind_v, &
                                   diag%flux_density_u,diag%flux_density_v)
        call div_h(diag%flux_density_u,diag%flux_density_v,diag%flux_density_div,grid)
      ! all other constituents
      else
        call scalar_times_vector_h_upstream( &
        state_scalar%rho(:,:,:,jc),state_vector%wind_u,state_vector%wind_v,diag%flux_density_u,diag%flux_density_v)
        call div_h_tracers(diag%flux_density_u,diag%flux_density_v,state_scalar%rho(:,:,:,jc), &
        state_vector%wind_u,state_vector%wind_v,diag%flux_density_div,grid)
      endif

      ! mass diffusion
      diff_switch = 0
      if (lmass_diff_h) then
        diff_switch = 1
        ! firstly, we need to calculate the mass diffusion coeffcients
        if (jc==1) then
          call mass_diffusion_coeffs(state_scalar,diag,grid)
        endif
        ! The diffusion of the mass density depends on its gradient.
        call grad(state_scalar%rho(:,:,:,jc),diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
        ! Now the diffusive mass flux density can be obtained.
        call scalar_times_vector_h(diag%scalar_diff_coeff_h,diag%u_placeholder,diag%v_placeholder, &
        diag%u_placeholder,diag%v_placeholder)
        ! The divergence of the diffusive mass flux density is the diffusive mass source rate.
        call div_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder,grid)
        ! vertical mass diffusion
        if (lmass_diff_v) then
          call scalar_times_vector_v(diag%scalar_diff_coeff_v,diag%w_placeholder,diag%w_placeholder)
          call add_vertical_div(diag%w_placeholder,diag%scalar_placeholder,grid)
        endif
      endif
      
      !$omp parallel workshare
      tend%rho(:,:,:,jc) = old_weight(jc)*tend%rho(:,:,:,jc) &
      + new_weight(jc)*(-diag%flux_density_div)+diff_switch*diag%scalar_placeholder
      !$omp end parallel workshare
    
      ! explicit virtual potential temperature density integration
      ! --------------------------------------------------
      ! calculating the virtual potential temperature density flux
      if (jc==n_condensed_constituents+1) then
        !$omp parallel workshare
        diag%scalar_placeholder = grid%theta_v_bg + state_scalar%theta_v_pert
        !$omp end parallel workshare
        call scalar_times_vector_h(diag%scalar_placeholder,diag%flux_density_u,diag%flux_density_v,&
        diag%flux_density_u,diag%flux_density_v)
        ! calculating the divergence of the virtual potential temperature flux density
        call div_h(diag%flux_density_u,diag%flux_density_v,diag%flux_density_div,grid)
        
        !$omp parallel workshare
        tend%rhotheta_v = -diag%flux_density_div &
        ! diabatic heating rates
        + (diag%heating_diss + diag%radiation_tendency + diag%heat_source_rates &
        + diag%temp_diff_heating) &
        /(c_d_p*(grid%exner_bg+state_scalar%exner_pert))
        !$omp end parallel workshare
        
      endif
      
    enddo
        
  end subroutine
  
  subroutine moisturizer(state,diag,grid)
  
    ! This subroutine manages the calculation of the phase transition rates.
    
    type(t_state),intent(inout) :: state ! the state with which to calculate the phase transition rates
    type(t_diag), intent(inout) :: diag  ! diagnostic quantities
    type(t_grid), intent(in)    :: grid  ! grid properties
      
    if (n_constituents>1) then
    
      ! calculating the source rates
      call calc_h2otracers_source_rates(state,diag,grid)
      
      ! condensates
      !$omp parallel workshare
      state%rho(:,:,:,1:n_condensed_constituents) = state%rho(:,:,:,1:n_condensed_constituents) &
      + dtime*diag%mass_source_rates(:,:,:,1:n_condensed_constituents)
      !$omp end parallel workshare
      
      ! water vapour
      !$omp parallel workshare
      state%rho(:,:,:,n_constituents) = state%rho(:,:,:,n_constituents) &
      + dtime*diag%mass_source_rates(:,:,:,n_constituents-1)
      !$omp end parallel workshare
      
    endif
  
  end subroutine moisturizer

end module explicit_scalar_tendencies

















