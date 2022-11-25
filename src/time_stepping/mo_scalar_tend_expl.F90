! This source file is part of the Limited-area GAME version (L-GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_scalar_tend_expl

  ! This module manages the calculation of the explicit component of the scalar tendencies.

  use mo_run_nml,              only: dtime,ny,nx,n_layers,theta_adv_order
  use mo_constants,            only: c_d_p,c_d_v
  use mo_definitions,          only: wp,t_grid,t_state,t_diag,t_tend
  use mo_multiplications,      only: scalar_times_vector_h,scalar_times_vector_h_upstream,scalar_times_vector_v
  use mo_divergence_operators, only: div_h,div_h_tracers,add_vertical_div
  use mo_constituents_nml,     only: n_condensed_constituents,n_constituents
  use mo_diff_nml,             only: ltemp_diff_h,ltemp_diff_v,lmass_diff_h,lmass_diff_v
  use mo_eff_diff_coeffs,      only: scalar_diffusion_coeffs
  use mo_gradient_operators,   only: grad_hor,grad_vert
  use mo_inner_product,        only: theta_v_adv_3rd_order
  use mo_derived,              only: c_v_mass_weighted_air

  implicit none
  
  contains
  
  subroutine scalar_tend_expl(grid,state_old,state_new,tend,diag,rk_step)
    
    ! This subroutine manages the calculation of the explicit part of the scalar tendencies.
    
    type(t_grid),  intent(in)           :: grid      ! model grid
    type(t_state), intent(in),   target :: state_old ! state variables at the old predictor-corrector substep
    type(t_state), intent(in),   target :: state_new ! state variables at the new predictor-corrector substep
    type(t_tend),  intent(inout)        :: tend      ! state which will contain the tendencies
    type(t_diag),  intent(inout)        :: diag      ! diagnostic quantities
    integer,       intent(in)           :: rk_step   ! RK substep index
    
    ! local variables
    integer                :: ji                         ! horizontal index
    integer                :: jk                         ! horizontal index
    integer                :: jl                         ! layer index
    integer                :: jc                         ! constituent index
    real(wp)               :: old_weight(n_constituents) ! time stepping weight
    real(wp)               :: new_weight(n_constituents) ! time stepping weight
    type(t_state), pointer :: state_scalar               ! state from which to use the scalar quantities
    
    ! setting the scalar state
    if (rk_step==1) then
      state_scalar => state_old
    else
      state_scalar => state_new
    endif
    
    ! setting the time stepping weights
    do jc=1,n_constituents
      new_weight(jc) = 1._wp
        if (rk_step==2 .and. jc/=n_condensed_constituents+1) then
          new_weight(jc) = 0.5_wp
        endif
      old_weight(jc) = 1._wp - new_weight(jc)
    enddo
    
    ! updating the scalar diffusion coefficient if required
    if (rk_step==1 .and. (lmass_diff_h .or. ltemp_diff_h)) then
      call scalar_diffusion_coeffs(state_scalar,diag,grid)
    endif
    
    ! Temperature diffusion gets updated here,but only at the first RK step and if heat conduction is switched on.
    if (ltemp_diff_h) then
      ! The diffusion of the temperature depends on its gradient.
      call grad_vert(diag%temperature,diag%w_placeholder,grid)
      call grad_hor(diag%temperature,diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
      ! Now the diffusive temperature flux density can be obtained.
      call scalar_times_vector_h(diag%temp_diffusion_coeff_numerical_h,diag%u_placeholder,diag%v_placeholder, &
                                 diag%flux_density_u,diag%flux_density_v)
      ! The divergence of the diffusive temperature flux density is the diffusive temperature heating.
      call div_h(diag%flux_density_u,diag%flux_density_v,diag%temp_diff_heating,grid)
      ! vertical temperature diffusion
      if (ltemp_diff_v) then
        call scalar_times_vector_v(diag%temp_diffusion_coeff_numerical_v,diag%w_placeholder,diag%flux_density_w)
        call add_vertical_div(diag%flux_density_w,diag%temp_diff_heating,grid)
      endif
    endif

    ! Mass diffusion gets updated at the first RK step if required.
    if (lmass_diff_h .and. rk_step==1) then
      ! loop over all constituents
      do jc=1,n_constituents
        ! The diffusion of the mass density depends on its gradient.
        call grad_vert(state_scalar%rho(:,:,:,jc),diag%w_placeholder,grid)
        call grad_hor(state_scalar%rho(:,:,:,jc),diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
        ! Now the diffusive mass flux density can be obtained.
        call scalar_times_vector_h(diag%mass_diffusion_coeff_numerical_h,diag%u_placeholder,diag%v_placeholder, &
                                   diag%flux_density_u,diag%flux_density_v)
        ! The divergence of the diffusive mass flux density is the diffusive mass source rate.
        call div_h(diag%flux_density_u,diag%flux_density_v,diag%mass_diff_tendency(:,:,:,jc),grid)
        ! vertical mass diffusion
        if (lmass_diff_v) then
          call scalar_times_vector_v(diag%mass_diffusion_coeff_numerical_v,diag%w_placeholder,diag%flux_density_w)
          call add_vertical_div(diag%flux_density_w,diag%mass_diff_tendency(:,:,:,jc),grid)
        endif
      enddo
    endif
    
    ! loop over all constituents
    do jc=1,n_constituents
      ! Explicit mass densities integration
      ! -----------------------------------
      ! calculating the divergence of the mass flux density
      ! main gaseous constituent
      if (jc==n_condensed_constituents+1) then
        call scalar_times_vector_h(state_scalar%rho(:,:,:,jc),state_new%wind_u,state_new%wind_v, &
                                   diag%flux_density_u,diag%flux_density_v)
        call div_h(diag%flux_density_u,diag%flux_density_v,diag%flux_density_div,grid)
      ! all other constituents
      else
        call scalar_times_vector_h_upstream(state_scalar%rho(:,:,:,jc),state_new%wind_u,state_new%wind_v, &
                                            diag%flux_density_u,diag%flux_density_v)
        call div_h_tracers(diag%flux_density_u,diag%flux_density_v,state_scalar%rho(:,:,:,jc), &
                           state_new%wind_u,state_new%wind_v,diag%flux_density_div,grid)
      endif
      
      !$omp parallel workshare
      tend%rho(:,:,:,jc) = old_weight(jc)*tend%rho(:,:,:,jc) &
      ! advection
      + new_weight(jc)*(-diag%flux_density_div &
      ! diffusion
      + diag%mass_diff_tendency(:,:,:,jc) &
      ! phase transitions
      + diag%phase_trans_rates(:,:,:,min(jc,n_condensed_constituents+1)))
      !$omp end parallel workshare
      
      ! Explicit virtual potential temperature density integration
      ! ----------------------------------------------------------
      ! calculating the virtual potential temperature density flux
      if (jc==n_condensed_constituents+1) then
        !$omp parallel workshare
        diag%scalar_placeholder = grid%theta_v_bg + state_scalar%theta_v_pert
        !$omp end parallel workshare
        if (theta_adv_order==2) then
          call scalar_times_vector_h(diag%scalar_placeholder,diag%flux_density_u,diag%flux_density_v, &
                                     diag%u_placeholder,diag%v_placeholder)
        elseif (theta_adv_order==3) then
          call theta_v_adv_3rd_order(state_new,diag,grid)
          !$omp parallel workshare
          diag%u_placeholder = diag%theta_v_u*diag%flux_density_u
          diag%v_placeholder = diag%theta_v_v*diag%flux_density_v
          !$omp end parallel workshare
        endif
        ! calculating the divergence of the virtual potential temperature flux density
        call div_h(diag%u_placeholder,diag%v_placeholder,diag%flux_density_div,grid)
        
        !$omp parallel do private(ji,jk,jl)
        do jl=1,n_layers
          do jk=1,nx
            do ji=1,ny
              tend%rhotheta_v(ji,jk,jl) = -diag%flux_density_div(ji,jk,jl) &
              ! diabatic heating rates
              ! weighting factor accounting for condensates
              + c_d_v*state_scalar%rho(ji,jk,jl,jc)/c_v_mass_weighted_air(state_scalar%rho,diag%temperature,ji,jk,jl)*( &
              ! dissipative heating
              + diag%heating_diss(ji,jk,jl) &
              ! tendency due to radiation
              + diag%radiation_tendency(ji,jk,jl) &
              ! tendency due to phase transitions
              + diag%phase_trans_heating_rate(ji,jk,jl) &
              ! tendency due to falling condensates
              + diag%condensates_sediment_heat(ji,jk,jl) &
              ! tendency due to temperature diffusion
              + diag%temp_diff_heating(ji,jk,jl)) &
              /(c_d_p*(grid%exner_bg(ji,jk,jl)+state_scalar%exner_pert(ji,jk,jl))) &
              ! tendency of due to phase transitions and mass diffusion
              + (diag%phase_trans_rates(ji,jk,jl,jc) + diag%mass_diff_tendency(ji,jk,jl,jc)) &
              *diag%scalar_placeholder(ji,jk,jl)
            enddo
          enddo
        enddo
        !$omp end parallel do
        
      endif
      
    enddo
    
  end subroutine scalar_tend_expl

end module mo_scalar_tend_expl

















