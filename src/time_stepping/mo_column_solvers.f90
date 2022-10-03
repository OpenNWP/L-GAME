! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_column_solvers

  ! This module contains the implicit vertical routines (implicit part of the HEVI scheme).

  use mo_run_nml,          only: ny,nx,n_layers,n_levels,dtime,toa,impl_weight,partial_impl_weight
  use mo_constituents_nml, only: n_condensed_constituents,n_constituents, &
                                 snow_velocity,rain_velocity,cloud_droplets_velocity,graupel_velocity
  use mo_definitions,      only: t_grid,t_state,t_tend,t_diag,wp
  use mo_diff_nml,         only: lklemp,klemp_damp_max,klemp_begin_rel
  use mo_surface_nml,      only: nsoillays,lprog_soil_temp,lsfc_sensible_heat_flux
  use mo_constants,        only: M_PI,r_d,c_d_v,c_d_p
  use mo_dictionary,       only: c_p_cond

  implicit none
  
  contains
  
  subroutine three_band_solver_ver(state_old,state_new,state_target,diag,tend,grid,rk_step)
  
    ! This subroutine is the main implicit vertical solver.

    ! input arguments and output
    type(t_state), intent(in)    :: state_old    ! state at the old timestep
    type(t_state), intent(inout) :: state_new    ! state at the new timestep
    type(t_state), intent(inout) :: state_target ! state into which to write the results
    type(t_diag),  intent(inout) :: diag         ! diagnostic quantities
    type(t_tend),  intent(inout) :: tend         ! explicit tendencies
    type(t_grid),  intent(in)    :: grid         ! model grid
    integer,       intent(in)    :: rk_step      ! Runge Kutta substep

    ! local variables
    integer  :: soil_switch                           ! soil switch: 0 if soil does not have to be calculated off, 1 if soil has to be calculated
    real(wp) :: c_vector(n_layers-2+nsoillays)        ! needed for the vertical solver
    real(wp) :: d_vector(n_layers-1+nsoillays)        ! needed for the vertical solver
    real(wp) :: e_vector(n_layers-2+nsoillays)        ! needed for the vertical solver
    real(wp) :: r_vector(n_layers-1+nsoillays)        ! needed for the vertical solver
    real(wp) :: rho_expl(n_layers)                    ! explicit mass density
    real(wp) :: rhotheta_v_expl(n_layers)             ! explicit virtual potential temperature density
    real(wp) :: exner_pert_expl(n_layers)             ! explicit Exner pressure perturbation
    real(wp) :: theta_v_pert_expl(n_layers)           ! explicit virtual potential temperature perturbation
    real(wp) :: rho_int_old(n_layers-1)               ! old interface mass density
    real(wp) :: rho_int_expl(n_layers-1)              ! explicit interface mass density
    real(wp) :: theta_v_int_new(n_layers-1)           ! preliminary new virtual potential temperature interface values
    real(wp) :: rho_int_new                           ! new density interface value
    real(wp) :: alpha_old(n_layers)                   ! alpha at the old time step
    real(wp) :: beta_old(n_layers)                    ! beta at the old time step
    real(wp) :: gamma_old(n_layers)                   ! gamma at the old time step
    real(wp) :: alpha_new(n_layers)                   ! alpha at the new time step
    real(wp) :: beta_new(n_layers)                    ! beta at the new time step
    real(wp) :: gamma_new(n_layers)                   ! gamma at the new time step
    real(wp) :: alpha(n_layers)                       ! alpha
    real(wp) :: beta(n_layers)                        ! beta
    real(wp) :: gammaa(n_layers)                      ! gamma
    real(wp) :: damping_start_height                  ! lower boundary height of the Klemp layer
    real(wp) :: damping_coeff                         ! damping coefficient of the Klemp layer
    real(wp) :: above_damping                         ! height above the lower boundary of the damping height
    real(wp) :: t_gas_lowest_layer_old                ! temperature of the gas in the lowest layer of the model atmosphere, old timestep
    real(wp) :: t_gas_lowest_layer_new                ! temperature of the gas in the lowest layer of the model atmosphere, new timestep
    real(wp) :: heat_flux_density_expl(nsoillays)     ! explicit heat_flux_density in the soil
    real(wp) :: solution_vector(n_layers-1+nsoillays) ! vector containing the solution of the linear problem to solve here
    integer  :: ji,jk,jl                              ! loop variables
    
    damping_start_height = klemp_begin_rel*toa
    
    ! calculating the sensible power flux density if soil is switched on
    if (lsfc_sensible_heat_flux) then
      !$omp parallel do private(ji,jk,t_gas_lowest_layer_old,t_gas_lowest_layer_new)
      do ji=1,ny
        do jk=1,nx

          ! gas temperature in the lowest layer
          t_gas_lowest_layer_old = (grid%exner_bg(ji,jk,n_layers)+state_old%exner_pert(ji,jk,n_layers)) &
          *(grid%theta_v_bg(ji,jk,n_layers)+state_old%theta_v_pert(ji,jk,n_layers))
          t_gas_lowest_layer_new = (grid%exner_bg(ji,jk,n_layers)+state_new%exner_pert(ji,jk,n_layers)) &
          *(grid%theta_v_bg(ji,jk,n_layers)+state_new%theta_v_pert(ji,jk,n_layers))

          ! the sensible power flux density
          diag%power_flux_density_sensible(ji,jk) = 0.5_wp*c_d_v*(state_new%rho(ji,jk,n_layers,n_condensed_constituents+1) &
          *(t_gas_lowest_layer_old - state_old%temperature_soil(ji,jk,1)) &
          + state_old%rho(ji,jk,n_layers,n_condensed_constituents+1) &
          *(t_gas_lowest_layer_new - state_new%temperature_soil(ji,jk,1)))/diag%scalar_flux_resistance(ji,jk)

          ! contribution of sensible heat to rhotheta_v
          tend%rhotheta_v(ji,jk,n_layers) = tend%rhotheta_v(ji,jk,n_layers) &
          -grid%area_z(ji,jk,n_levels)*diag%power_flux_density_sensible(ji,jk) &
          /((grid%exner_bg(ji,jk,n_layers)+state_new%exner_pert(ji,jk,n_layers))*c_d_p)/grid%volume(ji,jk,n_layers)
          
        enddo
      enddo
      !$omp end parallel do
    endif
	
    !$omp parallel do private(ji,jk,jl,c_vector,d_vector,e_vector,r_vector,solution_vector, &
    !$omp rho_expl,rhotheta_v_expl,exner_pert_expl,theta_v_pert_expl,rho_int_old, &
    !$omp rho_int_expl,theta_v_int_new,rho_int_new,alpha_old,beta_old,gamma_old,alpha_new, &
    !$omp beta_new,gamma_new,alpha,beta,gammaa,damping_coeff,above_damping,soil_switch)
    do ji=1,ny
      do jk=1,nx
      
        ! determining wether soil needs to be calculated
        soil_switch = 0
        if (lprog_soil_temp .and. grid%is_land(ji,jk)==1) then
          soil_switch=1
        endif
      
        ! explicit quantities
        do jl=1,n_layers
          ! explicit density
          rho_expl(jl) = state_old%rho(ji,jk,jl,n_condensed_constituents+1) &
          + dtime*tend%rho(ji,jk,jl,n_condensed_constituents+1)
          ! explicit virtual potential temperature density
          rhotheta_v_expl(jl) = state_old%rhotheta_v(ji,jk,jl) + dtime*tend%rhotheta_v(ji,jk,jl)
          if (rk_step==1) then
            ! old time step partial derivatives of theta_v and Pi (divided by the volume)
            alpha(jl) = -state_old%rhotheta_v(ji,jk,jl)/state_old%rho(ji,jk,jl,n_condensed_constituents+1)**2 &
            /grid%volume(ji,jk,jl)
            beta(jl)  = 1._wp/state_old%rho(ji,jk,jl,n_condensed_constituents+1)/grid%volume(ji,jk,jl)
            gammaa(jl) = r_d/(c_d_v*state_old%rhotheta_v(ji,jk,jl))* &
            (grid%exner_bg(ji,jk,jl)+state_old%exner_pert(ji,jk,jl))/grid%volume(ji,jk,jl)
          else
            ! old time step partial derivatives of theta_v and Pi
            alpha_old(jl) = -state_old%rhotheta_v(ji,jk,jl)/state_old%rho(ji,jk,jl,n_condensed_constituents+1)**2
            beta_old(jl)  = 1._wp/state_old%rho(ji,jk,jl,n_condensed_constituents+1)
            gamma_old(jl) = r_d/(c_d_v*state_old%rhotheta_v(ji,jk,jl))* &
            (grid%exner_bg(ji,jk,jl)+state_old%exner_pert(ji,jk,jl))
            ! new time step partial derivatives of theta_v and Pi
            alpha_new(jl) = -state_new%rhotheta_v(ji,jk,jl)/state_new%rho(ji,jk,jl,n_condensed_constituents+1)**2
            beta_new(jl)  = 1._wp/state_new%rho(ji,jk,jl,n_condensed_constituents+1)
            gamma_new(jl) = r_d/(c_d_v*state_new%rhotheta_v(ji,jk,jl)) &
            *(grid%exner_bg(ji,jk,jl)+state_new%exner_pert(ji,jk,jl))
            ! interpolation of partial derivatives of theta_v and Pi (divided by the volume)
            alpha(jl) = ((1._wp - partial_impl_weight)*alpha_old(jl)+partial_impl_weight*alpha_new(jl))/grid%volume(ji,jk,jl)
            beta(jl) = ((1._wp - partial_impl_weight)*beta_old (jl)+partial_impl_weight*beta_new(jl))/grid%volume(ji,jk,jl)
            gammaa(jl) = ((1._wp - partial_impl_weight)*gamma_old(jl)+partial_impl_weight*gamma_new(jl))/grid%volume(ji,jk,jl)
          endif
          ! explicit virtual potential temperature perturbation
          theta_v_pert_expl(jl) = state_old%theta_v_pert(ji,jk,jl) &
          + dtime*grid%volume(ji,jk,jl)*(alpha(jl)*tend%rho(ji,jk,jl,n_condensed_constituents+1) &
          + beta(jl)*tend%rhotheta_v(ji,jk,jl))
          ! explicit Exner pressure perturbation
          exner_pert_expl(jl) = state_old%exner_pert(ji,jk,jl)+dtime*grid%volume(ji,jk,jl)*gammaa(jl)*tend%rhotheta_v(ji,jk,jl)
        enddo
        
        ! interface values
        do jl=1,n_layers-1
          rho_int_old(jl) = 0.5_wp*(state_old%rho(ji,jk,jl,n_condensed_constituents+1) &
          + state_old%rho(ji,jk,jl+1,n_condensed_constituents+1))
          rho_int_expl(jl) = 0.5_wp*(rho_expl(jl)+rho_expl(jl+1))
          theta_v_int_new(jl) = 0.5_wp*(state_new%rhotheta_v(ji,jk,jl)/state_new%rho(ji,jk,jl,n_condensed_constituents+1) &
          + state_new%rhotheta_v(ji,jk,jl+1)/state_new%rho(ji,jk,jl+1,n_condensed_constituents+1))
        enddo
      
        ! filling up the coefficient vectors
        do jl=1,n_layers-1
          ! main diagonal
          d_vector(jl) = -theta_v_int_new(jl)**2*(gammaa(jl)+gammaa(jl+1)) &
          + 0.5_wp*(grid%exner_bg(ji,jk,jl)-grid%exner_bg(ji,jk,jl+1)) &
          *(alpha(jl+1)-alpha(jl)+theta_v_int_new(jl)*(beta(jl+1)-beta(jl))) &
          - (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))/(impl_weight*dtime**2*c_d_p*rho_int_old(jl)) &
          *(2._wp/grid%area_z(ji,jk,jl+1)+dtime*state_old%wind_w(ji,jk,jl+1)*0.5_wp &
          *(-1._wp/grid%volume(ji,jk,jl)+1._wp/grid%volume(ji,jk,jl+1)))
          ! right hand side
          r_vector(jl) = -(state_old%wind_w(ji,jk,jl+1)+dtime*tend%wind_w(ji,jk,jl+1))* &
          (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1)) &
          /(impl_weight*dtime**2*c_d_p) &
          + theta_v_int_new(jl)*(exner_pert_expl(jl)-exner_pert_expl(jl+1))/dtime &
          + 0.5_wp/dtime*(theta_v_pert_expl(jl)+theta_v_pert_expl(jl+1))*(grid%exner_bg(ji,jk,jl)-grid%exner_bg(ji,jk,jl+1)) &
          - (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))/(impl_weight*dtime**2*c_d_p) &
          *state_old%wind_w(ji,jk,jl+1)*rho_int_expl(jl)/rho_int_old(jl)
        enddo
        
        do jl=1,n_layers-2
          ! lower diagonal
          c_vector(jl) = theta_v_int_new(jl+1)*gammaa(jl+1)*theta_v_int_new(jl) &
          + 0.5_wp*(grid%exner_bg(ji,jk,jl+1)-grid%exner_bg(ji,jk,jl+2)) &
          *(alpha(jl+1)+beta(jl+1)*theta_v_int_new(jl)) &
          - (grid%z_scalar(ji,jk,jl+1)-grid%z_scalar(ji,jk,jl+2))/(impl_weight*dtime*c_d_p)*0.5_wp &
          *state_old%wind_w(ji,jk,jl+2)/(grid%volume(ji,jk,jl+1)*rho_int_old(jl+1))
          ! upper diagonal
          e_vector(jl) = theta_v_int_new(jl)*gammaa(jl+1)*theta_v_int_new(jl+1) &
          - 0.5_wp*(grid%exner_bg(ji,jk,jl)-grid%exner_bg(ji,jk,jl+1)) &
          *(alpha(jl+1)+beta(jl+1)*theta_v_int_new(jl+1)) &
          + (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))/(impl_weight*dtime*c_d_p)*0.5_wp &
          *state_old%wind_w(ji,jk,jl+1)/(grid%volume(ji,jk,jl+1)*rho_int_old(jl))
        enddo
	
        ! soil components of the matrix
        if (soil_switch==1) then
          ! calculating the explicit part of the heat flux density
          do jl=1,nsoillays-1
            heat_flux_density_expl(jl) &
            = -grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk)*(state_old%temperature_soil(ji,jk,jl) &
            - state_old%temperature_soil(ji,jk,jl+1)) &
            /(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
          enddo
          
          heat_flux_density_expl(nsoillays) &
          = -grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk)*(state_old%temperature_soil(ji,jk,nsoillays) &
          - grid%t_const_soil(ji,jk)) &
          /(2._wp*(grid%z_soil_center(nsoillays) - grid%z_t_const))
          
          ! calculating the explicit part of the temperature change
          r_vector(n_layers) &
          ! old temperature
          = state_old%temperature_soil(ji,jk,1) &
          ! sensible heat flux
          + (diag%power_flux_density_sensible(ji,jk) &
          ! latent heat flux
          + diag%power_flux_density_latent(ji,jk) &
          ! shortwave inbound radiation
          + diag%sfc_sw_in(ji,jk) &
          ! longwave outbound radiation
          - diag%sfc_lw_out(ji,jk) &
          ! heat conduction from below
          + 0.5_wp*heat_flux_density_expl(1)) &
          /((grid%z_soil_interface(1) - grid%z_soil_interface(2))*grid%sfc_rho_c(ji,jk))*dtime

          ! loop over all soil layers below the first layer
          do jl=2,nsoillays
            r_vector(n_layers-1+jl) &
            ! old temperature
            = state_old%temperature_soil(ji,jk,jl) &
            ! heat conduction from above
            + 0.5_wp*(-heat_flux_density_expl(jl-1) &
            ! heat conduction from below
            + heat_flux_density_expl(jl)) &
            /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji,jk))*dtime
          enddo
  
          ! the diagonal component
          do jl=1,nsoillays
            if (jl==1) then
              d_vector(jl+n_layers-1) = 1.0_wp+0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
              /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji,jk)) &
              *1._wp/(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
            elseif (jl==nsoillays) then
              d_vector(jl+n_layers-1) = 1.0_wp+0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
              /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji,jk)) &
              *1._wp/(grid%z_soil_center(jl-1) - grid%z_soil_center(jl))
            else
              d_vector(jl+n_layers-1) = 1.0_wp+0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
              /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji,jk)) &
              *(1._wp/(grid%z_soil_center(jl-1) - grid%z_soil_center(jl)) &
              + 1._wp/(grid%z_soil_center(jl) - grid%z_soil_center(jl+1)))
            endif
          enddo
          
          ! the off-diagonal components
          c_vector(n_layers-1) = 0._wp
          e_vector(n_layers-1) = 0._wp
          do jl=1,nsoillays-1
            c_vector(jl+n_layers-1) = -0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
            /((grid%z_soil_interface(jl+1) - grid%z_soil_interface(jl+2))*grid%sfc_rho_c(ji,jk)) &
            /(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
            e_vector(jl+n_layers-1) = -0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
            /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji,jk)) &
            /(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
          enddo
        endif
		
        ! calling the subroutine to solve the system of linear equations
        call thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,n_layers-1+soil_switch*nsoillays)
       
        ! Klemp swamp layer
        do jl=1,n_layers-1
          above_damping = grid%z_w(ji,jk,jl+1)-damping_start_height
          if (above_damping<0._wp .or. .not. lklemp) then
            damping_coeff = 0._wp
          else
            damping_coeff = klemp_damp_max*sin(0.5_wp*M_PI*above_damping/(toa-damping_start_height))**2
          endif
          solution_vector(jl) = solution_vector(jl)/(1._wp+damping_coeff*dtime)
        enddo
        
        ! results
        ! density, virtual potential temperature density
        do jl=2,n_layers-1
          state_target%rho(ji,jk,jl,n_condensed_constituents+1) = rho_expl(jl) &
          + dtime*(-solution_vector(jl-1)+solution_vector(jl))/grid%volume(ji,jk,jl)
          state_target%rhotheta_v(ji,jk,jl) = rhotheta_v_expl(jl) &
          + dtime*(-theta_v_int_new(jl-1)*solution_vector(jl-1)+theta_v_int_new(jl)*solution_vector(jl))/grid%volume(ji,jk,jl)
        enddo
        ! uppermost layer
        state_target%rho(ji,jk,1,n_condensed_constituents+1) = rho_expl(1)+dtime*solution_vector(1)/grid%volume(ji,jk,1)
        state_target%rhotheta_v(ji,jk,1) = rhotheta_v_expl(1)+dtime*theta_v_int_new(1)*solution_vector(1)/grid%volume(ji,jk,1)
        ! lowest layer
        state_target%rho(ji,jk,n_layers,n_condensed_constituents+1) = rho_expl(n_layers) &
        - dtime*solution_vector(n_layers-1)/grid%volume(ji,jk,n_layers)
        state_target%rhotheta_v(ji,jk,n_layers) = rhotheta_v_expl(n_layers) &
        - dtime*theta_v_int_new(n_layers-1)*solution_vector(n_layers-1)/grid%volume(ji,jk,n_layers)
        ! vertical velocity
        do jl=2,n_layers
          rho_int_new = 0.5_wp*(state_target%rho(ji,jk,jl-1,n_condensed_constituents+1) &
          + state_target%rho(ji,jk,jl,n_condensed_constituents+1))
          state_target%wind_w(ji,jk,jl) = (2._wp*solution_vector(jl-1)/grid%area_z(ji,jk,jl) &
          - rho_int_new*state_old%wind_w(ji,jk,jl))/rho_int_old(jl-1)
        enddo
        ! Exner pressure
        do jl=1,n_layers
          state_target%exner_pert(ji,jk,jl) = state_old%exner_pert(ji,jk,jl) &
          + grid%volume(ji,jk,jl)*gammaa(jl)*(state_target%rhotheta_v(ji,jk,jl)-state_old%rhotheta_v(ji,jk,jl))
        enddo
        
        ! soil temperature
        if (soil_switch==1) then
          do jl=1,nsoillays
            state_target%temperature_soil(ji,jk,jl) = solution_vector(n_layers-1+jl)
          enddo
        endif
		
      enddo
    enddo
    !$omp end parallel do
    
    ! virtual potential temperature perturbation at the new time step
    !$omp parallel workshare
    state_target%theta_v_pert(:,:,:) = state_target%rhotheta_v(:,:,:)/state_target%rho(:,:,:,n_condensed_constituents+1) &
    - grid%theta_v_bg(:,:,:)
    !$omp end parallel workshare

  end subroutine three_band_solver_ver
  
  subroutine three_band_solver_gen_densities(state_old,state_new,tend,diag,grid,rk_step)
  
    ! Vertical advection of generalized densities (of tracers) with 3-band matrices.
    ! mass densities, density x temperatures
    
    ! input arguments and output
    type(t_state), intent(in)    :: state_old ! state at the old timestep
    type(t_state), intent(inout) :: state_new ! state at the new timestep
    type(t_tend),  intent(in)    :: tend      ! explicit tendencies
    type(t_diag),  intent(inout) :: diag      ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid      ! model grid
    integer,       intent(in)    :: rk_step   ! predictor-corrector substep index
    
    ! local variables
    integer  :: n_relevant_constituents                   ! number of relevant constituents for a certain quantity
    real(wp) :: impl_weight,expl_weight                   ! time stepping weights
    real(wp) :: c_vector(n_layers-1)                      ! vector for solving the system of linear equations
    real(wp) :: d_vector(n_layers)                        ! vector for solving the system of linear equations
    real(wp) :: e_vector(n_layers-1)                      ! vector for solving the system of linear equations
    real(wp) :: r_vector(n_layers)                        ! vector for solving the system of linear equations
    real(wp) :: vertical_flux_vector_impl(n_layers-1)     ! vertical flux at the new timestep
    real(wp) :: vertical_flux_vector_rhs(n_layers-1)      ! vertical flux at the old timestep
    real(wp) :: vertical_enthalpy_flux_vector(n_layers-1) ! vertical enthalpy flux density vector
    real(wp) :: solution_vector(n_layers)                 ! solution of the system of linear equations
    real(wp) :: density_old_at_interface,added_mass       ! abbreviations
    real(wp) :: temperature_old_at_interface              ! temperature in a level at the old PC substep
    integer  :: jc,ji,jk,jl                               ! loop indices
    
    ! setting the time stepping weights
    impl_weight = 0.5_wp
    expl_weight = 1._wp - impl_weight
    
    ! firstly the number of relevant constituents needs to be determined
    n_relevant_constituents = n_constituents ! the main gaseous constituent is excluded later

    ! loop over all relevant constituents
    do jc=1,n_relevant_constituents
      ! This is done do all tracers apart from the main gaseous constituent.
      if (jc/=n_condensed_constituents+1) then
        
        ! loop over all columns
        !$omp parallel do private(ji,jk,jl,vertical_flux_vector_impl,vertical_flux_vector_rhs,density_old_at_interface,c_vector, &
        !$omp d_vector,e_vector,r_vector,solution_vector,added_mass,vertical_enthalpy_flux_vector,temperature_old_at_interface)
        do ji=1,ny
          do jk=1,nx

            ! diagnozing the vertical fluxes
            do jl=1,n_layers-1
              ! resetting the vertical enthalpy flux density divergence
              if (rk_step==1 .and. jc==1) then
                diag%condensates_sediment_heat(ji,jk,jl) = 0._wp
              endif
              vertical_flux_vector_impl(jl) = state_old%wind_w(ji,jk,jl+1)
              vertical_flux_vector_rhs(jl) = state_new%wind_w(ji,jk,jl+1)
              
              ! for condensed constituents, a sink velocity must be added.
              ! precipitation
              ! snow
              if (jc<=n_condensed_constituents/4) then
                vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - snow_velocity
                vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - snow_velocity
              ! rain
              elseif (jc<=n_condensed_constituents/2) then
                vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - rain_velocity
                vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - rain_velocity
              ! clouds
              elseif (jc<=n_condensed_constituents-1) then
                vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - cloud_droplets_velocity
                vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - cloud_droplets_velocity
              ! graupel
              elseif (jc==n_condensed_constituents) then
                vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - graupel_velocity
                vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - graupel_velocity
              endif
              ! multiplying the vertical velocity by the area
              vertical_flux_vector_impl(jl) = grid%area_z(ji,jk,jl+1)*vertical_flux_vector_impl(jl)
              vertical_flux_vector_rhs(jl) = grid%area_z(ji,jk,jl+1)*vertical_flux_vector_rhs(jl)
              ! old density at the interface
              if (vertical_flux_vector_rhs(jl)>=0._wp) then
                density_old_at_interface = state_old%rho(ji,jk,jl+1,jc)
                temperature_old_at_interface = diag%temperature(ji,jk,jl+1)
              else
                density_old_at_interface = state_old%rho(ji,jk,jl,jc)
                temperature_old_at_interface = diag%temperature(ji,jk,jl+1)
              endif
              vertical_flux_vector_rhs(jl) = density_old_at_interface*vertical_flux_vector_rhs(jl)
              vertical_enthalpy_flux_vector(jl) = c_p_cond(jc,temperature_old_at_interface) &
                                                  *temperature_old_at_interface*vertical_flux_vector_rhs(jl)
            enddo
            if (rk_step==1 .and. jc==1) then
              diag%condensates_sediment_heat(ji,jk,n_layers) = 0._wp
            endif

            ! Now we proceed to solving the vertical tridiagonal problems.

            ! filling up the original vectors
            do jl=1,n_layers-1
              if (vertical_flux_vector_impl(jl)>=0._wp) then
                c_vector(jl) = 0._wp
                e_vector(jl) = -impl_weight*dtime/grid%volume(ji,jk,jl)*vertical_flux_vector_impl(jl)
              else
                c_vector(jl) = impl_weight*dtime/grid%volume(ji,jk,jl+1)*vertical_flux_vector_impl(jl)
                e_vector(jl) = 0._wp
              endif
            enddo
            do jl=1,n_layers
              if (jl==1) then
                if (vertical_flux_vector_impl(1)>=0._wp) then
                  d_vector(jl) = 1._wp
                else
                  d_vector(jl) = 1._wp - impl_weight*dtime/grid%volume(ji,jk,jl)*vertical_flux_vector_impl(1)
                endif
              elseif (jl==n_layers) then
                if (vertical_flux_vector_impl(jl-1)>=0._wp) then
                  d_vector(jl) = 1._wp + impl_weight*dtime/grid%volume(ji,jk,jl)*vertical_flux_vector_impl(jl-1)
                else
                  d_vector(jl) = 1._wp
                endif
                ! precipitation
                ! snow
                if (jc<=n_condensed_constituents/4) then
                  d_vector(jl) = d_vector(jl) + impl_weight*snow_velocity*dtime &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                ! rain
                elseif (jc<=n_condensed_constituents/2) then
                  d_vector(jl) = d_vector(jl) + impl_weight*rain_velocity*dtime &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                ! clouds
                elseif (jc<=n_condensed_constituents) then
                  d_vector(jl) = d_vector(jl) + impl_weight*cloud_droplets_velocity*dtime &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                endif
              else
                d_vector(jl) = 1._wp
                if (vertical_flux_vector_impl(jl-1)>=0._wp) then
                  d_vector(jl) = d_vector(jl) + impl_weight*dtime/grid%volume(ji,jk,jl)*vertical_flux_vector_impl(jl-1)
                endif
                if (vertical_flux_vector_impl(jl)<0._wp) then
                  d_vector(jl) = d_vector(jl) - impl_weight*dtime/grid%volume(ji,jk,jl)*vertical_flux_vector_impl(jl)
                endif
              endif
              ! the explicit component
              r_vector(jl) = state_old%rho(ji,jk,jl,jc) + dtime*tend%rho(ji,jk,jl,jc)
              ! adding the explicit part of the vertical flux divergence
              if (jl==1) then
                r_vector(jl) = r_vector(jl) + expl_weight*dtime*vertical_flux_vector_rhs(jl)/grid%volume(ji,jk,jl)
                if (rk_step==1 .and. jc<=n_condensed_constituents) then
                  diag%condensates_sediment_heat(ji,jk,jl) = diag%condensates_sediment_heat(ji,jk,jl) &
                  + vertical_enthalpy_flux_vector(jl)/grid%volume(ji,jk,jl)
                endif
              elseif (jl==n_layers) then
                r_vector(jl) = r_vector(jl) - expl_weight*dtime*vertical_flux_vector_rhs(jl-1)/grid%volume(ji,jk,jl)
                if (rk_step==1 .and. jc<=n_condensed_constituents) then
                  diag%condensates_sediment_heat(ji,jk,jl) = diag%condensates_sediment_heat(ji,jk,jl) &
                  - vertical_enthalpy_flux_vector(jl-1)/grid%volume(ji,jk,jl)
                endif
                ! precipitation
                ! snow
                if (jc<=n_condensed_constituents/4) then
                  r_vector(jl) = r_vector(jl) - expl_weight*snow_velocity*dtime*state_old%rho(ji,jk,jl,jc) &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                ! rain
                elseif (jc<=n_condensed_constituents/2) then
                  r_vector(jl) = r_vector(jl) - expl_weight*rain_velocity*dtime*state_old%rho(ji,jk,jl,jc) &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                  if (rk_step==1) then
                    diag%condensates_sediment_heat(ji,jk,jl) = diag%condensates_sediment_heat(ji,jk,jl) &
                    - snow_velocity &
                    *diag%temperature(ji,jk,n_layers)*c_p_cond(jc,diag%temperature(ji,jk,n_layers)) &
                    *state_old%rho(ji,jk,n_layers,jc)*grid%area_z(ji,jk,n_levels)/grid%volume(ji,jk,jl)
                  endif
                ! clouds
                elseif (jc<=n_condensed_constituents) then
                  r_vector(jl) = r_vector(jl) - expl_weight*cloud_droplets_velocity*dtime*state_old%rho(ji,jk,jl,jc) &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                  if (rk_step==1) then
                    diag%condensates_sediment_heat(ji,jk,jl) = diag%condensates_sediment_heat(ji,jk,jl) &
                    -cloud_droplets_velocity &
                    *diag%temperature(ji,jk,n_layers)*c_p_cond(jc,diag%temperature(ji,jk,n_layers)) &
                    *state_old%rho(ji,jk,n_layers,jc)*grid%area_z(ji,jk,n_levels)/grid%volume(ji,jk,jl)
                  endif
                endif
              else
                r_vector(jl) = r_vector(jl) + expl_weight*dtime*(-vertical_flux_vector_rhs(jl-1)+vertical_flux_vector_rhs(jl)) &
                               /grid%volume(ji,jk,jl)
                if (rk_step==1 .and. jc<=n_condensed_constituents) then
                  diag%condensates_sediment_heat(ji,jk,jl) = diag%condensates_sediment_heat(ji,jk,jl) &
                  + (-vertical_enthalpy_flux_vector(jl-1) + vertical_enthalpy_flux_vector(jl))/grid%volume(ji,jk,jl)
                endif
              endif
            enddo

            ! calling the algorithm to solve the system of linear equations
            call thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,n_layers)
      
            ! this should account for round-off errors only
            do jl=1,n_layers
              if (solution_vector(jl)<0._wp) then
                solution_vector(jl) = 0._wp
              endif
            enddo

            ! writing the result into the new state
            do jl=1,n_layers
              state_new%rho(ji,jk,jl,jc) = solution_vector(jl)
            enddo
          
          enddo ! column index
        enddo ! line index
        !$omp end parallel do
      endif
    enddo ! constituent
  
  end subroutine three_band_solver_gen_densities
  
  subroutine thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,solution_length)

    ! This subroutine solves a system of linear equations with a three-band matrix.
    
    real(wp), intent(in)    :: c_vector(:)        ! lower diagonal vector
    real(wp), intent(in)    :: d_vector(:)        ! main diagonal vector
    real(wp), intent(in)    :: e_vector(:)        ! upper diagonal vector
    real(wp), intent(in)    :: r_vector(:)        ! right hand side vector
    real(wp), intent(inout) :: solution_vector(:) ! vector containing the solution
    integer,  intent(in)    :: solution_length    ! length of the solution vector
    
    ! local variables
    real(wp) :: e_prime_vector(solution_length-1) ! help vector for solving the matrix equation
    real(wp) :: r_prime_vector(solution_length)   ! help vector for solving the matrix equation
    integer  :: jl                                ! loop index
    
    ! downward sweep (matrix)
    e_prime_vector(1) = e_vector(1)/d_vector(1)
    do jl=2,solution_length-1
      e_prime_vector(jl) = e_vector(jl)/(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
    enddo
    ! downward sweep (right-hand side)
    r_prime_vector(1) = r_vector(1)/d_vector(1)
    do jl=2,solution_length
      r_prime_vector(jl) = (r_vector(jl) - r_prime_vector(jl-1)*c_vector(jl-1)) &
      /(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
    enddo
    
    ! upward sweep (final solution)
    solution_vector(solution_length) = r_prime_vector(solution_length)
    do jl=solution_length-1,1,-1
      solution_vector(jl) = r_prime_vector(jl) - e_prime_vector(jl)*solution_vector(jl+1)
    enddo
  
  end subroutine thomas_algorithm

end module mo_column_solvers











