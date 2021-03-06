! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module column_solvers

  ! This module contains the implicit vertical routines (implicit part of the HEVI scheme).

  use run_nml,          only: nlins,ncols,wp,nlays,dtime,toa,impl_weight,partial_impl_weight
  use constituents_nml, only: no_of_condensed_constituents,no_of_constituents, &
                              snow_velocity,rain_velocity,cloud_droplets_velocity
  use definitions,      only: t_grid,t_state,t_tend,t_diag
  use dictionary,       only: spec_heat_capacities_v_gas,spec_heat_capacities_p_gas,specific_gas_constants
  use diff_nml,         only: lklemp,klemp_damp_max,klemp_begin_rel
  use surface_nml,      only: nsoillays,lprog_soil_temp,lsfc_sensible_heat_flux
  use constants,        only: M_PI

  implicit none
  
  private
  
  public :: three_band_solver_ver
  public :: three_band_solver_gen_densities
  
  contains
  
  subroutine three_band_solver_ver(state_old,state_new,diag,tend,grid,rk_step)
  
    ! This subroutine is the main implicit vertical solver.

    ! input arguments and output
    type(t_state), intent(in)    :: state_old ! state at the old timestep
    type(t_state), intent(inout) :: state_new ! state at the new timestep
    type(t_diag),  intent(inout) :: diag      ! diagnostic quantities
    type(t_tend),  intent(inout) :: tend      ! explicit tendencies
    type(t_grid),  intent(in)    :: grid      ! model grid
    integer,       intent(in)    :: rk_step   ! Runge Kutta substep

    ! local variables
    integer  :: soil_switch                        ! soil switch: 0 if soil does not have to be calculated off, 1 if soil has to be calculated
    real(wp) :: c_vector(nlays-2+nsoillays)        ! needed for the vertical solver
    real(wp) :: d_vector(nlays-1+nsoillays)        ! needed for the vertical solver
    real(wp) :: e_vector(nlays-2+nsoillays)        ! needed for the vertical solver
    real(wp) :: r_vector(nlays-1+nsoillays)        ! needed for the vertical solver
    real(wp) :: rho_expl(nlays)                    ! explicit mass density
    real(wp) :: rhotheta_v_expl(nlays)             ! explicit virtual potential temperature density
    real(wp) :: exner_pert_expl(nlays)             ! explicit Exner pressure perturbation
    real(wp) :: theta_v_pert_expl(nlays)           ! explicit virtual potential temperature perturbation
    real(wp) :: rho_int_old(nlays-1)               ! old interface mass density
    real(wp) :: rho_int_expl(nlays-1)              ! explicit interface mass density
    real(wp) :: theta_v_int_new(nlays-1)           ! preliminary new virtual potential temperature interface values
    real(wp) :: rho_int_new                        ! new density interface value
    real(wp) :: alpha_old(nlays)                   ! alpha at the old time step
    real(wp) :: beta_old(nlays)                    ! beta at the old time step
    real(wp) :: gamma_old(nlays)                   ! gamma at the old time step
    real(wp) :: alpha_new(nlays)                   ! alpha at the new time step
    real(wp) :: beta_new(nlays)                    ! beta at the new time step
    real(wp) :: gamma_new(nlays)                   ! gamma at the new time step
    real(wp) :: alpha(nlays)                       ! alpha
    real(wp) :: beta(nlays)                        ! beta
    real(wp) :: gammaa(nlays)                      ! gamma
    real(wp) :: c_v                                ! specific heat capacity at constant volume
    real(wp) :: c_p                                ! specific heat capacity at constant pressure
    real(wp) :: r_d                                ! individual gas constant of dry air
    real(wp) :: damping_start_height               ! lower boundary height of the Klemp layer
    real(wp) :: damping_coeff                      ! damping coefficient of the Klemp layer
    real(wp) :: above_damping                      ! height above the lower boundary of the damping height
    real(wp) :: t_gas_lowest_layer_old             ! temperature of the gas in the lowest layer of the model atmosphere, old timestep
    real(wp) :: t_gas_lowest_layer_new             ! temperature of the gas in the lowest layer of the model atmosphere, new timestep
    real(wp) :: heat_flux_density_expl(nsoillays)  ! explicit heat_flux_density in the soil
    real(wp) :: solution_vector(nlays-1+nsoillays) ! vector containing the solution of the linear problem to solve here
    integer  :: ji,jk,jl                           ! loop variables

    c_v = spec_heat_capacities_v_gas(0)
    c_p = spec_heat_capacities_p_gas(0)
    r_d = specific_gas_constants(0)
    damping_start_height = klemp_begin_rel*toa
    
    ! calculating the sensible power flux density if soil is switched on
    if (lsfc_sensible_heat_flux) then
      !$omp parallel do private(ji,jk,t_gas_lowest_layer_old,t_gas_lowest_layer_new)
      do ji=1,nlins
        do jk=1,ncols

          ! gas temperature in the lowest layer
          t_gas_lowest_layer_old = (grid%exner_bg(ji,jk,nlays)+state_old%exner_pert(ji,jk,nlays)) &
          *(grid%theta_v_bg(ji,jk,nlays)+state_old%theta_v_pert(ji,jk,nlays))
          t_gas_lowest_layer_new = (grid%exner_bg(ji,jk,nlays)+state_new%exner_pert(ji,jk,nlays)) &
          *(grid%theta_v_bg(ji,jk,nlays)+state_new%theta_v_pert(ji,jk,nlays))

          ! the sensible power flux density
          diag%power_flux_density_sensible(ji,jk) = 0.5_wp*c_v*(state_new%rho(ji,jk,nlays,no_of_condensed_constituents+1) &
          *(t_gas_lowest_layer_old - state_old%temperature_soil(ji,jk,1)) &
          + state_old%rho(ji,jk,nlays,no_of_condensed_constituents+1) &
          *(t_gas_lowest_layer_new - state_new%temperature_soil(ji,jk,1)))/diag%scalar_flux_resistance(ji,jk)

          ! contribution of sensible heat to rhotheta_v
          tend%rhotheta_v(ji,jk,nlays) = tend%rhotheta_v(ji,jk,nlays) &
          -grid%area_z(ji,jk,nlays+1)*diag%power_flux_density_sensible(ji,jk) &
          /((grid%exner_bg(ji,jk,nlays)+state_new%exner_pert(ji,jk,nlays))*c_p)/grid%volume(ji,jk,nlays)
          
        enddo
      enddo
      !$omp end parallel do
    endif
	
    !$omp parallel do private(ji,jk,jl,c_vector,d_vector,e_vector,r_vector,solution_vector, &
    !$omp rho_expl,rhotheta_v_expl,exner_pert_expl,theta_v_pert_expl,rho_int_old, &
    !$omp rho_int_expl,theta_v_int_new,rho_int_new,alpha_old,beta_old,gamma_old,alpha_new, &
    !$omp beta_new,gamma_new,alpha,beta,gammaa,damping_coeff,above_damping,soil_switch)
    do ji=1,nlins
      do jk=1,ncols
      
        ! determining wether soil needs to be calculated
        soil_switch = 0
        if (lprog_soil_temp .and. grid%is_land(ji,jk)==1) then
          soil_switch=1
        endif
      
        ! explicit quantities
        do jl=1,nlays
          ! explicit density
          rho_expl(jl) = state_old%rho(ji,jk,jl,no_of_condensed_constituents+1) &
          + dtime*tend%rho(ji,jk,jl,no_of_condensed_constituents+1)
          ! explicit virtual potential temperature density
          rhotheta_v_expl(jl) = state_old%rhotheta_v(ji,jk,jl) + dtime*tend%rhotheta_v(ji,jk,jl)
          if (rk_step==1) then
            ! old time step partial derivatives of theta_v and Pi (divided by the volume)
            alpha(jl) = -state_old%rhotheta_v(ji,jk,jl)/state_old%rho(ji,jk,jl,no_of_condensed_constituents+1)**2 &
            /grid%volume(ji,jk,jl)
            beta(jl)  = 1._wp/state_old%rho(ji,jk,jl,no_of_condensed_constituents+1)/grid%volume(ji,jk,jl)
            gammaa(jl) = r_d/(c_v*state_old%rhotheta_v(ji,jk,jl))* &
            (grid%exner_bg(ji,jk,jl)+state_old%exner_pert(ji,jk,jl))/grid%volume(ji,jk,jl)
          else
            ! old time step partial derivatives of theta_v and Pi
            alpha_old(jl) = -state_old%rhotheta_v(ji,jk,jl)/state_old%rho(ji,jk,jl,no_of_condensed_constituents+1)**2
            beta_old(jl)  = 1._wp/state_old%rho(ji,jk,jl,no_of_condensed_constituents+1)
            gamma_old(jl) = r_d/(c_v*state_old%rhotheta_v(ji,jk,jl))* &
            (grid%exner_bg(ji,jk,jl)+state_old%exner_pert(ji,jk,jl))
            ! new time step partial derivatives of theta_v and Pi
            alpha_new(jl) = -state_new%rhotheta_v(ji,jk,jl)/state_new%rho(ji,jk,jl,no_of_condensed_constituents+1)**2
            beta_new(jl)  = 1._wp/state_new%rho(ji,jk,jl,no_of_condensed_constituents+1)
            gamma_new(jl) = r_d/(c_v*state_new%rhotheta_v(ji,jk,jl)) &
            *(grid%exner_bg(ji,jk,jl)+state_new%exner_pert(ji,jk,jl))
            ! interpolation of partial derivatives of theta_v and Pi (divided by the volume)
            alpha(jl) = ((1._wp - partial_impl_weight)*alpha_old(jl)+partial_impl_weight*alpha_new(jl))/grid%volume(ji,jk,jl)
            beta(jl) = ((1._wp - partial_impl_weight)*beta_old (jl)+partial_impl_weight*beta_new(jl))/grid%volume(ji,jk,jl)
            gammaa(jl) = ((1._wp - partial_impl_weight)*gamma_old(jl)+partial_impl_weight*gamma_new(jl))/grid%volume(ji,jk,jl)
          endif
          ! explicit virtual potential temperature perturbation
          theta_v_pert_expl(jl) = state_old%theta_v_pert(ji,jk,jl) &
          + dtime*grid%volume(ji,jk,jl)*(alpha(jl)*tend%rho(ji,jk,jl,no_of_condensed_constituents+1) &
          + beta(jl)*tend%rhotheta_v(ji,jk,jl))
          ! explicit Exner pressure perturbation
          exner_pert_expl(jl) = state_old%exner_pert(ji,jk,jl)+dtime*grid%volume(ji,jk,jl)*gammaa(jl)*tend%rhotheta_v(ji,jk,jl)
        enddo
        
        ! interface values
        do jl=1,nlays-1
          rho_int_old(jl) = 0.5_wp*(state_old%rho(ji,jk,jl,no_of_condensed_constituents+1) &
          + state_old%rho(ji,jk,jl+1,no_of_condensed_constituents+1))
          rho_int_expl(jl) = 0.5_wp*(rho_expl(jl)+rho_expl(jl+1))
          theta_v_int_new(jl) = 0.5_wp*(state_new%rhotheta_v(ji,jk,jl)/state_new%rho(ji,jk,jl,no_of_condensed_constituents+1) &
          + state_new%rhotheta_v(ji,jk,jl+1)/state_new%rho(ji,jk,jl+1,no_of_condensed_constituents+1))
        enddo
      
        ! filling up the coefficient vectors
        do jl=1,nlays-1
          ! main diagonal
          d_vector(jl) = -theta_v_int_new(jl)**2*(gammaa(jl)+gammaa(jl+1)) &
          + 0.5_wp*(grid%exner_bg(ji,jk,jl)-grid%exner_bg(ji,jk,jl+1)) &
          *(alpha(jl+1)-alpha(jl)+theta_v_int_new(jl)*(beta(jl+1)-beta(jl))) &
          - (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))/(impl_weight*dtime**2*c_p*rho_int_old(jl)) &
          *(2._wp/grid%area_z(ji,jk,jl+1)+dtime*state_old%wind_w(ji,jk,jl+1)*0.5_wp &
          *(-1._wp/grid%volume(ji,jk,jl)+1._wp/grid%volume(ji,jk,jl+1)))
          ! right hand side
          r_vector(jl) = -(state_old%wind_w(ji,jk,jl+1)+dtime*tend%wind_w(ji,jk,jl+1))* &
          (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1)) &
          /(impl_weight*dtime**2*c_p) &
          + theta_v_int_new(jl)*(exner_pert_expl(jl)-exner_pert_expl(jl+1))/dtime &
          + 0.5_wp/dtime*(theta_v_pert_expl(jl)+theta_v_pert_expl(jl+1))*(grid%exner_bg(ji,jk,jl)-grid%exner_bg(ji,jk,jl+1)) &
          - (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))/(impl_weight*dtime**2*c_p) &
          *state_old%wind_w(ji,jk,jl+1)*rho_int_expl(jl)/rho_int_old(jl)
        enddo
        
        do jl=1,nlays-2
          ! lower diagonal
          c_vector(jl) = theta_v_int_new(jl+1)*gammaa(jl+1)*theta_v_int_new(jl) &
          + 0.5_wp*(grid%exner_bg(ji,jk,jl+1)-grid%exner_bg(ji,jk,jl+2)) &
          *(alpha(jl+1)+beta(jl+1)*theta_v_int_new(jl)) &
          - (grid%z_scalar(ji,jk,jl+1)-grid%z_scalar(ji,jk,jl+2))/(impl_weight*dtime*c_p)*0.5_wp &
          *state_old%wind_w(ji,jk,jl+2)/(grid%volume(ji,jk,jl+1)*rho_int_old(jl+1))
          ! upper diagonal
          e_vector(jl) = theta_v_int_new(jl)*gammaa(jl+1)*theta_v_int_new(jl+1) &
          - 0.5_wp*(grid%exner_bg(ji,jk,jl)-grid%exner_bg(ji,jk,jl+1)) &
          *(alpha(jl+1)+beta(jl+1)*theta_v_int_new(jl+1)) &
          + (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))/(impl_weight*dtime*c_p)*0.5_wp &
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
          r_vector(nlays) &
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
            r_vector(nlays-1+jl) &
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
              d_vector(jl+nlays-1) = 1.0_wp+0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
              /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji,jk)) &
              *1._wp/(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
            elseif (jl==nsoillays) then
              d_vector(jl+nlays-1) = 1.0_wp+0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
              /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji,jk)) &
              *1._wp/(grid%z_soil_center(jl-1) - grid%z_soil_center(jl))
            else
              d_vector(jl+nlays-1) = 1.0_wp+0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
              /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji,jk)) &
              *(1._wp/(grid%z_soil_center(jl-1) - grid%z_soil_center(jl)) &
              + 1._wp/(grid%z_soil_center(jl) - grid%z_soil_center(jl+1)))
            endif
          enddo
          
          ! the off-diagonal components
          c_vector(nlays-1) = 0._wp
          e_vector(nlays-1) = 0._wp
          do jl=1,nsoillays-1
            c_vector(jl+nlays-1) = -0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
            /((grid%z_soil_interface(jl+1) - grid%z_soil_interface(jl+2))*grid%sfc_rho_c(ji,jk)) &
            /(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
            e_vector(jl+nlays-1) = -0.5_wp*dtime*grid%sfc_rho_c(ji,jk)*grid%t_conduc_soil(ji,jk) &
            /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji,jk)) &
            /(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
          enddo
        endif
		
		! calling the subroutine to solve the system of linear equations
        call thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,nlays-1+soil_switch*nsoillays)
       
        ! Klemp swamp layer
        do jl=1,nlays-1
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
        do jl=2,nlays-1
          state_new%rho(ji,jk,jl,no_of_condensed_constituents+1) = rho_expl(jl) &
          + dtime*(-solution_vector(jl-1)+solution_vector(jl))/grid%volume(ji,jk,jl)
          state_new%rhotheta_v(ji,jk,jl) = rhotheta_v_expl(jl) &
          + dtime*(-theta_v_int_new(jl-1)*solution_vector(jl-1)+theta_v_int_new(jl)*solution_vector(jl))/grid%volume(ji,jk,jl)
        enddo
        ! uppermost layer
        state_new%rho(ji,jk,1,no_of_condensed_constituents+1) = rho_expl(1)+dtime*solution_vector(1)/grid%volume(ji,jk,1)
        state_new%rhotheta_v(ji,jk,1) = rhotheta_v_expl(1)+dtime*theta_v_int_new(1)*solution_vector(1)/grid%volume(ji,jk,1)
        ! lowest layer
        state_new%rho(ji,jk,nlays,no_of_condensed_constituents+1) = rho_expl(nlays) &
        - dtime*solution_vector(nlays-1)/grid%volume(ji,jk,nlays)
        state_new%rhotheta_v(ji,jk,nlays) = rhotheta_v_expl(nlays) &
        - dtime*theta_v_int_new(nlays-1)*solution_vector(nlays-1)/grid%volume(ji,jk,nlays)
        ! vertical velocity
        do jl=2,nlays
          rho_int_new = 0.5_wp*(state_new%rho(ji,jk,jl-1,no_of_condensed_constituents+1) &
          + state_new%rho(ji,jk,jl,no_of_condensed_constituents+1))
          state_new%wind_w(ji,jk,jl) = (2._wp*solution_vector(jl-1)/grid%area_z(ji,jk,jl) &
          - rho_int_new*state_old%wind_w(ji,jk,jl))/rho_int_old(jl-1)
        enddo
        ! Exner pressure
        do jl=1,nlays
          state_new%exner_pert(ji,jk,jl) = state_old%exner_pert(ji,jk,jl) &
          + grid%volume(ji,jk,jl)*gammaa(jl)*(state_new%rhotheta_v(ji,jk,jl)-state_old%rhotheta_v(ji,jk,jl))
        enddo
        
        ! soil temperature
        if (soil_switch==1) then
          do jl=1,nsoillays
            state_new%temperature_soil(ji,jk,jl) = solution_vector(nlays-1+jl)
          enddo
        endif
		
      enddo
    enddo
    !$omp end parallel do
    
    ! virtual potential temperature perturbation at the new time step
    !$omp parallel workshare
    state_new%theta_v_pert(:,:,:) = state_new%rhotheta_v(:,:,:)/state_new%rho(:,:,:,no_of_condensed_constituents+1) &
    - grid%theta_v_bg(:,:,:)
    !$omp end parallel workshare

  end subroutine three_band_solver_ver
  
  subroutine three_band_solver_gen_densities(state_old,state_new,tend,grid)
  
    ! Vertical advection of generalized densities (of tracers) with 3-band matrices.
    ! mass densities, density x temperatures
    
    ! input arguments and output
    type(t_state), intent(in)    :: state_old ! state at the old timestep
    type(t_state), intent(inout) :: state_new ! state at the new timestep
    type(t_tend),  intent(in)    :: tend      ! explicit tendencies
    type(t_grid),  intent(in)    :: grid      ! model grid
    
    ! local variables
    integer  :: no_of_relevant_constituents         ! number of relevant constituents for a certain quantity
    real(wp) :: impl_weight,expl_weight             ! time stepping weights
    real(wp) :: c_vector(nlays-1)                   ! vector for solving the system of linear equations
    real(wp) :: d_vector(nlays)                     ! vector for solving the system of linear equations
    real(wp) :: e_vector(nlays-1)                   ! vector for solving the system of linear equations
    real(wp) :: r_vector(nlays)                     ! vector for solving the system of linear equations
    real(wp) :: vertical_flux_vector_impl(nlays-1)  ! vertical flux at the new timestep
    real(wp) :: vertical_flux_vector_rhs(nlays-1)   ! vertical flux at the old timestep
    real(wp) :: solution_vector(nlays)              ! solution of the system of linear equations
    real(wp) :: density_old_at_interface,added_mass ! abbreviations
    integer  :: j_constituent,ji,jk,jl              ! loop indices
    
    ! setting the time stepping weights
    impl_weight = 0.5_wp
    expl_weight = 1._wp - impl_weight
    
    ! firstly the number of relevant constituents needs to be determined
    no_of_relevant_constituents = 0
    ! mass densities
    no_of_relevant_constituents = no_of_constituents ! the main gaseous constituent is excluded later

    ! loop over all relevant constituents
    do j_constituent=1,no_of_relevant_constituents
      ! This is done do all tracers apart from the main gaseous constituent.
      if (j_constituent/=no_of_condensed_constituents+1) then
        
        ! loop over all columns
        !$omp parallel do private(ji,jk,jl,vertical_flux_vector_impl,vertical_flux_vector_rhs,density_old_at_interface,c_vector, &
        !$omp d_vector,e_vector,r_vector,solution_vector,added_mass)
        do ji=1,nlins
          do jk=1,ncols

            ! diagnozing the vertical fluxes
            do jl=1,nlays-1
              vertical_flux_vector_impl(jl) = state_old%wind_w(ji,jk,jl+1)
              vertical_flux_vector_rhs(jl) = state_new%wind_w(ji,jk,jl+1)
              
              ! for condensed constituents, a sink velocity must be added.
              ! precipitation
              ! snow
              if (j_constituent<=no_of_condensed_constituents/4) then
                vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - snow_velocity
                vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - snow_velocity
              ! rain
              elseif (j_constituent<=no_of_condensed_constituents/2) then
                vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - rain_velocity
                vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - rain_velocity
              ! clouds
              elseif (j_constituent<=no_of_condensed_constituents) then
                vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - cloud_droplets_velocity
                vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - cloud_droplets_velocity
              endif
              ! multiplying the vertical velocity by the area
              vertical_flux_vector_impl(jl) = grid%area_z(ji,jk,jl+1)*vertical_flux_vector_impl(jl)
              vertical_flux_vector_rhs(jl) = grid%area_z(ji,jk,jl+1)*vertical_flux_vector_rhs(jl)
              ! old density at the interface
              if (vertical_flux_vector_rhs(jl)>=0._wp) then
                density_old_at_interface = state_old%rho(ji,jk,jl+1,j_constituent)
              else
                density_old_at_interface = state_old%rho(ji,jk,jl,j_constituent)
              endif
              vertical_flux_vector_rhs(jl) = density_old_at_interface*vertical_flux_vector_rhs(jl)
            enddo

            ! Now we proceed to solving the vertical tridiagonal problems.

            ! filling up the original vectors
            do jl=1,nlays-1
              if (vertical_flux_vector_impl(jl)>=0._wp) then
                c_vector(jl) = 0._wp
                e_vector(jl) = -impl_weight*dtime/grid%volume(ji,jk,jl)*vertical_flux_vector_impl(jl)
              else
                c_vector(jl) = impl_weight*dtime/grid%volume(ji,jk,jl+1)*vertical_flux_vector_impl(jl)
                e_vector(jl) = 0._wp
              endif
            enddo
            do jl=1,nlays
              if (jl==1) then
                if (vertical_flux_vector_impl(1)>=0._wp) then
                  d_vector(jl) = 1._wp
                else
                  d_vector(jl) = 1._wp - impl_weight*dtime/grid%volume(ji,jk,jl)*vertical_flux_vector_impl(1)
                endif
              elseif (jl==nlays) then
                if (vertical_flux_vector_impl(jl-1)>=0._wp) then
                  d_vector(jl) = 1._wp + impl_weight*dtime/grid%volume(ji,jk,jl)*vertical_flux_vector_impl(jl-1)
                else
                  d_vector(jl) = 1._wp
                endif
                ! precipitation
                ! snow
                if (j_constituent<=no_of_condensed_constituents/4) then
                  d_vector(jl) = d_vector(jl) + impl_weight*snow_velocity*dtime &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                ! rain
                elseif (j_constituent<=no_of_condensed_constituents/2) then
                  d_vector(jl) = d_vector(jl) + impl_weight*rain_velocity*dtime &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                ! clouds
                elseif (j_constituent<=no_of_condensed_constituents) then
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
              r_vector(jl) = state_old%rho(ji,jk,jl,j_constituent) + dtime*tend%rho(ji,jk,jl,j_constituent)
              ! adding the explicit part of the vertical flux divergence
              if (jl==1) then
                r_vector(jl) = r_vector(jl) + expl_weight*dtime*vertical_flux_vector_rhs(jl)/grid%volume(ji,jk,jl)
              elseif (jl==nlays) then
                r_vector(jl) = r_vector(jl) - expl_weight*dtime*vertical_flux_vector_rhs(jl-1)/grid%volume(ji,jk,jl)
                ! precipitation
                ! snow
                if (j_constituent<=no_of_condensed_constituents/4) then
                  r_vector(jl) = r_vector(jl) - expl_weight*snow_velocity*dtime*state_old%rho(ji,jk,jl,j_constituent) &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                ! rain
                elseif (j_constituent<=no_of_condensed_constituents/2) then
                  r_vector(jl) = r_vector(jl) - expl_weight*rain_velocity*dtime*state_old%rho(ji,jk,jl,j_constituent) &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                ! clouds
                elseif (j_constituent<=no_of_condensed_constituents) then
                  r_vector(jl) = r_vector(jl) - expl_weight*cloud_droplets_velocity*dtime*state_old%rho(ji,jk,jl,j_constituent) &
                  *grid%area_z(ji,jk,jl+1)/grid%volume(ji,jk,jl)
                endif
              else
                r_vector(jl) = r_vector(jl) + expl_weight*dtime*(-vertical_flux_vector_rhs(jl-1)+vertical_flux_vector_rhs(jl)) &
                /grid%volume(ji,jk,jl)
              endif
            enddo

            ! calling the algorithm to solve the system of linear equations
            call thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,nlays)
      
            ! this should account for round-off errors only
            do jl=1,nlays
              if (solution_vector(jl)<0._wp) then
                solution_vector(jl) = 0._wp
              endif
            enddo

            ! writing the result into the new state
            do jl=1,nlays
              state_new%rho(ji,jk,jl,j_constituent) = solution_vector(jl)
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

end module column_solvers











