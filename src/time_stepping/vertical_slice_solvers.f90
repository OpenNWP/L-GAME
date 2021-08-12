! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module vertical_slice_solvers

  ! This module contains the implicit vertical routines (implicit part of the HEVI scheme).

  use run_nml,        only: nlins,ncols,wp,nlays,dtime,p_0,toa
  use definitions,    only: t_grid,t_state,t_tend
  use thermodynamics, only: spec_heat_cap_diagnostics_v,spec_heat_cap_diagnostics_p,gas_constant_diagnostics
  use diff_nml,       only: lklemp,klemp_damp_max,klemp_begin_rel

  implicit none
  
  private
  
  public :: three_band_solver_ver
  
  contains
  
  subroutine three_band_solver_ver(state_old,state_new,tend,grid,rk_step)

    type(t_state), intent(in)    :: state_old ! state at the old timestep
    type(t_state), intent(inout) :: state_new ! state at the new timestep
    type(t_tend),  intent(in)    :: tend      ! explicit tendencies
    type(t_grid),  intent(in)    :: grid      ! model grid
    integer,       intent(in)    :: rk_step   ! Runge Kutta substep

    ! local variables
    real(wp)                 :: c_vector(nlays-2)       ! needed for the vertical solver
    real(wp)                 :: d_vector(nlays-1)       ! needed for the vertical solver
    real(wp)                 :: e_vector(nlays-2)       ! needed for the vertical solver
    real(wp)                 :: r_vector(nlays-1)       ! needed for the vertical solver
    real(wp)                 :: solution(nlays-1)       ! covariant mass flux density at the interfaces (solution)
    real(wp)                 :: rho_expl(nlays)         ! explicit mass density
    real(wp)                 :: rhotheta_expl(nlays)    ! explicit potential temperature density
    real(wp)                 :: exner_pert_expl(nlays)  ! explicit Exner pressure perturbation
    real(wp)                 :: theta_pert_expl(nlays)  ! explicit potential temperature perturbation
    real(wp)                 :: rho_int_old(nlays-1)    ! old interface mass density
    real(wp)                 :: rho_int_expl(nlays-1)   ! explicit interface mass density
    real(wp)                 :: theta_int_new(nlays-1)  ! preliminary new potential temperature interface values
    integer                  :: ji,jk,jl                ! loop variables
    real(wp)                 :: rho_int_new             ! new density interface value
    real(wp)                 :: alpha_old(nlays)        ! alpha at the old time step
    real(wp)                 :: beta_old(nlays)         ! beta at the old time step
    real(wp)                 :: gamma_old(nlays)        ! gamma at the old time step
    real(wp)                 :: alpha_new(nlays)        ! alpha at the new time step
    real(wp)                 :: beta_new(nlays)         ! beta at the new time step
    real(wp)                 :: gamma_new(nlays)        ! gamma at the new time step
    real(wp)                 :: alpha(nlays)            ! alpha
    real(wp)                 :: beta(nlays)             ! beta
    real(wp)                 :: gammaa(nlays)           ! gamma
    real(wp)                 :: impl_weight             ! implicit weight
    real(wp)                 :: c_v                     ! specific heat capacity at constant volume
    real(wp)                 :: c_p                     ! specific heat capacity at constant pressure
    real(wp)                 :: r_d                     ! individual gas constant of dry air
    real(wp)                 :: damping_start_height    ! lower boundary height of the Klemp layer
    real(wp)                 :: damping_coeff           ! damping coefficient of the Klemp layer
    real(wp)                 :: z_above_damping         ! height above the lower boundary of the damping height

    c_v = spec_heat_cap_diagnostics_v(1)
    c_p = spec_heat_cap_diagnostics_p(1)
    r_d = gas_constant_diagnostics(1)
    damping_start_height = klemp_begin_rel*toa

    ! setting the implicit weight
    impl_weight = c_v/c_p

    do ji=1,nlins
      do jk=1,ncols
      
        ! explicit quantities
        do jl=1,nlays
          ! explicit density
          rho_expl(jl) = state_old%rho(ji+1,jk+1,jl) + dtime*tend%rho(ji,jk,jl)
          ! explicit potential temperature density
          rhotheta_expl(jl) = state_old%rhotheta(ji+1,jk+1,jl) + dtime*tend%rhotheta(ji,jk,jl)
          if (rk_Step == 1) then
            ! old time step partial derivatives of rho*theta and Pi (divided by the volume)
            alpha(jl) = -state_old%rhotheta(ji+1,jk+1,jl)/state_old%rho(ji+1,jk+1,jl)**2/grid%volume(ji,jk,jl)
            beta(jl)  = 1._wp/state_old%rho(ji+1,jk+1,jl)/grid%volume(ji,jk,jl)
            gammaa(jl) = r_d/(c_v*state_old%rhotheta(ji+1,jk+1,jl))* &
            (grid%exner_bg(ji+1,jk+1,jl)+state_old%exner_pert(ji+1,jk+1,jl))/grid%volume(ji,jk,jl)
          else
            ! old time step partial derivatives of rho*theta and Pi
            alpha_old(jl) = -state_old%rhotheta(ji+1,jk+1,jl)/state_old%rho(ji+1,jk+1,jl)**2
            beta_old(jl)  = 1._wp/state_old%rho(ji+1,jk+1,jl)
            gamma_old(jl) = r_d/(c_v*state_old%rhotheta(ji+1,jk+1,jl))* &
            (grid%exner_bg(ji+1,jk+1,jl)+state_old%exner_pert(ji+1,jk+1,jl))
            ! new time step partial derivatives of rho*theta and Pi
            alpha_new(jl) = -state_new%rhotheta(ji+1,jk+1,jl)/state_new%rho(ji+1,jk+1,jl)**2
            beta_new(jl)  = 1._wp/state_new%rho(ji+1,jk+1,jl)
            gamma_new(jl) = r_d/(c_v*state_new%rhotheta(ji+1,jk+1,jl)) &
            *(grid%exner_bg(ji+1,jk+1,jl)+state_new%exner_pert(ji+1,jk+1,jl))
            ! interpolation of partial derivatives of rho times theta and Pi
            alpha (jl) = ((1._wp - impl_weight)*alpha_old(jl) + impl_weight*alpha_new(jl))/grid%volume(ji,jk,jl)
            beta  (jl) = ((1._wp - impl_weight)*beta_old (jl) + impl_weight*beta_new (jl))/grid%volume(ji,jk,jl)
            gammaa(jl) = ((1._wp - impl_weight)*gamma_old(jl) + impl_weight*gamma_new(jl))/grid%volume(ji,jk,jl)
          endif
          ! explicit potential temperature perturbation
          theta_pert_expl(jl) = state_old%theta_pert(ji+1,jk+1,jl) + dtime*grid%volume(ji,jk,jl)*(alpha(jl)*tend%rho(ji,jk,jl) &
          + beta(jl)*tend%rhotheta(ji,jk,jl))
          ! explicit Exner pressure perturbation
          exner_pert_expl(jl) = state_old%exner_pert(ji+1,jk+1,jl) + dtime*grid%volume(ji,jk,jl)*gammaa(jl)*tend%rhotheta(ji,jk,jl)
        enddo
        
        ! interface values
        do jl=1,nlays-1
          rho_int_old(jl) = 0.5_wp*(state_old%rho(ji+1,jk+1,jl)+state_old%rho(ji+1,jk+1,jl+1))
          rho_int_expl(jl) = 0.5_wp*(rho_expl(jl)+rho_expl(jl+1))
          theta_int_new(jl) = 0.5_wp*(state_new%rhotheta(ji+1,jk+1,jl)/state_new%rho(ji+1,jk+1,jl) &
          + state_new%rhotheta(ji+1,jk+1,jl+1)/state_new%rho(ji+1,jk+1,jl+1))
        enddo
      
        ! filling up the coefficient vectors
        do jl=1,nlays-1
          ! main diagonal
          d_vector(jl) = -theta_int_new(jl)**2*(gammaa(jl)+gammaa(jl+1)) &
          + 0.5_wp*(grid%exner_bg(ji+1,jk+1,jl)-grid%exner_bg(ji+1,jk+1,jl+1)) &
          *(alpha(jl+1)-alpha(jl)+theta_int_new(jl)*(beta(jl+1)-beta(jl))) &
          - (grid%z_geo_scal(ji+1,jk+1,jl)-grid%z_geo_scal(ji+1,jk+1,jl+1))/(impl_weight*dtime**2*c_p*rho_int_old(jl)) &
          *(2._wp/grid%area_z(ji,jk,jl+1)-dtime*state_old%wind_w(ji+1,jk+1,jl+1)*0.5_wp &
          *(1._wp/grid%volume(ji,jk,jl)+1._wp/grid%volume(ji,jk,jl+1)))
          ! Klemp swamp layer
          z_above_damping = grid%z_geo_w(ji+1,jk+1,jl+1)-damping_start_height
          if (z_above_damping < 0._wp .or. .not. lklemp) then
            damping_coeff = 0._wp
          else
            damping_coeff = klemp_damp_max*sin(0.5_wp*4*atan(1.d0)*z_above_damping/(toa-damping_start_height))**2
          endif
          d_vector(jl) = (1._wp + damping_coeff*dtime)*d_vector(jl)
          ! right hand side
          r_vector(jl) = -(state_old%wind_w(ji+1,jk+1,jl+1)+dtime*tend%wind_w(ji,jk,jl+1))* &
          (grid%z_geo_scal(ji+1,jk+1,jl)-grid%z_geo_scal(ji+1,jk+1,jl+1)) &
          /(impl_weight*dtime**2*c_p) &
          + theta_int_new(jl)*(exner_pert_expl(jl)-exner_pert_expl(jl+1))/dtime &
          + 0.5_wp/dtime*(theta_pert_expl(jl)+theta_pert_expl(jl+1))*(grid%exner_bg(ji+1,jk+1,jl)-grid%exner_bg(ji+1,jk+1,jl+1)) &
          - (grid%z_geo_scal(ji+1,jk+1,jl)-grid%z_geo_scal(ji+1,jk+1,jl+1))/(impl_weight*dtime**2*c_p) &
          *state_old%wind_w(ji+1,jk+1,jl+1)*rho_int_expl(jl)/rho_int_old(jl)
        enddo
        
        do jl=1,nlays-2
          ! lower diagonal
          c_vector(jl) = theta_int_new(jl+1)*gammaa(jl+1)*theta_int_new(jl) &
          + 0.5_wp*(grid%exner_bg(ji+1,jk+1,jl+1)-grid%exner_bg(ji+1,jk+1,jl+2)) &
          *(alpha(jl+1)+beta(jl+1)*theta_int_new(jl)) &
          - (grid%z_geo_scal(ji+1,jk+1,jl+1)-grid%z_geo_scal(ji+1,jk+1,jl+2))/(impl_weight*dtime*c_p)*0.5_wp &
          *state_old%wind_w(ji+1,jk+1,jl+2)/(grid%volume(ji,jk,jl+1)*rho_int_old(jl+1))
          ! upper diagonal
          e_vector(jl) = theta_int_new(jl)*gammaa(jl+1)*theta_int_new(jl+1) &
          - 0.5_wp*(grid%exner_bg(ji+1,jk+1,jl)-grid%exner_bg(ji+1,jk+1,jl+1)) &
          *(alpha(jl+1)+beta(jl+1)*theta_int_new(jl+1)) &
          + (grid%z_geo_scal(ji+1,jk+1,jl)-grid%z_geo_scal(ji+1,jk+1,jl+1))/(impl_weight*dtime*c_p)*0.5_wp &
          *state_old%wind_w(ji+1,jk+1,jl+1)/(grid%volume(ji,jk,jl+1)*rho_int_old(jl))
        enddo
    
        call thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution,nlays-1)
        
        ! results
        ! density, potential temperature density
        do jl=2,nlays-1
          state_new%rho(ji+1,jk+1,jl) = rho_expl(jl) + dtime*(-solution(jl-1)+solution(jl))
          state_new%rhotheta(ji+1,jk+1,jl) = rhotheta_expl(jl) &
          + dtime*(-theta_int_new(jl-1)*solution(jl-1)+theta_int_new(jl)*solution(jl))
        enddo
        ! uppermost layer
        state_new%rho(ji+1,jk+1,1) = rho_expl(1) + dtime*solution(1)
        state_new%rhotheta(ji+1,jk+1,1) = rhotheta_expl(1) + dtime*theta_int_new(1)*solution(1)
        ! lowest layer
        state_new%rho(ji+1,jk+1,nlays) = rho_expl(nlays) - dtime*solution(nlays-1)
        state_new%rhotheta(ji+1,jk+1,nlays) = rhotheta_expl(nlays) - dtime*theta_int_new(nlays-1)*solution(nlays-1)
        ! vertical velocity
        do jl=2,nlays
          rho_int_new = 0.5_wp*(state_new%rho(ji+1,jk+1,jl-1)+state_new%rho(ji+1,jk+1,jl))
          state_new%wind_w(ji+1,jk+1,jl)  = (2._wp*solution(jl-1)/grid%area_z(ji,jk,jl) &
          - rho_int_new*state_old%wind_w(ji+1,jk+1,jl))/rho_int_old(jl-1)
        enddo
        ! Exner pressure
        do jl=1,nlays
          state_new%exner_pert(ji+1,jk+1,jl) = exner_pert_expl(jl) &
          + grid%volume(ji,jk,jl)*gammaa(jl)*(state_new%rhotheta(ji+1,jk+1,jl)-rhotheta_expl(jl))
        enddo
        
      enddo
    enddo
    
    ! potential temperature perturbation at the new time step
    state_new%theta_pert(:,:,:) = state_new%rhotheta(:,:,:)/state_new%rho(:,:,:) - grid%theta_bg(:,:,:)

  end subroutine three_band_solver_ver
  
  subroutine thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,solution_length)

    ! This subroutine solves a system of linear equations with a three-band matrix.
    
    real(wp), intent(in)    :: c_vector(:)
    real(wp), intent(in)    :: d_vector(:)
    real(wp), intent(in)    :: e_vector(:)
    real(wp), intent(in)    :: r_vector(:)
    real(wp), intent(inout) :: solution_vector(:) ! vector containing the solution
    integer,  intent(in)    :: solution_length    ! length of the solution vector
    
    ! local variables
    real(wp) :: e_prime_vector(solution_length-1) ! help vector for solving the matrix equation
    real(wp) :: r_prime_vector(solution_length)   ! help vector for solving the matrix equation
    integer  :: jl                                ! loop index
    
    ! downward sweep (matrix)
    if (d_vector(1) /= 0._wp) then
      e_prime_vector(1) = e_vector(1)/d_vector(1)
    else
      e_prime_vector(1) = 0._wp
    endif
    do jl=2,solution_length-1
      if (d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1) /= 0) then
        e_prime_vector(jl) = e_vector(jl)/(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
      else
        e_prime_vector(jl) = 0._wp
      endif
    enddo
    ! downward sweep (right-hand side)
    if (d_vector(1) /= 0) then
      r_prime_vector(1) = r_vector(1)/d_vector(1)
    else
      r_prime_vector(1) = 0._wp
    endif
    do jl=2,solution_length
      if (d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1) /= 0) then
        r_prime_vector(jl) = (r_vector(jl) - r_prime_vector(jl-1)*c_vector(jl-1)) &
        /(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
      else
        r_prime_vector(jl) = 0._wp
      endif
    enddo
    
    ! upward sweep (final solution)
    solution_vector(solution_length) = r_prime_vector(solution_length)
    do jl=solution_length-1,1,-1
      solution_vector(jl) = r_prime_vector(jl) - e_prime_vector(jl)*solution_vector(jl+1)
    enddo
  
  end subroutine thomas_algorithm

end module vertical_slice_solvers











