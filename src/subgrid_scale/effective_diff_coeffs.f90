! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module effective_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use run_nml,              only: nlins,ncols,nlays,dtime
  use definitions,          only: wp,t_state,t_diag,t_irrev,t_grid
  use diff_nml,             only: diff_h_smag_div,diff_h_smag_rot,lmom_diff_h,ltemp_diff_h
  use derived_quantities,   only: calc_diffusion_coeff
  use constituents_nml,     only: no_of_condensed_constituents
  use tke,                  only: tke_update
  use divergence_operators, only: div_h
  use derived_quantities,   only: density_gas,spec_heat_cap_diagnostics_v
  use bc_nml,               only: lperiodic
  
  implicit none
  
  private
  
  public :: hor_div_viscosity
  public :: hor_curl_viscosity
  public :: vert_hor_mom_viscosity
  public :: vert_vert_mom_viscosity
  public :: temp_diffusion_coeffs
  public :: mass_diffusion_coeffs
  
  contains
  
  subroutine hor_div_viscosity(state,diag,divergence_h,irrev,grid)
  
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal divergent movements.
    
    ! input arguments and output
    type(t_state), intent(in)    :: state               ! the state variables of the model atmosphere
    type(t_diag),  intent(in)    :: diag                ! diagnostic quantities
    real(wp),      intent(in)    :: divergence_h(:,:,:) ! divergence of the horizontal wind field
    type(t_irrev), intent(inout) :: irrev               ! irreversible quantities
    type(t_grid),  intent(in)    :: grid                ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! computing the eddy viscosity
    irrev%viscosity_coeff_div = diff_h_smag_div*grid%mean_velocity_area*abs(divergence_h)
    
    ! calculation of the molecular diffusion coefficient
    !$omp parallel do private(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          irrev%viscosity_molecular(ji,jk,jl) = calc_diffusion_coeff(diag%temperature(ji,jk,jl), &
          state%rho(ji,jk,jl,no_of_condensed_constituents+1))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! adding the molecular diffusion coefficient
    !$omp parallel workshare
    irrev%viscosity_coeff_div = irrev%viscosity_molecular + irrev%viscosity_coeff_div
    !$omp end parallel workshare
    
    ! multiplying by the density
    !$omp parallel do private(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          irrev%viscosity_coeff_div(ji,jk,jl) = density_gas(state,ji,jk,jl)*irrev%viscosity_coeff_div(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine hor_div_viscosity
  
  subroutine hor_curl_viscosity(state,diag,irrev,grid)
  
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal curl movements.
  
    ! input arguments and output
    type(t_state), intent(in)    :: state ! the state variables of the model atmosphere
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! initialization with zeros
    !$omp parallel workshare
     irrev%viscosity_coeff_curl_dual = 0._wp
    !$omp end parallel workshare
    
    ! molecular component
    !$omp parallel do private(ji,jk,jl)
    do jl=1,nlays
      do ji=2,nlins
        do jk=2,ncols
          irrev%viscosity_coeff_curl_dual(ji,jk,jl) = 0.25_wp*sum(irrev%viscosity_molecular(ji-1:ji,jk-1:jk,jl))
        enddo
      enddo
      
      ! periodic boundary conditions
      if (lperiodic) then
      
        do jk=2,ncols
          irrev%viscosity_coeff_curl_dual(1,jk,jl) = 0.25_wp*( &
          irrev%viscosity_molecular(1,jk-1,jl) + irrev%viscosity_molecular(nlins,jk-1,jl) &
          + irrev%viscosity_molecular(1,jk,jl) + irrev%viscosity_molecular(nlins,jk,jl))
          irrev%viscosity_coeff_curl_dual(nlins+1,jk,jl) = irrev%viscosity_coeff_curl_dual(1,jk,jl)
        enddo
        do ji=2,nlins
          irrev%viscosity_coeff_curl_dual(ji,1,jl) = 0.25_wp*( &
          irrev%viscosity_molecular(ji-1,ncols,jl) + irrev%viscosity_molecular(ji-1,1,jl) &
          + irrev%viscosity_molecular(ji,ncols,jl) + irrev%viscosity_molecular(ji,1,jl))
          irrev%viscosity_coeff_curl_dual(ji,ncols+1,jl) = irrev%viscosity_coeff_curl_dual(ji,1,jl)
        enddo
      
        ! corners
        irrev%viscosity_coeff_curl_dual(1,1,jl) = 0.25*(irrev%viscosity_molecular(1,1,jl) + irrev%viscosity_molecular(nlins,1,jl) &
        + irrev%viscosity_molecular(1,ncols,jl) + irrev%viscosity_molecular(nlins,ncols,jl))
        irrev%viscosity_coeff_curl_dual(1,ncols+1,jl) = irrev%viscosity_coeff_curl_dual(1,1,jl)
        irrev%viscosity_coeff_curl_dual(nlins+1,1,jl) = irrev%viscosity_coeff_curl_dual(1,1,jl)
        irrev%viscosity_coeff_curl_dual(nlins+1,ncols+1,jl) = irrev%viscosity_coeff_curl_dual(1,1,jl)
        
      endif
    
    enddo
    !$omp end parallel do
    
    ! turbulent component
    !$omp parallel do private(ji,jk,jl)
    do ji=1,nlins+1
      do jk=1,ncols+1
        do jl=1,nlays
          irrev%viscosity_coeff_curl_dual(ji,jk,jl) = irrev%viscosity_coeff_curl_dual(ji,jk,jl) &
          + diff_h_smag_rot*grid%mean_velocity_area*abs(diag%zeta_z(ji,jk,jl))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! multiplication by the density
    !$omp parallel do private(ji,jk,jl)
    do jl=1,nlays
      do ji=2,nlins
        do jk=2,ncols
          irrev%viscosity_coeff_curl_dual(ji,jk,jl) = 0.25_wp*(density_gas(state,ji-1,jk-1,jl)+density_gas(state,ji-1,jk,jl) &
          + density_gas(state,ji,jk-1,jl) + density_gas(state,ji,jk,jl))*irrev%viscosity_coeff_curl_dual(ji,jk,jl)
        enddo
      enddo
      
      ! periodic boundary conditions
      if (lperiodic) then
      
        do jk=2,ncols
          irrev%viscosity_coeff_curl_dual(1,jk,jl) = irrev%viscosity_coeff_curl_dual(1,jk,jl) &
          *(0.25_wp*(density_gas(state,1,jk-1,jl) + density_gas(state,nlins,jk-1,jl) &
          + density_gas(state,1,jk,jl) + density_gas(state,nlins,jk,jl)))
          irrev%viscosity_coeff_curl_dual(nlins+1,jk,jl) = irrev%viscosity_coeff_curl_dual(1,jk,jl)
        enddo
        do ji=2,nlins
          irrev%viscosity_coeff_curl_dual(ji,1,jl) = irrev%viscosity_coeff_curl_dual(ji,1,jl) &
          *(0.25_wp*(density_gas(state,ji-1,ncols,jl)+density_gas(state,ji-1,1,jl) &
          + density_gas(state,ji,ncols,jl)+density_gas(state,ji,1,jl)))
          irrev%viscosity_coeff_curl_dual(ji,ncols+1,jl) = irrev%viscosity_coeff_curl_dual(ji,1,jl)
        enddo
      
        ! corners
        irrev%viscosity_coeff_curl_dual(1,1,jl) = irrev%viscosity_coeff_curl_dual(1,1,jl) &
        *(0.25_wp*(density_gas(state,nlins,ncols,jl)+density_gas(state,nlins,1,jl) &
        + density_gas(state,1,ncols,jl)+density_gas(state,1,1,jl)))
        irrev%viscosity_coeff_curl_dual(1,ncols+1,jl) = irrev%viscosity_coeff_curl_dual(1,1,jl)
        irrev%viscosity_coeff_curl_dual(nlins+1,1,jl) = irrev%viscosity_coeff_curl_dual(1,1,jl)
        irrev%viscosity_coeff_curl_dual(nlins+1,ncols+1,jl) = irrev%viscosity_coeff_curl_dual(1,1,jl)
        
      endif
    
    enddo
    !$omp end parallel do
    
    ! averaging the curl diffusion coefficient to the cell centers
    !$omp parallel do private(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          irrev%viscosity_coeff_curl(ji,jk,jl) = 0.25_wp*sum(irrev%viscosity_coeff_curl_dual(ji:ji+1,jk:jk+1,jl))
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine hor_curl_viscosity
  
  subroutine vert_hor_mom_viscosity(state,diag,irrev,grid)
  
    ! This subroutine computes the effective viscosity (eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
	! This quantity is located at the half level edges.
	! To obey the symmetry of the stress tensor, the same coefficient must be used for the horizontal diffusion of vertical velocity.
  
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl               ! loop indices
	
	!  updating the TKE
    call tke_update(state,diag,irrev,grid)
    
    !$omp parallel do private(ji,jk,jl)
    do ji=1,nlins
      do jl=2,nlays
        do jk=2,ncols
          irrev%vert_hor_viscosity_u(ji,jk,jl) = &
          ! molecular component
          0.25_wp*( &
          irrev%viscosity_molecular(ji,jk-1,jl-1) + irrev%viscosity_molecular(ji,jk,jl-1) &
          + irrev%viscosity_molecular(ji,jk-1,jl) + irrev%viscosity_molecular(ji,jk,jl)) &
          ! turbulent component
          + 0.25_wp*( &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,jk-1,jl-1)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,jk,jl-1)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,jk-1,jl)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,jk,jl)))
          
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          irrev%vert_hor_viscosity_u(ji,1,jl) = &
          ! molecular component
          0.25_wp*( &
          irrev%viscosity_molecular(ji,ncols,jl-1) + irrev%viscosity_molecular(ji,1,jl-1) &
          + irrev%viscosity_molecular(ji,ncols,jl) + irrev%viscosity_molecular(ji,1,jl)) &
          ! turbulent component
          + 0.25_wp*( &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,ncols,jl-1)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,1,jl-1)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,ncols,jl)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,1,jl)))
          
          irrev%vert_hor_viscosity_u(ji,ncols+1,jl) = irrev%vert_hor_viscosity_u(ji,1,jl)
          
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl)
    do jk=1,ncols
      do jl=2,nlays
        do ji=2,nlins
          irrev%vert_hor_viscosity_v(ji,jk,jl) = &
          ! molecular component
          0.25_wp*( &
          irrev%viscosity_molecular(ji-1,jk,jl-1) + irrev%viscosity_molecular(ji,jk,jl-1) &
          + irrev%viscosity_molecular(ji-1,jk,jl) + irrev%viscosity_molecular(ji,jk,jl)) &
          ! turbulent component
          + 0.25_wp*( &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji-1,jk,jl-1)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,jk,jl-1)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji-1,jk,jl)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(ji,jk,jl)))
          
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          irrev%vert_hor_viscosity_v(1,jk,jl) = &
          ! molecular component
          0.25_wp*( &
          irrev%viscosity_molecular(nlins,jk,jl-1) + irrev%viscosity_molecular(1,jk,jl-1) &
          + irrev%viscosity_molecular(nlins,jk,jl) + irrev%viscosity_molecular(1,jk,jl)) &
          ! turbulent component
          + 0.25_wp*( &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(nlins,jk,jl-1)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(1,jk,jl-1)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(nlins,jk,jl)) + &
          tke2vertical_diff_coeff(irrev%viscosity_molecular(1,jk,jl)))
          
          irrev%vert_hor_viscosity_v(nlins+1,jk,jl) = irrev%vert_hor_viscosity_v(1,jk,jl)
          
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    ! multiplication by the density
    !$omp parallel do private(ji,jk,jl)
    do ji=1,nlins
      do jl=2,nlays
        do jk=2,ncols
          irrev%vert_hor_viscosity_u(ji,jk,jl) = irrev%vert_hor_viscosity_u(ji,jk,jl) &
          ! molecular component
          *0.25_wp*( &
          density_gas(state,ji,jk-1,jl-1) + density_gas(state,ji,jk,jl-1) &
          + density_gas(state,ji,jk-1,jl) + density_gas(state,ji,jk,jl))
          
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          irrev%vert_hor_viscosity_u(ji,1,jl) = irrev%vert_hor_viscosity_u(ji,1,jl) &
          ! molecular component
          *0.25_wp*( &
          density_gas(state,ji,ncols,jl-1) + density_gas(state,ji,1,jl-1) &
          + density_gas(state,ji,ncols,jl) + density_gas(state,ji,1,jl))
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl)
    do jk=1,ncols
      do jl=2,nlays
        do ji=2,nlins
          irrev%vert_hor_viscosity_v(ji,jk,jl) = irrev%vert_hor_viscosity_v(ji,jk,jl) &
          ! molecular component
          *0.25_wp*( &
          density_gas(state,ji-1,jk,jl-1) + density_gas(state,ji,jk,jl-1) &
          + irrev%viscosity_molecular(ji-1,jk,jl) + density_gas(state,ji,jk,jl))
          
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          irrev%vert_hor_viscosity_v(1,jk,jl) = irrev%vert_hor_viscosity_v(1,jk,jl) &
          ! molecular component
          *0.25_wp*( &
          density_gas(state,nlins,jk,jl-1) + density_gas(state,1,jk,jl-1) &
          + density_gas(state,nlins,jk,jl) + density_gas(state,1,jk,jl))
          
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel workshare
    ! for now, we set the vertical diffusion coefficient at the TOA equal to the vertical diffusion coefficient in the layer below
    irrev%vert_hor_viscosity_u(:,:,1) = irrev%vert_hor_viscosity_u(:,:,2)
    ! for now, we set the vertical diffusion coefficient at the surface equal to the vertical diffusion coefficient in the layer above
    irrev%vert_hor_viscosity_u(:,:,nlays+1) = irrev%vert_hor_viscosity_u(:,:,nlays)
    !$omp end parallel workshare
    !$omp parallel workshare
    ! for now, we set the vertical diffusion coefficient at the TOA equal to the vertical diffusion coefficient in the layer below
    irrev%vert_hor_viscosity_v(:,:,1) = irrev%vert_hor_viscosity_v(:,:,2)
    ! for now, we set the vertical diffusion coefficient at the surface equal to the vertical diffusion coefficient in the layer above
    irrev%vert_hor_viscosity_v(:,:,nlays+1) = irrev%vert_hor_viscosity_v(:,:,nlays)
    !$omp end parallel workshare
  
  end subroutine vert_hor_mom_viscosity
  
  subroutine vert_vert_mom_viscosity(state,diag,irrev)
  
    ! This subroutine multiplies scalar_placeholder (containing dw/dz) by the diffusion coefficient acting on w because of w.
   
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    
    ! local variables
    real(wp) :: mom_diff_coeff        ! the diffusion coefficient
    integer  :: ji,jk,jl              ! loop indices
    
    !$omp parallel do private(ji,jk,jl,mom_diff_coeff)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
    
          mom_diff_coeff &
          ! molecular viscosity
          = irrev%viscosity_molecular(ji,jk,jl) &
          ! turbulent component
          + tke2vertical_diff_coeff(irrev%tke(ji,jk,jl))

          diag%scalar_placeholder(ji,jk,jl) = density_gas(state,ji,jk,jl)*mom_diff_coeff*diag%scalar_placeholder(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine vert_vert_mom_viscosity
  
  subroutine temp_diffusion_coeffs(state,diag,irrev,grid)
  
    ! This function computes the viscous temperature diffusion coefficient (including eddies).
  
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl  ! loop variables
    real(wp) :: c_g_v    ! specific heat capacity
    
    ! The eddy viscosity coefficient and the TKE only has to be calculated if it has not yet been done.
    if (.not. lmom_diff_h) then
    
      call div_h(state%wind_u,state%wind_v,diag%scalar_placeholder,grid)
      call hor_div_viscosity(state,diag,diag%scalar_placeholder,irrev,grid)
      call hor_curl_viscosity(state,diag,irrev,grid)
      call tke_update(state,diag,irrev,grid)
      
      ! molecular viscosity
      !$omp parallel do private(ji,jk,jl)
      do ji=1,nlins
        do jk=1,ncols
          do jl=1,nlays
            irrev%viscosity_molecular(ji,jk,jl) = calc_diffusion_coeff(diag%temperature(ji,jk,jl), &
            state%rho(ji,jk,jl,no_of_condensed_constituents+1))
          enddo
        enddo
      enddo
      !$omp end parallel do
    
    endif
    
    !$omp parallel do private(ji,jk,jl,c_g_v)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          c_g_v = spec_heat_cap_diagnostics_v(state,ji,jk,jl)
          ! horizontal diffusion coefficient
          irrev%scalar_diff_coeff_h(ji,jk,jl) = c_g_v &
          *0.5_wp*(irrev%viscosity_coeff_div(ji,jk,jl) + irrev%viscosity_coeff_curl(ji,jk,jl))
          ! vertical diffusion coefficient
          irrev%scalar_diff_coeff_v(ji,jk,jl) &
          ! molecular component
          = density_gas(state,ji,jk,jl)*c_g_v*(irrev%viscosity_molecular(ji,jk,jl) &
          ! turbulent component
          + tke2vertical_diff_coeff(irrev%tke(ji,jk,jl)))
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine temp_diffusion_coeffs
  
  subroutine mass_diffusion_coeffs(state,diag,irrev,grid)
  
    ! This subroutine computes the viscous tracer diffusion coefficient (including eddies).
  
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl ! loop variables
    
    ! The eddy viscosity coefficient and the TKE only has to be calculated if it has not yet been done.
    if (.not. lmom_diff_h .and. .not. ltemp_diff_h) then
    
      call div_h(state%wind_u,state%wind_v,diag%scalar_placeholder,grid)
      call hor_div_viscosity(state,diag,diag%scalar_placeholder,irrev,grid)
      call hor_curl_viscosity(state,diag,irrev,grid)
      call tke_update(state,diag,irrev,grid)
      
      ! molecular viscosity
      !$omp parallel do private(ji,jk,jl)
      do ji=1,nlins
        do jk=1,ncols
          do jl=1,nlays
            irrev%viscosity_molecular(ji,jk,jl) = calc_diffusion_coeff(diag%temperature(ji,jk,jl), &
            state%rho(ji,jk,jl,no_of_condensed_constituents+1))
          enddo
        enddo
      enddo
      !$omp end parallel do
    
    endif
    
    !$omp parallel do private(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          ! horizontal diffusion coefficient
          irrev%scalar_diff_coeff_h(ji,jk,jl) &
          = 0.5_wp*(irrev%viscosity_coeff_div(ji,jk,jl) + irrev%viscosity_coeff_curl(ji,jk,jl)) &
          /density_gas(state,ji,jk,jl)
          ! vertical diffusion coefficient
          irrev%scalar_diff_coeff_v(ji,jk,jl) &
          ! molecular component
          = irrev%viscosity_molecular(ji,jk,jl) &
          ! turbulent component
          + tke2vertical_diff_coeff(irrev%tke(ji,jk,jl))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine mass_diffusion_coeffs
  
  function tke2vertical_diff_coeff(tke)
    
    ! This function returns the vertical kinematic eddy viscosity as a function of the specific TKE.
	
    ! input
    real(wp), intent(in) :: tke                     ! specific turbulent kinetic energy (TKE)
    ! output
    real(wp)             :: tke2vertical_diff_coeff ! the result (vertical eddy viscosity im m^2/s)
    
    ! local variable
    real(wp) :: prop_constant ! semi-empirical constant
	
    prop_constant = 0.4_wp ! unit: m
    ! calculating the result
    tke2vertical_diff_coeff = prop_constant*tke**0.5_wp
	
  end function tke2vertical_diff_coeff
  
end module effective_diff_coeffs







