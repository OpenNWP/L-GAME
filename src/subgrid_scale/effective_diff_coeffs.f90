! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module effective_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use run_nml,            only: nlins,ncols,nlays
  use definitions,        only: wp,t_state,t_diag,t_irrev,t_grid
  use diff_nml,           only: diff_h_smag_div
  use derived_quantities, only: calc_diffusion_coeff
  use constituents_nml,   only: no_of_condensed_constituents
  use tke,                only: tke_update
  
  implicit none
  
  private
  
  public :: hori_div_viscosity
  public :: hori_curl_viscosity
  public :: vert_hori_mom_viscosity
  public :: vert_vert_mom_viscosity
  public :: temp_diffusion_coeffs
  public :: mass_diffusion_coeffs
  
  contains
  
  subroutine hori_div_viscosity(state,diag,divergence_h,irrev,grid)
  
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal divergent movements.
    
    ! input arguments and output
    type(t_state), intent(in)    :: state               ! the state variables of the model atmosphere
    type(t_diag),  intent(in)    :: diag                ! diagnostic quantities
    real(wp),      intent(in)    :: divergence_h(:,:,:) ! divergence of the horizontal wind field
    type(t_irrev), intent(inout) :: irrev               ! irreversible quantities
    type(t_grid),  intent(in)    :: grid                ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! computing the Eddy viscosity
    irrev%viscosity_coeff = diff_h_smag_div*grid%mean_velocity_area*abs(divergence_h)
    
    ! calculation of the molecular diffusion coefficient
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          irrev%viscosity_molecular(ji,jk,jl) = calc_diffusion_coeff(diag%temperature_gas(ji,jk,jl), &
          state%rho(ji,jk,jl,no_of_condensed_constituents+1))
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! adding the molecular diffusion coefficient
    irrev%viscosity_coeff = irrev%viscosity_molecular + irrev%viscosity_coeff
    
    ! restricting the values to a maximum to ensure stability
    irrev%viscosity_coeff = min(irrev%viscosity_coeff,irrev%max_diff_h_coeff_turb)
  
  end subroutine hori_div_viscosity
  
  subroutine hori_curl_viscosity()
  
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal curl movements.
  
  end subroutine hori_curl_viscosity
  
  subroutine vert_hori_mom_viscosity(state,diag,irrev,grid)
  
    ! This subroutine computes the effective viscosity (Eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
	! This quantity is located at the half level edges.
	! To obey the symmetry of the stress tensor, the same coefficient must be used for the horizontal diffusion of vertical velocity.
  
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
	
	!  updating the TKE
    call tke_update(state,diag,irrev,grid)
  
  end subroutine vert_hori_mom_viscosity
  
  subroutine vert_vert_mom_viscosity()
  
    ! This subroutine multiplies scalar_field_placeholder (containing dw/dz) by the diffusion coefficient acting on w because of w.
  
  end subroutine vert_vert_mom_viscosity
  
  subroutine temp_diffusion_coeffs()
  
    ! This subroutine computes the viscous temperature diffusion coefficient (including eddies).
  
  end subroutine temp_diffusion_coeffs
  
  subroutine mass_diffusion_coeffs()
  
    ! This subroutine computes the viscous tracer diffusion coefficient (including eddies).
  
  end subroutine mass_diffusion_coeffs
  
  function tke2vertical_diff_coeff(tke)
    
    ! This function returns the vertical kinematic Eddy viscosity as a function of the specific TKE.
	
    ! input
    real(wp), intent(in) :: tke         ! specific turbulent kinetic energy (TKE)
    ! output
    real(wp) :: tke2vertical_diff_coeff ! the result (vertical Eddy viscosity im m^2/s)
    
    ! local variable
    real(wp) :: prop_constant ! semi-empirical constant
	
    prop_constant = 0.4_wp ! unit: m
    ! calculating the result
    tke2vertical_diff_coeff = prop_constant*tke**0.5_wp
	
  end function tke2vertical_diff_coeff
  
end module effective_diff_coeffs







