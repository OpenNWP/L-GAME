! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module effective_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use definitions, only: wp
  
  implicit none
  
  private
  
  public :: hori_div_viscosity
  public :: hori_curl_viscosity
  public :: vert_hori_mom_viscosity
  public :: vert_vert_mom_viscosity
  public :: temp_diffusion_coeffs
  public :: mass_diffusion_coeffs
  
  contains
  
  subroutine hori_div_viscosity()
  
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal divergent movements.
  
  end subroutine hori_div_viscosity
  
  subroutine hori_curl_viscosity()
  
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal curl movements.
  
  end subroutine hori_curl_viscosity
  
  subroutine vert_hori_mom_viscosity()
  
    ! This subroutine computes the effective viscosity (Eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
	! This quantity is located at the half level edges.
	! To obey the symmetry of the stress tensor, the same coefficient must be used for the horizontal diffusion of vertical velocity.
  
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
    real(wp), intent(in) :: tke
	
    ! output
    real(wp)             :: tke2vertical_diff_coeff
    
    ! local variable
    real(wp)             :: prop_constant
	
    prop_constant = 0.4_wp ! unit: m
    tke2vertical_diff_coeff = prop_constant*tke**0.5_wp
	
  end function tke2vertical_diff_coeff
  
end module effective_diff_coeffs







