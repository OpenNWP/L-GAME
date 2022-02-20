! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module effective_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use definitions, only: wp
  
  implicit none
  
  private
  
  contains
  
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
