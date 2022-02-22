! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module constants

  use definitions, only: wp

  ! This is a collection of some constants.
  
  implicit none
  
  public
  
  real(wp) :: re = 6371000.789927_wp      ! Earth radius
  real(wp) :: k_B = 1.380649e-23_wp       ! Boltzmann's constant
  real(wp) :: T_0 = 273.15_wp             ! 273.15 K
  real(wp) :: EPSILON_SECURITY = 1e-10_wp ! security constant
  
end module constants








