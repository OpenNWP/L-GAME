! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module constants

  ! This is a collection of some quantities that are hardly ever changed.
  
  use definitions, only: wp
  
  implicit none
  
  public
  
  real(wp) :: re = 6371000.789927_wp      ! Earth radius
  real(wp) :: k_B = 1.380649e-23_wp       ! Boltzmann's constant
  real(wp) :: T_0 = 273.15_wp             ! 273.15 K
  real(wp) :: density_water = 1024._wp    ! typical density of water
  real(wp) :: p_0 = 100000._wp            ! reference pressure
  real(wp) :: omega = 7.292115e-5_wp      ! angular frequency of Earth rotation
  real(wp) :: gravity = 9.80616_wp        ! average surface gravity value
  
  ! non-physical constants
  ! ----------------------
  real(wp) :: M_PI = 4._wp*atan(1._wp)    ! pi
  real(wp) :: EPSILON_SECURITY = 1e-10_wp ! security constant
  
  ! some properties of the standard atmosphere
  ! ------------------------------------------
  real(wp) :: lapse_rate = 0.0065_wp      ! lapse_rate within the troposphere
  real(wp) :: surface_temp = 288.15_wp    ! the temperature at the surface
  real(wp) :: tropo_height = 12000._wp    ! the tropopause height
  real(wp) :: inv_height = 20000._wp      ! height where the temperature inversion begins
  real(wp) :: t_grad_inv = 0.001_wp       ! temperature gradient above the inversion
  real(wp) :: p_0_standard = 101325._wp   ! reference pressure of the standard atmosphere

end module constants








