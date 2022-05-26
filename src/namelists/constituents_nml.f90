! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module constituents_nml

  ! This namelist defines the constituents of the model atmosphere.

  use definitions, only: wp
  use run_nml,     only: lmoist

  implicit none
  
  integer  :: no_of_gaseous_constituents   ! number of constituents of the gas phase
  integer  :: no_of_condensed_constituents ! number of condensed constituents
  integer  :: no_of_constituents           ! the total number of constituents
  real(wp) :: snow_velocity                ! sedimentation velocity of snow
  real(wp) :: rain_velocity                ! sedimentation velocity of rain
  real(wp) :: cloud_droplets_velocity      ! sedimentation velocity of cloud droplets

  contains

  subroutine constituents_nml_setup()
    
    no_of_condensed_constituents = 4
    no_of_gaseous_constituents = 2
    ! the dry case
    if (.not. lmoist) then
      no_of_condensed_constituents = 0
      no_of_gaseous_constituents = 1
    endif
    snow_velocity = 5._wp
    rain_velocity = 10._wp
    cloud_droplets_velocity = .01_wp
    no_of_constituents = no_of_condensed_constituents + no_of_gaseous_constituents
    
  end subroutine constituents_nml_setup
  
end module constituents_nml






