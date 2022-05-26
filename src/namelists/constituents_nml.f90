! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module constituents_nml

  ! This namelist defines the constituents of the model atmosphere.

  use definitions, only: wp

  implicit none
  
  integer  :: no_of_gaseous_constituents   ! number of constituents of the gas phase
  integer  :: no_of_condensed_constituents ! number of condensed constituents
  integer  :: no_of_constituents           ! the total number of constituents
  real(wp) :: snow_velocity                ! sedimentation velocity of snow
  real(wp) :: rain_velocity                ! sedimentation velocity of rain
  real(wp) :: cloud_droplets_velocity      ! sedimentation velocity of cloud droplets  
  
  public :: no_of_constituents
    
  namelist /constituents/no_of_condensed_constituents,no_of_gaseous_constituents, &
  snow_velocity,rain_velocity,cloud_droplets_velocity

  contains

  subroutine constituents_nml_setup()
  
    ! local variables
    integer        :: fileunit
    
    ! default values
    no_of_condensed_constituents = 4
    no_of_gaseous_constituents = 2
    snow_velocity = 5._wp
    rain_velocity = 10._wp
    cloud_droplets_velocity = .01_wp
    
    ! open and read namelist file
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=constituents, unit=fileunit)
    
    close(fileunit)
    
    no_of_constituents = no_of_condensed_constituents + no_of_gaseous_constituents
    
  end subroutine constituents_nml_setup
  
end module constituents_nml






