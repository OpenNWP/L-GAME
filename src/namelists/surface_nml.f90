! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module surface_nml

  use definitions, only: wp
  
  implicit none
  
  logical  :: lsoil                   ! soil switch
  integer  :: nsoillays               ! number of soil layers
  real(wp) :: snow_velocity           ! velocity of the falling snow
  real(wp) :: rain_velocity           ! velocity of the falling snow
  real(wp) :: cloud_droplets_velocity ! velocity of the falling snow
  
  
  namelist /surface/lsoil,nsoillays
  
  contains
  
  subroutine surface_nml_setup

    ! local variables
    integer :: fileunit
    
    ! default values
    lsoil                   = .true.
    nsoillays               = 5
    snow_velocity           = 10._wp
    rain_velocity           = 10._wp
    cloud_droplets_velocity = .01_wp
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=surface, unit=fileunit)
        
    close(fileunit)
  
  end subroutine surface_nml_setup

end module surface_nml
