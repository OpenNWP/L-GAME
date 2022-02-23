! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module rad_nml

  ! In this namelist, the radiation is configured.
  
  use definitions, only: wp
  
  implicit none
  
  logical  :: lrad      ! thickness of the horizontal swamp layer
  real(wp) :: dtime_rad ! radiation timestep
  
  namelist /rad/lrad,dtime_rad
  
  contains
  
  subroutine rad_nml_setup

    ! local variables
    integer :: fileunit
    
    ! default values
    lrad = .true.
    dtime_rad = 3600._wp
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=rad, unit=fileunit)
        
    close(fileunit)
  
  end subroutine rad_nml_setup

end module rad_nml
