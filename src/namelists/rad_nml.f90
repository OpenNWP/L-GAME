! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module rad_nml

  ! In this namelist, the radiation is configured.
  
  use definitions, only: wp
  use run_nml,     only: dy
  
  implicit none
  
  logical            :: lrad                        ! thickness of the horizontal swamp layer
  real(wp)           :: dtime_rad                   ! radiation timestep
  character(len=128) :: rrtmgp_coefficients_file_sw ! name of the short wave data file
  character(len=128) :: rrtmgp_coefficients_file_lw ! name of the long wave data file
  character(len=128) :: cloud_coefficients_file_sw  ! name of the short wave cloud optics file
  character(len=128) :: cloud_coefficients_file_lw  ! name of the long wave cloud optics file
  
  namelist /rad/lrad,dtime_rad
  
  contains
  
  subroutine rad_nml_setup

    ! local variables
    integer :: fileunit
    
    ! default values
    lrad = .true.
    dtime_rad = 60._wp*dy/1000._wp
    rrtmgp_coefficients_file_sw = "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g112-210809.nc"
    rrtmgp_coefficients_file_lw = "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g128-210809.nc"
    cloud_coefficients_file_sw = "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc"
    cloud_coefficients_file_lw = "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc"
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=rad, unit=fileunit)
        
    close(fileunit)
  
  end subroutine rad_nml_setup

end module rad_nml
