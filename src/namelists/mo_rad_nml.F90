! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_rad_nml
  
  ! In this namelist, the radiation is configured.
  
  use mo_definitions, only: wp
  use mo_run_nml,     only: dy
  
  implicit none
  
  logical            :: lrad                        ! thickness of the horizontal swamp layer
  real(wp)           :: dtime_rad                   ! radiation time step
  character(len=128) :: rrtmgp_coefficients_file_sw ! name of the shortwave data file
  character(len=128) :: rrtmgp_coefficients_file_lw ! name of the longwave data file
  character(len=128) :: cloud_coefficients_file_sw  ! name of the shortwave cloud optics file
  character(len=128) :: cloud_coefficients_file_lw  ! name of the longwave cloud optics file
  
  namelist /rad/lrad,dtime_rad,rrtmgp_coefficients_file_sw,rrtmgp_coefficients_file_lw, &
                cloud_coefficients_file_sw,cloud_coefficients_file_lw
  
  contains
  
  subroutine rad_nml_setup
    
    ! local variables
    integer :: fileunit ! file unit of the namelist file
    
    ! default values
    lrad = .true.
    dtime_rad = 60._wp*dy/1000._wp
    rrtmgp_coefficients_file_sw = "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-sw-g112-210809.nc"
    rrtmgp_coefficients_file_lw = "/home/max/code/rte-rrtmgp/rrtmgp/data/rrtmgp-data-lw-g128-210809.nc"
    cloud_coefficients_file_sw = "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc"
    cloud_coefficients_file_lw = "/home/max/code/rte-rrtmgp/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc"
    
    ! Open and read namelist file.
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=rad,unit=fileunit)
    
    close(fileunit)
    
  end subroutine rad_nml_setup
  
end module mo_rad_nml
