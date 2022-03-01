! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module io_nml

  ! This nameslist configures the IO behaviour of the model.

  use definitions, only: wp
  
  implicit none
  
  integer           :: dt_write_min      ! output interval in minutes
  logical           :: lread_oro         ! wether or not to read the orography from a file
  logical           :: lwrite_grid       ! wether or not to write grid properties to a file
  character(len=64) :: grid_filename     ! filename of the grid to read or write
  real(wp)          :: dt_write          ! output interval in seconds
  character(len=64) :: restart_filename  ! filename from which to read the inital state in case restart mode is on
  logical           :: lread_land_sea    ! switch for reading the land-sea mask
  character(len=64) :: land_sea_filename ! filename of the land-sea mask
  
  namelist /io/dt_write_min,lread_oro,lwrite_grid,grid_filename,restart_filename,lread_land_sea,land_sea_filename

  contains

  subroutine io_nml_setup
  
    ! local variables
    integer :: fileunit
    
    dt_write_min = 60
    lread_oro = .true.
    lwrite_grid = .false.
    grid_filename = "grid.nc"
    restart_filename = "init.nc"
    lread_land_sea = .true.
    land_sea_filename = "is_land.nc"
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=io, unit=fileunit)
    
    close(fileunit)
    
    ! calculating the output timestep in seconds
    dt_write = 60._wp*dt_write_min
  
  end subroutine io_nml_setup
  
end module io_nml









