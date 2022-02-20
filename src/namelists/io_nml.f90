! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module io_nml

  use definitions, only: wp
  
  implicit none
  
  integer           :: dt_write_min       ! output interval in minutes
  logical           :: read_grid          ! wether or not to read the grid from a file
  logical           :: write_grid         ! wether or not to write the grid to a file
  character(len=64) :: grid_filename      ! filename of the grid to read or write
  real(wp)          :: dt_write           ! output interval in seconds
  
  namelist /io/dt_write_min,read_grid,write_grid,grid_filename

  contains

  subroutine io_nml_setup
  
    ! local variables
    integer :: fileunit
    
    dt_write_min    = 60
    read_grid       = .false.
    write_grid      = .false.
    grid_filename   = "grid.nc"
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=io, unit=fileunit)
        
    close(fileunit)
    
    ! calculating the output timestep in seconds
    dt_write        = 60._wp*dt_write_min
  
  end subroutine io_nml_setup
  
end module io_nml









