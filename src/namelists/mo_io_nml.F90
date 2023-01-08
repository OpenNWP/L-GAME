! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_io_nml
  
  ! This nameslist configures the IO behaviour of the model.
  
  use mo_definitions, only: wp
  
  implicit none
  
  integer           :: dt_write_min      ! output interval in minutes
  logical           :: lread_geo         ! wether or not to read the surface properties from a file
  logical           :: lwrite_grid       ! wether or not to write grid properties to a file
  character(len=64) :: grid_filename     ! filename of the grid to read or write
  real(wp)          :: dt_write          ! output interval in seconds
  character(len=64) :: restart_filename  ! filename from which to read the inital state in case restart mode is on
  character(len=64) :: land_sea_filename ! filename of the land-sea mask
  logical           :: lcompute_geo      ! switch for computing the surface properties
  character(len=64) :: oro_raw_filename  ! filename from which to read the raw orography
  logical           :: lwrite_integrals  ! If set to true, fundamental integrals of the atmosphere will be written out at every time step.
  
  namelist /io/dt_write_min,lread_geo,lwrite_grid,grid_filename,restart_filename,land_sea_filename, &
               lcompute_geo,lwrite_integrals
  
  contains
  
  subroutine io_nml_setup
    
    ! local variables
    integer :: fileunit ! file unit of the namelist file
    
    dt_write_min = 60
    lread_geo = .true.
    lwrite_grid = .false.
    grid_filename = "grid.nc"
    restart_filename = "init.nc"
    land_sea_filename = "land_fraction.nc"
    lcompute_geo = .false.
    oro_raw_filename = "etopo.nc"
    lwrite_integrals = .false.
    
    ! Open and read namelist file.
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=io,unit=fileunit)
    
    close(fileunit)
    
    ! calculating the output time step in seconds
    dt_write = 60._wp*dt_write_min
    
    ! sanity check
    if (lcompute_geo .and. lread_geo) then
      write(*,*) "Error: lcompute_geo and lread_geo cannot both be true at the same time."
      call exit(1)
    endif
    
  end subroutine io_nml_setup
  
end module mo_io_nml









