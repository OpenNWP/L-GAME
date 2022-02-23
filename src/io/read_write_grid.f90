! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module read_write_grid

  ! This module reads the grid from a file or writes the grid to a file. This is useful for efficiency.
  
  use netcdf
  use definitions,       only: t_grid
  use set_initial_state, only: nc_check
  use io_nml,            only: grid_filename
  use run_nml,           only: nlins,ncols
  
  implicit none
  
  private
  
  public :: write_grid
  public :: read_grid
  
  contains
  
  subroutine write_grid(grid)
  
    type(t_grid), intent(in) :: grid
    
    ! local variables
    integer                   :: ncid                      ! ID of the NetCDF file
    character(len=64)         :: filename                  ! output filename
    integer                   :: x_dimid                   ! ID of the x dimension
    integer                   :: y_dimid                   ! ID of the y dimension
    integer                   :: dimids(2)                 ! dimensions of horizontal fields
    integer                   :: varid_lat                 ! variable ID of the latitudes
    integer                   :: varid_lon                 ! variable ID of the longitudes
    
    filename = "../../grids/" // trim(grid_filename)
    
    ! creating the NetCDF file
    call nc_check(nf90_create(trim(filename),NF90_CLOBBER,ncid))
    
    call nc_check(nf90_put_att(ncid,NF90_GLOBAL,"Description","This is a grid of L-GAME."))
    
    ! defining the dimensions
    call nc_check(nf90_def_dim(ncid,"lon_model",ncols,x_dimid))
    call nc_check(nf90_def_dim(ncid,"lat_model",nlins,y_dimid))

    ! setting the dimension ID arrays
    ! 2D
    dimids(1) = y_dimid
    dimids(2) = x_dimid
    
    call nc_check(nf90_def_var(ncid,"lat_geo",NF90_REAL,dimids,varid_lat))
    call nc_check(nf90_put_att(ncid,varid_lat,"Description","latitude"))
    call nc_check(nf90_put_att(ncid,varid_lat,"Unit","radians"))
    call nc_check(nf90_def_var(ncid,"lon_geo",NF90_REAL,dimids,varid_lon))
    call nc_check(nf90_put_att(ncid,varid_lon,"Description","longitude"))
    call nc_check(nf90_put_att(ncid,varid_lon,"Unit","radians"))
    
    ! ending the definition section
    call nc_check(nf90_enddef(ncid))
    
    ! writing the data to the NetCDF file
    call nc_check(nf90_put_var(ncid,varid_lat,grid%lat_geo_scalar))
    call nc_check(nf90_put_var(ncid,varid_lon,grid%lon_geo_scalar))
    
    ! closing the NetCDF file
    call nc_check(nf90_close(ncid))
  
  end subroutine write_grid
  
  subroutine read_grid()
  
    ! This subroutine reads important grid properties from a file.
  
  end subroutine read_grid

end module read_write_grid








