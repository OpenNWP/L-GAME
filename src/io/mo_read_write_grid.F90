! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_read_write_grid

  ! This module reads the grid from a file or writes the grid to a file. This is useful for efficiency.
  
  use netcdf
  use mo_definitions,       only: t_grid,wp
  use mo_set_initial_state, only: nc_check
  use mo_io_nml,            only: grid_filename
  use mo_run_nml,           only: ny,nx,n_layers,n_levels
  
  implicit none
  
  contains
  
  subroutine write_grid(grid)
    
    ! This subroutine writes important grid properties to a file so they can be reused elsewhere.
    
    ! input arguments
    type(t_grid), intent(in) :: grid ! grid properties
    
    ! local variables
    integer           :: ncid                   ! ID of the netCDF file
    character(len=64) :: filename               ! output filename
    integer           :: x_dimid                ! ID of the x dimension
    integer           :: y_dimid                ! ID of the y dimension
    integer           :: x_dimidp1              ! ID of the x dimension + one point
    integer           :: y_dimidp1              ! ID of the y dimension + one point
    integer           :: dimids(2)              ! dimensions of horizontal fields
    integer           :: dimids_u(2)            ! dimensions of horizontal u-vector fields
    integer           :: dimids_v(2)            ! dimensions of horizontal v-vector fields
    integer           :: varid_lat              ! variable ID of the latitudes
    integer           :: varid_lon              ! variable ID of the longitudes
    integer           :: varid_lat_u            ! variable ID of the latitudes of the u-vectors
    integer           :: varid_lon_u            ! variable ID of the longitudes of the u-vectors
    integer           :: varid_lat_v            ! variable ID of the latitudes of the v-vectors
    integer           :: varid_lon_v            ! variable ID of the longitudes of the v-vectors
    integer           :: varid_land_fraction    ! variable ID of the land fractionn
    integer           :: varid_lake_fraction    ! variable ID of the lake fraction
    integer           :: varid_roughness_length ! variable ID of the roughness length
    integer           :: varid_sfc_albedo       ! variable ID of the surface albedo
    integer           :: varid_sfc_rho_c        ! variable ID of the surface volumetric heat capacity
    integer           :: varid_t_conduc_soil    ! variable ID of the surface temperature conductivity of the soil
    integer           :: varid_oro              ! variable ID of the orography
    integer           :: varid_oro_smoothed     ! variable ID of the smoothed orography
    integer           :: varid_dir_geo_u        ! variable ID of the direction of u-vectors
    integer           :: varid_dir_geo_v        ! variable ID of the direction of v-vectors
    integer           :: varid_dir_geo_u_scalar ! variable ID of the direction of u-vectors at the scalar data points
    
    filename = "../../grids/" // trim(grid_filename)
    
    ! creating the netCDF file
    call nc_check(nf90_create(trim(filename),NF90_CLOBBER,ncid))
    
    call nc_check(nf90_put_att(ncid,NF90_GLOBAL,"Description","This is a grid of L-GAME."))
    
    ! defining the dimensions
    call nc_check(nf90_def_dim(ncid,"lon_model",nx,x_dimid))
    call nc_check(nf90_def_dim(ncid,"lat_model",ny,y_dimid))
    call nc_check(nf90_def_dim(ncid,"lon_model_plus1",nx+1,x_dimidp1))
    call nc_check(nf90_def_dim(ncid,"lat_model_plus1",ny+1,y_dimidp1))

    ! setting the dimension ID arrays
    ! 2D
    dimids(1) = y_dimid
    dimids(2) = x_dimid
    dimids_u(1) = y_dimid
    dimids_u(2) = x_dimidp1
    dimids_v(1) = y_dimidp1
    dimids_v(2) = x_dimid
    
    call nc_check(nf90_def_var(ncid,"lat_geo",NF90_REAL,dimids,varid_lat))
    call nc_check(nf90_put_att(ncid,varid_lat,"Description","latitude"))
    call nc_check(nf90_put_att(ncid,varid_lat,"Unit","radians"))
    
    call nc_check(nf90_def_var(ncid,"lon_geo",NF90_REAL,dimids,varid_lon))
    call nc_check(nf90_put_att(ncid,varid_lon,"Description","longitude"))
    call nc_check(nf90_put_att(ncid,varid_lon,"Unit","radians"))
    
    call nc_check(nf90_def_var(ncid,"lat_geo_u",NF90_REAL,dimids_u,varid_lat_u))
    call nc_check(nf90_put_att(ncid,varid_lat_u,"Description","latitude of u-vectors"))
    call nc_check(nf90_put_att(ncid,varid_lat_u,"Unit","radians"))
    
    call nc_check(nf90_def_var(ncid,"lon_geo_u",NF90_REAL,dimids_u,varid_lon_u))
    call nc_check(nf90_put_att(ncid,varid_lon_u,"Description","longitude of u-vectors"))
    call nc_check(nf90_put_att(ncid,varid_lon_u,"Unit","radians"))
    
    call nc_check(nf90_def_var(ncid,"lat_geo_v",NF90_REAL,dimids_v,varid_lat_v))
    call nc_check(nf90_put_att(ncid,varid_lat_v,"Description","latitude of v-vectors"))
    call nc_check(nf90_put_att(ncid,varid_lat_v,"Unit","radians"))
    
    call nc_check(nf90_def_var(ncid,"lon_geo_v",NF90_REAL,dimids_v,varid_lon_v))
    call nc_check(nf90_put_att(ncid,varid_lon_v,"Description","longitude of v-vectors"))
    call nc_check(nf90_put_att(ncid,varid_lon_v,"Unit","radians"))
    
    call nc_check(nf90_def_var(ncid,"land_fraction",NF90_REAL,dimids,varid_land_fraction))
    call nc_check(nf90_put_att(ncid,varid_land_fraction,"Description","land fraction"))
    call nc_check(nf90_put_att(ncid,varid_land_fraction,"Unit","1"))
    
    call nc_check(nf90_def_var(ncid,"lake_fraction",NF90_REAL,dimids,varid_lake_fraction))
    call nc_check(nf90_put_att(ncid,varid_lake_fraction,"Description","lake fraction"))
    call nc_check(nf90_put_att(ncid,varid_lake_fraction,"Unit","1"))
    
    call nc_check(nf90_def_var(ncid,"roughness_length",NF90_REAL,dimids,varid_roughness_length))
    call nc_check(nf90_put_att(ncid,varid_roughness_length,"Description","roughness length"))
    call nc_check(nf90_put_att(ncid,varid_roughness_length,"Unit","m"))
    
    call nc_check(nf90_def_var(ncid,"sfc_albedo",NF90_REAL,dimids,varid_sfc_albedo))
    call nc_check(nf90_put_att(ncid,varid_sfc_albedo,"Description","lake fraction"))
    call nc_check(nf90_put_att(ncid,varid_sfc_albedo,"Unit","1"))
    
    call nc_check(nf90_def_var(ncid,"sfc_rho_c",NF90_REAL,dimids,varid_sfc_rho_c))
    call nc_check(nf90_put_att(ncid,varid_sfc_rho_c,"Description","temperature conductivity"))
    call nc_check(nf90_put_att(ncid,varid_sfc_rho_c,"Unit","J/(K*m**3)"))
    
    call nc_check(nf90_def_var(ncid,"t_conduc_soil",NF90_REAL,dimids,varid_t_conduc_soil))
    call nc_check(nf90_put_att(ncid,varid_t_conduc_soil,"Description","lake fraction"))
    call nc_check(nf90_put_att(ncid,varid_t_conduc_soil,"Unit","m**2/2"))
    
    call nc_check(nf90_def_var(ncid,"oro",NF90_REAL,dimids,varid_oro))
    call nc_check(nf90_put_att(ncid,varid_oro,"Description","orography"))
    call nc_check(nf90_put_att(ncid,varid_oro,"Unit","m"))
    
    call nc_check(nf90_def_var(ncid,"oro_smoothed",NF90_REAL,dimids,varid_oro_smoothed))
    call nc_check(nf90_put_att(ncid,varid_oro_smoothed,"Description","smoothed orography"))
    call nc_check(nf90_put_att(ncid,varid_oro_smoothed,"Unit","m"))
    
    call nc_check(nf90_def_var(ncid,"u_dir",NF90_REAL,dimids_u,varid_dir_geo_u))
    call nc_check(nf90_put_att(ncid,varid_dir_geo_u,"Description","direction of u-vectors"))
    call nc_check(nf90_put_att(ncid,varid_dir_geo_u,"Unit","radians"))
    
    call nc_check(nf90_def_var(ncid,"v_dir",NF90_REAL,dimids_v,varid_dir_geo_v))
    call nc_check(nf90_put_att(ncid,varid_dir_geo_v,"Description","direction of v-vectors"))
    call nc_check(nf90_put_att(ncid,varid_dir_geo_v,"Unit","radians"))
    
    call nc_check(nf90_def_var(ncid,"u_dir_center",NF90_REAL,dimids,varid_dir_geo_u_scalar))
    call nc_check(nf90_put_att(ncid,varid_dir_geo_u_scalar,"Description","direction of u-vectors at the scalar data points"))
    call nc_check(nf90_put_att(ncid,varid_dir_geo_u_scalar,"Unit","radians"))
    
    ! ending the definition section
    call nc_check(nf90_enddef(ncid))
    
    ! writing the data to the netCDF file
    call nc_check(nf90_put_var(ncid,varid_lat,grid%lat_geo_scalar))
    call nc_check(nf90_put_var(ncid,varid_lon,grid%lon_geo_scalar))
    call nc_check(nf90_put_var(ncid,varid_lat_u,grid%lat_geo_u))
    call nc_check(nf90_put_var(ncid,varid_lon_u,grid%lon_geo_u))
    call nc_check(nf90_put_var(ncid,varid_lat_v,grid%lat_geo_v))
    call nc_check(nf90_put_var(ncid,varid_lon_v,grid%lon_geo_v))
    call nc_check(nf90_put_var(ncid,varid_land_fraction,grid%land_fraction))
    call nc_check(nf90_put_var(ncid,varid_lake_fraction,grid%lake_fraction))
    call nc_check(nf90_put_var(ncid,varid_roughness_length,grid%roughness_length))
    call nc_check(nf90_put_var(ncid,varid_sfc_albedo,grid%sfc_albedo))
    call nc_check(nf90_put_var(ncid,varid_sfc_rho_c,grid%sfc_rho_c))
    call nc_check(nf90_put_var(ncid,varid_t_conduc_soil,grid%t_conduc_soil))
    call nc_check(nf90_put_var(ncid,varid_oro,grid%oro))
    call nc_check(nf90_put_var(ncid,varid_oro_smoothed,grid%oro_smoothed))
    call nc_check(nf90_put_var(ncid,varid_dir_geo_u,grid%dir_geo_u))
    call nc_check(nf90_put_var(ncid,varid_dir_geo_v,grid%dir_geo_v))
    call nc_check(nf90_put_var(ncid,varid_dir_geo_u_scalar,grid%dir_geo_u_scalar))
    
    ! closing the netCDF file
    call nc_check(nf90_close(ncid))
    
  end subroutine write_grid
  
  subroutine read_geo(grid)
    
    ! This subroutine reads surface properties from a file.
    
    type(t_grid), intent(inout) :: grid ! grid properties
    
    ! local variables
    integer           :: ncid                   ! ID of the netCDF file
    character(len=64) :: filename               ! input filename
    integer           :: varid_land_fraction    ! variable ID of the land fraction
    integer           :: varid_lake_fraction    ! variable ID of the lake fraction
    integer           :: varid_roughness_length ! variable ID of the roughness length
    integer           :: varid_sfc_albedo       ! variable ID of the surface albedo
    integer           :: varid_sfc_rho_c        ! variable ID of the surface volumetric heap capacity
    integer           :: varid_t_conduc_soil    ! variable ID of the surface temperature conductivity
    integer           :: varid_oro              ! variable ID of the orography
    integer           :: varid_oro_smoothed     ! variable ID of the smoothed orography
    
    ! the filename of the grid file including the relative path
    filename = "../../grids/" // trim(grid_filename)
    
    ! opening the netCDF file
    call nc_check(nf90_open(trim(filename),NF90_CLOBBER,ncid))
    
    ! reading the variable IDs
    call nc_check(nf90_inq_varid(ncid,"land_fraction",varid_land_fraction))
    call nc_check(nf90_inq_varid(ncid,"lake_fraction",varid_lake_fraction))
    call nc_check(nf90_inq_varid(ncid,"roughness_length",varid_roughness_length))
    call nc_check(nf90_inq_varid(ncid,"sfc_albedo",varid_sfc_albedo))
    call nc_check(nf90_inq_varid(ncid,"sfc_rho_c",varid_sfc_rho_c))
    call nc_check(nf90_inq_varid(ncid,"t_conduc_soil",varid_t_conduc_soil))
    call nc_check(nf90_inq_varid(ncid,"oro",varid_oro))
    call nc_check(nf90_inq_varid(ncid,"oro_smoothed",varid_oro_smoothed))
    
    ! reading the arrays
    call nc_check(nf90_get_var(ncid,varid_land_fraction,grid%land_fraction))
    call nc_check(nf90_get_var(ncid,varid_lake_fraction,grid%lake_fraction))
    call nc_check(nf90_get_var(ncid,varid_roughness_length,grid%roughness_length))
    call nc_check(nf90_get_var(ncid,varid_sfc_albedo,grid%sfc_albedo))
    call nc_check(nf90_get_var(ncid,varid_sfc_rho_c,grid%sfc_rho_c))
    call nc_check(nf90_get_var(ncid,varid_t_conduc_soil,grid%t_conduc_soil))
    call nc_check(nf90_get_var(ncid,varid_oro,grid%oro))
    call nc_check(nf90_get_var(ncid,varid_oro_smoothed,grid%oro_smoothed))
    
    ! closing the netCDF file
    call nc_check(nf90_close(ncid))
    
  end subroutine read_geo

end module mo_read_write_grid








