! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module write_out

  ! This module handles everything dealing with IO.

  use definitions,       only: t_state,wp,t_diag,t_grid
  use netcdf
  use run_nml,           only: nlins,ncols,nlays,scenario,p_0,run_id
  use constituents_nml,  only: no_of_condensed_constituents
  use dictionary,        only: specific_gas_constants,spec_heat_capacities_p_gas
  use grid_generator,    only: bg_temp,bg_pres,geopot
  use set_initial_state, only: nc_check

  implicit none
  
  private
  
  public :: write_output
  
  contains
  
  subroutine write_output(state,diag,time_since_init_min,grid)
    ! reads out the state of the model atmosphere
    ! at a single timestep to a netcdf file
    
    type(t_state), intent(in)    :: state               ! state to write out
    type(t_diag),  intent(inout) :: diag                ! diagnostic quantities
    integer,       intent(in)    :: time_since_init_min ! time in minutes since init
    type(t_grid),  intent(in)    :: grid                ! model grid
    
    ! local variables
    integer                   :: ncid                      ! ID of the netcdf file
    integer                   :: x_dimid                   ! ID of the x dimension
    integer                   :: y_dimid                   ! ID of the y dimension
    integer                   :: z_dimid                   ! ID of the z dimension
    integer                   :: dimids_2d(2)              ! dimensions of horizontal fields
    integer                   :: dimids_3d(3)              ! dimensions of 3D fields
    integer                   :: varid_lat                 ! variable ID of the lat coordinates
    integer                   :: varid_lon                 ! variable ID of the lon coordinates
    integer                   :: varid_z                   ! variable ID of the z coordinates
    integer                   :: varid_p                   ! variable ID of the 3D pressure field
    integer                   :: varid_t                   ! variable ID of the 3D temperature field
    integer                   :: varid_u                   ! variable ID of the 3D u wind field
    integer                   :: varid_v                   ! variable ID of the 3D v wind field
    integer                   :: varid_w                   ! variable ID of the 3D w wind field
    character(len=64)         :: filename                  ! output filename
    character(len=64)         :: time_since_init_min_str   ! time_since_init_min as string
    integer                   :: ji,jk,jl                  ! line indices
    real(wp)                  :: upper_weight(nlins,ncols) ! interpolation weights
    
    ! creating the netcdf file
    write(time_since_init_min_str,*) time_since_init_min
    time_since_init_min_str = adjustl(time_since_init_min_str)
    filename = trim(run_id) // "+" // trim(time_since_init_min_str) // "min.nc"
    call nc_check(nf90_create(trim(filename),NF90_CLOBBER,ncid))
    
    call nc_check(nf90_put_att(ncid,NF90_GLOBAL,"Description","This is output of L-GAME."))
    call nc_check(nf90_put_att(ncid,NF90_GLOBAL,"Run_ID",trim(run_id)))
    call nc_check(nf90_put_att(ncid,NF90_GLOBAL,"Time_since_initialization_in_minutes",trim(time_since_init_min_str)))
    
    ! defining the dimensions
    call nc_check(nf90_def_dim(ncid,"lon_model",ncols,x_dimid))
    call nc_check(nf90_def_dim(ncid,"lat_model",nlins,y_dimid))
    call nc_check(nf90_def_dim(ncid,"z",nlays,z_dimid))

    ! setting the dimension ID arrays
    ! 2D
    dimids_2d(1) = y_dimid
    dimids_2d(2) = x_dimid
    ! 3D
    dimids_3d(1) = y_dimid
    dimids_3d(2) = x_dimid
    dimids_3d(3) = z_dimid

    ! Define the variable. The type of the variable in this case is
    ! NF90_INT (4-byte integer).
    call nc_check(nf90_def_var(ncid,"lat_model",NF90_REAL,y_dimid,varid_lat))
    call nc_check(nf90_put_att(ncid,varid_lat, &
    "Description","latitudes of the grid points (in the frame of reference of the model)"))
    call nc_check(nf90_put_att(ncid,varid_lat,"Unit","rad"))
    call nc_check(nf90_def_var(ncid,"lon_model",NF90_REAL,x_dimid,varid_lon))
    call nc_check(nf90_put_att(ncid,varid_lon, &
    "Description","longitudes of the grid points (in the frame of reference of the model)"))
    call nc_check(nf90_put_att(ncid,varid_lon,"Unit","rad"))
    call nc_check(nf90_def_var(ncid,"z",NF90_REAL,dimids_3d,varid_z))
    call nc_check(nf90_put_att(ncid,varid_z,"Description","z coordinates of the grid points (MSL)"))
    call nc_check(nf90_put_att(ncid,varid_z,"Unit","m"))
    call nc_check(nf90_def_var(ncid,"T",NF90_REAL,dimids_3d,varid_t))
    call nc_check(nf90_put_att(ncid,varid_t,"Description","air temperature"))
    call nc_check(nf90_put_att(ncid,varid_t,"Unit","K"))
    call nc_check(nf90_def_var(ncid,"p",NF90_REAL,dimids_3d,varid_p))
    call nc_check(nf90_put_att(ncid,varid_p,"Description","air pressure"))
    call nc_check(nf90_put_att(ncid,varid_p,"Unit","Pa"))
    call nc_check(nf90_def_var(ncid,"u",NF90_REAL,dimids_3d,varid_u))
    call nc_check(nf90_put_att(ncid,varid_u,"Description","zonal wind (in the frame of reference of the model)"))
    call nc_check(nf90_put_att(ncid,varid_u,"Unit","m/s"))
    call nc_check(nf90_def_var(ncid,"v",NF90_REAL,dimids_3d,varid_v))
    call nc_check(nf90_put_att(ncid,varid_v,"Description","meridional wind (in the frame of reference of the model)"))
    call nc_check(nf90_put_att(ncid,varid_v,"Unit","m/s"))
    call nc_check(nf90_def_var(ncid,"w",NF90_REAL,dimids_3d,varid_w))
    call nc_check(nf90_put_att(ncid,varid_w,"Description","vertical wind"))
    call nc_check(nf90_put_att(ncid,varid_w,"Unit","m/s"))
    
    ! ending the definition section
    call nc_check(nf90_enddef(ncid))
    ! latitude coordinates of the grid points
    call nc_check(nf90_put_var(ncid,varid_lat,grid%lat_scalar(2:nlins+1)))
    ! longitude coordinates of the grid points
    call nc_check(nf90_put_var(ncid,varid_lon,grid%lon_scalar(2:ncols+1)))
    ! z coordinates of the grid points
    call nc_check(nf90_put_var(ncid,varid_z,grid%z_geo_scal(2:nlins+1,2:ncols+1,:)))
    
    ! 3D temperature
    diag%scalar_placeholder(2:nlins+1,2:ncols+1,:) =  (grid%theta_bg(2:nlins+1,2:ncols+1,:) &
    + state%theta_pert(2:nlins+1,2:ncols+1,:)) &
    *(grid%exner_bg(2:nlins+1,2:ncols+1,:) + state%exner_pert(2:nlins+1,2:ncols+1,:))
    call nc_check(nf90_put_var(ncid,varid_t,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
    
    ! writing the data to the output file
    ! 3D pressure
    diag%scalar_placeholder(2:nlins+1,2:ncols+1,:) = state%rho(2:nlins+1,2:ncols+1,:,no_of_condensed_constituents+1) &
    *specific_gas_constants(0)*diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)
    call nc_check(nf90_put_var(ncid,varid_p,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
    
    ! 3D u wind
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=1,ncols
      diag%scalar_placeholder(2:nlins+1,jk+1,:) = 0.5_wp*(state%wind_u(2:nlins+1,jk,:)+state%wind_u(2:nlins+1,jk+1,:))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call nc_check(nf90_put_var(ncid,varid_u,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
    
    ! 3D v wind
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji)
    do ji=1,nlins
      diag%scalar_placeholder(ji+1,2:ncols+1,:) = 0.5_wp*(state%wind_v(ji,2:ncols+1,:)+state%wind_v(ji+1,2:ncols+1,:))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call nc_check(nf90_put_var(ncid,varid_v,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
    
    ! 3D w wind
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl,upper_weight)
    do jl=1,nlays
      upper_weight(:,:) = (grid%z_geo_scal(2:nlins+1,2:ncols+1,jl) -&
      grid%z_geo_w(2:nlins+1,2:ncols+1,jl+1))/(grid%z_geo_w(2:nlins+1,2:ncols+1,jl)  &
      - grid%z_geo_w(2:nlins+1,2:ncols+1,jl+1))
      diag%scalar_placeholder(2:nlins+1,2:ncols+1,jl) = &
      upper_weight(:,:)*state%wind_w(2:nlins+1,2:ncols+1,jl)+(1._wp-upper_weight(:,:))*state%wind_w(2:nlins+1,2:ncols+1,jl+1)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    call nc_check(nf90_put_var(ncid,varid_w,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
    
    ! closing the netcdf file
    call nc_check(nf90_close(ncid))
    
  end subroutine write_output

end module write_out







