! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_write_out

  ! This module handles everything dealing with IO.

  use mo_definitions,       only: t_state,wp,t_diag,t_grid
  use netcdf
  use mo_run_nml,           only: ny,nx,n_layers,scenario,run_id
  use mo_constituents_nml,  only: n_condensed_constituents,n_gaseous_constituents,n_constituents
  use mo_derived,           only: rel_humidity
  use mo_set_initial_state, only: bg_temp,bg_pres,geopot,nc_check

  implicit none
  
  contains
  
  subroutine write_output(state,diag,time_since_init_min,grid)
    
    ! This subroutine writes the state of the model atmosphere at a single time step to a NetCDF file.
    
    type(t_state), intent(in)    :: state               ! state to write out
    type(t_diag),  intent(inout) :: diag                ! diagnostic quantities
    integer,       intent(in)    :: time_since_init_min ! time in minutes since init
    type(t_grid),  intent(in)    :: grid                ! model grid
    
    ! local variables
    logical           :: lcontains_nan           ! used for detecting if the model crashed
    integer           :: ncid                    ! ID of the NetCDF file
    integer           :: x_dimid                 ! ID of the x-dimension
    integer           :: y_dimid                 ! ID of the y-dimension
    integer           :: z_dimid                 ! ID of the z-dimension
    integer           :: constituents_dimid      ! ID of the constituents dimension
    integer           :: single_double_dimid     ! ID of the dimension cotaining one double
    integer           :: dimids_2d(2)            ! dimensions of horizontal fields
    integer           :: dimids_3d(3)            ! dimensions of 3D fields
    integer           :: dimids_3d_rho(4)        ! dimensions of 3D mass density fields
    integer           :: varid_lat               ! variable ID of the lat coordinates
    integer           :: varid_lon               ! variable ID of the lon coordinates
    integer           :: varid_lat_center        ! variable ID of the latitude of the center
    integer           :: varid_lon_center        ! variable ID of the longitude of the center
    integer           :: varid_z                 ! variable ID of the z coordinates
    integer           :: varid_rho               ! variable ID of the 3D mass density fields
    integer           :: varid_r                 ! variable ID of the 3D relative humidity
    integer           :: varid_t                 ! variable ID of the 3D temperature field
    integer           :: varid_u                 ! variable ID of the 3D u wind field
    integer           :: varid_v                 ! variable ID of the 3D v wind field
    integer           :: varid_w                 ! variable ID of the 3D w wind field
    character(len=64) :: filename                ! output filename
    character(len=64) :: time_since_init_min_str ! time_since_init_min as string
    integer           :: ji,jk,jl                ! line indices
    real(wp)          :: upper_weight(ny,nx)     ! interpolation weights
    
    write(*,*) "Writing output ..."
    
    ! check if the model has crashed
    !$omp parallel workshare
    lcontains_nan = any(isnan(state%exner_pert))
    !$omp end parallel workshare
    if (lcontains_nan) then
      write(*,*) "Congratulations, the model crashed."
      call exit(1)
    endif
    
    ! creating the NetCDF file
    write(time_since_init_min_str,*) time_since_init_min
    time_since_init_min_str = adjustl(time_since_init_min_str)
    filename = trim(run_id) // "+" // trim(time_since_init_min_str) // "min.nc"
    call nc_check(nf90_create(trim(filename),NF90_CLOBBER,ncid))
    
    call nc_check(nf90_put_att(ncid,NF90_GLOBAL,"Description","This is output of L-GAME."))
    call nc_check(nf90_put_att(ncid,NF90_GLOBAL,"Run_ID",trim(run_id)))
    call nc_check(nf90_put_att(ncid,NF90_GLOBAL,"Time_since_initialization_in_minutes",trim(time_since_init_min_str)))
    
    ! defining the dimensions
    call nc_check(nf90_def_dim(ncid,"single_double",1,single_double_dimid))
    call nc_check(nf90_def_dim(ncid,"lon_model",nx,x_dimid))
    call nc_check(nf90_def_dim(ncid,"lat_model",ny,y_dimid))
    call nc_check(nf90_def_dim(ncid,"z",n_layers,z_dimid))
    call nc_check(nf90_def_dim(ncid,"jc",n_constituents,constituents_dimid))

    ! setting the dimension ID arrays
    ! 2D
    dimids_2d(1) = y_dimid
    dimids_2d(2) = x_dimid
    ! 3D
    dimids_3d(1) = y_dimid
    dimids_3d(2) = x_dimid
    dimids_3d(3) = z_dimid
    ! 3D mass density fields
    dimids_3d_rho(1) = y_dimid
    dimids_3d_rho(2) = x_dimid
    dimids_3d_rho(3) = z_dimid
    dimids_3d_rho(4) = constituents_dimid

    ! Define the variable. The type of the variable in this case is
    ! NF90_INT (4-byte integer).
    call nc_check(nf90_def_var(ncid,"lat_center",NF90_REAL,single_double_dimid,varid_lat_center))
    call nc_check(nf90_put_att(ncid,varid_lat_center, &
    "Description","latitude of the center of the model grid"))
    call nc_check(nf90_put_att(ncid,varid_lat_center,"Unit","rad"))
    
    call nc_check(nf90_def_var(ncid,"lon_center",NF90_REAL,single_double_dimid,varid_lon_center))
    call nc_check(nf90_put_att(ncid,varid_lon_center, &
    "Description","longitude of the model grid"))
    call nc_check(nf90_put_att(ncid,varid_lon_center,"Unit","rad"))
    
    call nc_check(nf90_def_var(ncid,"lat_model",NF90_REAL,y_dimid,varid_lat))
    call nc_check(nf90_put_att(ncid,varid_lat,"Unit","rad"))
    call nc_check(nf90_put_att(ncid,varid_lat, &
    "Description","latitudes of the gridpoints (in the frame of reference of the model)"))
    
    call nc_check(nf90_def_var(ncid,"lon_model",NF90_REAL,x_dimid,varid_lon))
    call nc_check(nf90_put_att(ncid,varid_lon, &
    "Description","longitudes of the gridpoints (in the frame of reference of the model)"))
    call nc_check(nf90_put_att(ncid,varid_lon,"Unit","rad"))
    
    call nc_check(nf90_def_var(ncid,"z",NF90_REAL,dimids_3d,varid_z))
    call nc_check(nf90_put_att(ncid,varid_z,"Description","z coordinates of the gridpoints (MSL)"))
    call nc_check(nf90_put_att(ncid,varid_z,"Unit","m"))
    
    call nc_check(nf90_def_var(ncid,"rho",NF90_REAL,dimids_3d_rho,varid_rho))
    call nc_check(nf90_put_att(ncid,varid_rho,"Description","mass densities"))
    call nc_check(nf90_put_att(ncid,varid_rho,"Unit","kg/m^3"))
    
    if (n_gaseous_constituents>1) then
      call nc_check(nf90_def_var(ncid,"r",NF90_REAL,dimids_3d,varid_r))
      call nc_check(nf90_put_att(ncid,varid_r,"Description","relative humidity"))
      call nc_check(nf90_put_att(ncid,varid_r,"Unit","%"))
    endif
    
    call nc_check(nf90_def_var(ncid,"T",NF90_REAL,dimids_3d,varid_t))
    call nc_check(nf90_put_att(ncid,varid_t,"Description","air temperature"))
    call nc_check(nf90_put_att(ncid,varid_t,"Unit","K"))
    
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
    ! latitude of the center
    call nc_check(nf90_put_var(ncid,varid_lat_center,grid%lat_center))
    ! longitude of the center
    call nc_check(nf90_put_var(ncid,varid_lon_center,grid%lon_center))
    ! latitude coordinates of the gridpoints
    call nc_check(nf90_put_var(ncid,varid_lat,grid%lat_scalar(:)))
    ! longitude coordinates of the gridpoints
    call nc_check(nf90_put_var(ncid,varid_lon,grid%lon_scalar(:)))
    ! z coordinates of the gridpoints
    call nc_check(nf90_put_var(ncid,varid_z,grid%z_scalar))
    
    ! writing the data to the output file
    ! 3D mass densities
    call nc_check(nf90_put_var(ncid,varid_rho,state%rho))
    
    ! 3D temperature
    !$omp parallel workshare
    diag%scalar_placeholder = (grid%theta_v_bg + state%theta_v_pert) &
    *(grid%exner_bg + state%exner_pert)
    !$omp end parallel workshare
    call nc_check(nf90_put_var(ncid,varid_t,diag%scalar_placeholder))
    
    if (n_gaseous_constituents>1) then
      ! relative humidity
      !$omp parallel do private(ji,jk,jl)
      do jl=1,n_layers
        do jk=1,nx
          do ji=1,ny
            diag%scalar_placeholder(ji,jk,jl) = rel_humidity(state%rho(ji,jk,jl,n_condensed_constituents+2), &
            diag%scalar_placeholder(ji,jk,jl))
          enddo
        enddo
      enddo
      !$omp end parallel do
      call nc_check(nf90_put_var(ncid,varid_r,diag%scalar_placeholder))
    endif
    
    ! 3D u wind
    !$omp parallel do private(jk)
    do jk=1,nx
      diag%scalar_placeholder(:,jk,:) = 0.5_wp*(state%wind_u(:,jk,:)+state%wind_u(:,jk+1,:))
    enddo
    !$omp end parallel do
    call nc_check(nf90_put_var(ncid,varid_u,diag%scalar_placeholder))
    
    ! 3D v wind
    !$omp parallel do private(ji)
    do ji=1,ny
      diag%scalar_placeholder(ji,:,:) = 0.5_wp*(state%wind_v(ji,:,:)+state%wind_v(ji+1,:,:))
    enddo
    !$omp end parallel do
    call nc_check(nf90_put_var(ncid,varid_v,diag%scalar_placeholder))
    
    ! 3D w wind
    !$omp parallel do private(jl,upper_weight)
    do jl=1,n_layers
      upper_weight = (grid%z_scalar(:,:,jl) -&
      grid%z_w(:,:,jl+1))/(grid%z_w(:,:,jl) - grid%z_w(:,:,jl+1))
      diag%scalar_placeholder(:,:,jl) = &
      upper_weight*state%wind_w(:,:,jl)+(1._wp-upper_weight)*state%wind_w(:,:,jl+1)
    enddo
    !$omp end parallel do
    call nc_check(nf90_put_var(ncid,varid_w,diag%scalar_placeholder))
    
    ! closing the NetCDF file
    call nc_check(nf90_close(ncid))
    
    write(*,*) "Output written."
    
  end subroutine write_output

end module mo_write_out








