! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

! This file contains the calculation of the grid properties.

module mo_grid_generator

  use netcdf
  use mo_definitions,        only: wp,t_grid
  use mo_run_nml,            only: ny,nx,n_layers,n_levels,dy,dx,toa,n_oro_layers,stretching_parameter,scenario,lat_center, &
                                   lon_center,lplane,n_flat_layers,eff_hor_res
  use mo_diff_nml,           only: klemp_begin_rel
  use mo_constants,          only: r_e,rho_h2o,t_0,M_PI,p_0,omega,gravity,p_0_standard, &
                                   surface_temp,tropo_height,inv_height,t_grad_inv,lapse_rate, &
                                   r_d,c_d_p,epsilon_security
  use mo_surface_nml,        only: nsoillays,orography_id,lsleve
  use mo_gradient_operators, only: grad_hor_cov,grad_hor,grad_vert
  use mo_io_nml,             only: lwrite_grid,lread_geo,oro_raw_filename
  use mo_read_write_grid,    only: write_grid,read_geo
  use mo_set_initial_state,  only: bg_temp,bg_pres,geopot,nc_check
  use mo_bc_nml,             only: lperiodic

  implicit none
  
  integer  :: n_damping_levels ! number of levels in which the swamp layer is active
  
  contains
  
  subroutine grid_setup(grid)
    
    ! This module computes the grid quantities.
    
    type(t_grid), intent(inout) :: grid ! the model grid
    
    ! local variables
    integer                       :: ext_fileunit                  ! file unit of an external data file
    integer                       :: ncid                          ! netCDF file ID
    integer                       :: nlon_ext                      ! number of longitude points of the external data grid
    integer                       :: nlat_ext                      ! number of latitude points of the external data grid
    integer                       :: nlon_ext_sst                  ! number of longitude points of the SST grid
    integer                       :: nlat_ext_sst                  ! number of latitude points of the SST grid
    integer                       :: lat_index_ext                 ! latitude index of a grid point of the external data
    integer                       :: lon_index_ext                 ! longitude index of a grid point of the external data
    integer                       :: left_index_ext                ! helper variable for interpolating external data to the GAME grid
    integer                       :: right_index_ext               ! helper variable for interpolating external data to the GAME grid
    integer                       :: lower_index_ext               ! helper variable for interpolating external data to the GAME grid
    integer                       :: upper_index_ext               ! helper variable for interpolating external data to the GAME grid
    integer                       :: n_points_ext_domain           ! helper variable for interpolating external data to the GAME grid
    integer                       :: jm_used                       ! helper variable for interpolating external data to the GAME grid
    integer                       :: jn_used                       ! helper variable for interpolating external data to the GAME grid
    integer                       :: etopo_oro_id                  ! variable ID of the input orography
    integer                       :: lsmask_id                     ! netCDF ID of the land-sea mask of the NCEP NSST grid
    integer                       :: ghcn_cams_id                  ! netCDF ID of the input 2-m-temperature mean (from GHCN-CAMS)
    integer(2),       allocatable :: lake_depth_gldb_raw(:,:)      ! GLDB lake depth data as read from file
    character(len=1), allocatable :: glcc_raw(:,:)                 ! GLCC raw data
    integer,          allocatable :: glcc(:,:)                     ! GLCC data
    integer                       :: ji                            ! horizontal index
    integer                       :: jk                            ! horizontal index
    integer                       :: jl                            ! layer or level index
    integer                       :: jm                            ! helper index
    integer                       :: jn                            ! helper index
    real(wp)                      :: lat_left_upper                ! latitude coordinate of upper left corner
    real(wp)                      :: lon_left_upper                ! longitude coordinate of upper left corner
    real(wp)                      :: dlat                          ! mesh size in y direction as angle
    real(wp)                      :: dlon                          ! mesh size in x direction as angle
    real(wp)                      :: A                             ! variable for calculating the vertical grid (height of a level without orography)
    real(wp)                      :: B                             ! variable for calculating the vertical grid (orography weighting factor)
    real(wp)                      :: sigma_z                       ! variable for calculating the vertical grid (A/toa)
    real(wp)                      :: z_rel                         ! variable for calculating the vertical grid (z/toa with equidistant levels)
    real(wp)                      :: vertical_vector_pre(n_levels) ! heights of the levels in a column
    real(wp)                      :: base_area                     ! variable for calculating the vertical grid
    real(wp)                      :: lower_z,upper_z,lower_length  ! variables needed for area calculations
    real(wp)                      :: height_mountain               ! height of Gaussian mountain (needed for test case)
    real(wp)                      :: x_coord                       ! help variable needed for the Schär test
    real(wp)                      :: rescale_factor                ! soil grid rescaling factor
    real(wp)                      :: sigma_soil                    ! sigma of the soil grid
    real(wp)                      :: density_soil                  ! typical density of soil
    real(wp)                      :: c_p_soil                      ! typical c_p of soil
    real(wp)                      :: c_p_water                     ! typical c_p of water
    real(wp)                      :: albedo_water                  ! albedo of water
    real(wp)                      :: albedo_soil                   ! albedo of soil
    real(wp)                      :: albedo_ice                    ! albedo of ice
    real(wp)                      :: t_conductivity_water          ! temperature conductivity of water
    real(wp)                      :: t_conductivity_soil           ! temperature conductivity of soil
    real(wp)                      :: lat_lower_center              ! variable for calculating the TRSK weights
    real(wp)                      :: lat_upper_center              ! variable for calculating the TRSK weights
    real(wp)                      :: rot_y(3,3)                    ! rotation matrix around the global y-axis
    real(wp)                      :: rot_z(3,3)                    ! rotation matrix around the global z-axis
    real(wp)                      :: rot(3,3)                      ! complete rotation matrix
    real(wp)                      :: r_old(3)                      ! positional vector before rotation
    real(wp)                      :: r_new(3)                      ! positional vector after rotation
    real(wp)                      :: basis_old(3)                  ! old local basis vector
    real(wp)                      :: basis_new(3)                  ! new local basis vector
    real(wp)                      :: local_i(3)                    ! local i-vector
    real(wp)                      :: local_j(3)                    ! local j-vector
    real(wp)                      :: x_basis_local,y_basis_local   ! local Cartesian components of the local basis vector
    real(wp)                      :: lat_local,lon_local,max_z     ! helper variables    
    real(wp)                      :: lon_geo_scalar_used           ! helper variable for interpolating external data to the model grid
    real(wp), allocatable         :: lake_depth_gldb(:,:)          ! GLDB lake depth data
    integer,  allocatable         :: etopo_oro(:,:)                ! ETOPO orography
    integer,  allocatable         :: invalid_counter(:,:)          ! counts invalid values encountered in an interpolation
    integer,  allocatable         :: ncep_nsst_lsmask(:,:)         ! NCEP NSST land-sea mask
    real(wp), allocatable         :: ghcn_cams(:,:,:)              ! GHCN-CAMS data (2-m-temperature mean)
    real(wp)                      :: toa_oro                       ! top of terrain-following coordinates
    real(wp)                      :: delta_lat_ext                 ! latitude resolution of the external data grid
    real(wp)                      :: delta_lon_ext                 ! longitude resolution of the external data grid
    real(wp)                      :: delta_lat_ext_sst             ! latitude resolution of the SST grid
    real(wp)                      :: delta_lon_ext_sst             ! longitude resolution of the SST grid
    real(wp)                      :: fractions_sum                 ! sum of land fraction and lake fraction
    real(wp)                      :: max_oro                       ! maximum orography value
    real(wp)                      :: dq_value                      ! data quality value
    
    ! Horizontal grid properties
    ! --------------------------
    
    ! setting the center and direction of the grid
    grid%lat_center = lat_center
    grid%lon_center = lon_center
    
    ! setting the latitude and longitude coordinates of the scalar gridpoints
    ! setting the dy of the model grid
    dlat = dy/r_e
    dlon = dx/r_e
    lat_left_upper = (ny-1._wp)/2._wp*dlat
    lon_left_upper = -(nx-1._wp)/2._wp*dlon
    !$omp parallel do private(ji)
    do ji=1,ny
      grid%lat_scalar(ji) = lat_left_upper - dlat*(ji-1._wp)
      ! this will be modified later
      grid%lat_geo_scalar(ji,:) = grid%lat_scalar(ji)
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(jk)
    do jk=1,nx
      grid%lon_scalar(jk) = lon_left_upper + dlon*(jk-1._wp)
      ! this will be modified later
      grid%lon_geo_scalar(:,jk) = grid%lon_scalar(jk)
    enddo
    !$omp end parallel do
    
    ! this will be modified later
    !$omp parallel do private(ji,jk)
    do jk=1,nx+1
      do ji=1,ny
        grid%lat_geo_u(ji,jk) = grid%lat_scalar(ji)
        if (jk==nx+1) then
          grid%lon_geo_u(ji,jk) = grid%lon_scalar(jk-1) + 0.5_wp*dlon
        else
          grid%lon_geo_u(ji,jk) = grid%lon_scalar(jk) - 0.5_wp*dlon
        endif
      enddo
    enddo
    !$omp end parallel do
    
    ! this will be modified later
    !$omp parallel do private(ji,jk)
    do jk=1,nx
      do ji=1,ny+1
        if (ji==ny+1) then
          grid%lat_geo_v(ji,jk) = grid%lat_scalar(ji-1) - 0.5_wp*dlat
        else
          grid%lat_geo_v(ji,jk) = grid%lat_scalar(ji) + 0.5_wp*dlat
        endif
        grid%lon_geo_v(ji,jk) =  grid%lon_scalar(jk)
      enddo
    enddo
    !$omp end parallel do
    
    ! this will be modified later
    !$omp parallel workshare
    grid%dir_geo_u(:,:) = 0._wp
    !$omp end parallel workshare
    
    ! this will be modified later
    !$omp parallel workshare
    grid%dir_geo_v(:,:) = 0.5_wp*M_PI
    !$omp end parallel workshare
    
    ! this will be modified later
    !$omp parallel workshare
    grid%dir_geo_u_scalar(:,:) = 0._wp
    !$omp end parallel workshare
    
    ! performing the rotation
    ! setting up the rotation matrix
    rot_y(1,1) = cos(lat_center)
    rot_y(1,2) = 0._wp
    rot_y(1,3) = -sin(lat_center)
    rot_y(2,1) = 0._wp
    rot_y(2,2) = 1._wp
    rot_y(2,3) = 0._wp
    rot_y(3,1) = sin(lat_center)
    rot_y(3,2) = 0._wp
    rot_y(3,3) = cos(lat_center)
    rot_z(1,1) = cos(lon_center)
    rot_z(1,2) = -sin(lon_center)
    rot_z(1,3) = 0._wp
    rot_z(2,1) = sin(lon_center)
    rot_z(2,2) = cos(lon_center)
    rot_z(2,3) = 0._wp
    rot_z(3,1) = 0._wp
    rot_z(3,2) = 0._wp
    rot_z(3,3) = 1._wp
    rot = matmul(rot_z,rot_y)
    
    ! calculating the geographic coordinates of the gridpoints
    ! scalar points
    !$omp parallel do private(ji,jk,r_old,r_new,basis_old,basis_new,local_i,local_j,x_basis_local,y_basis_local)
    do jk=1,nx
      do ji=1,ny
        ! calculating the local i-vector before the rotation
        call calc_local_i(grid%lon_geo_scalar(ji,jk),basis_old)
        call find_global_normal(grid%lat_geo_scalar(ji,jk),grid%lon_geo_scalar(ji,jk),r_old)
        r_new = matmul(rot,r_old)
        call find_geos(r_new,grid%lat_geo_scalar(ji,jk),grid%lon_geo_scalar(ji,jk))
        ! calculating the local basis elements
        call calc_local_i(grid%lon_geo_scalar(ji,jk),local_i)
        call calc_local_j(grid%lat_geo_scalar(ji,jk),grid%lon_geo_scalar(ji,jk),local_j)
        ! rotating the basis vector
        basis_new = matmul(rot,basis_old)
        ! calculting th components of the local basis vector
        x_basis_local = dot_product(local_i,basis_new)
        y_basis_local = dot_product(local_j,basis_new)
        ! calculating the direction of this vector
        grid%dir_geo_u_scalar(ji,jk) = atan2(y_basis_local,x_basis_local)
      enddo
    enddo
    !$omp end parallel do
    
    ! u-vector points, including directions
    !$omp parallel do private(ji,jk,r_old,r_new,basis_old,basis_new,local_i,local_j,x_basis_local,y_basis_local)
    do jk=1,nx+1
      do ji=1,ny
        ! calculating the local i-vector before the rotation
        call calc_local_i(grid%lon_geo_u(ji,jk),basis_old)
        ! rotating the gridpoint itself
        call find_global_normal(grid%lat_geo_u(ji,jk),grid%lon_geo_u(ji,jk),r_old)
        r_new = matmul(rot,r_old)
        call find_geos(r_new,grid%lat_geo_u(ji,jk),grid%lon_geo_u(ji,jk))
        ! calculating the local basis elements
        call calc_local_i(grid%lon_geo_u(ji,jk),local_i)
        call calc_local_j(grid%lat_geo_u(ji,jk),grid%lon_geo_u(ji,jk),local_j)
        ! rotating the basis vector
        basis_new = matmul(rot,basis_old)
        ! calculting th components of the local basis vector
        x_basis_local = dot_product(local_i,basis_new)
        y_basis_local = dot_product(local_j,basis_new)
        ! calculating the direction of this vector
        grid%dir_geo_u(ji,jk) = atan2(y_basis_local,x_basis_local)
      enddo
    enddo
    !$omp end parallel do
    
    ! v-vector points, including directions
    !$omp parallel do private(ji,jk,r_old,r_new,basis_old,basis_new,local_i,local_j,x_basis_local,y_basis_local)
    do jk=1,nx
      do ji=1,ny+1
        ! calculating the local j-vector before the rotation
        call calc_local_j(grid%lat_geo_v(ji,jk),grid%lon_geo_v(ji,jk),basis_old)
        ! rotating the gridpoint itself
        call find_global_normal(grid%lat_geo_v(ji,jk),grid%lon_geo_v(ji,jk),r_old)
        r_new = matmul(rot,r_old)
        call find_geos(r_new,grid%lat_geo_v(ji,jk),grid%lon_geo_v(ji,jk))
        ! calculating the local basis elements
        call calc_local_i(grid%lon_geo_v(ji,jk),local_i)
        call calc_local_j(grid%lat_geo_v(ji,jk),grid%lon_geo_v(ji,jk),local_j)
        ! rotating the basis vector
        basis_new = matmul(rot,basis_old)
        ! calculting th components of the local basis vector
        x_basis_local = dot_product(local_i,basis_new)
        y_basis_local = dot_product(local_j,basis_new)
        ! calculating the direction of this vector
        grid%dir_geo_v(ji,jk) = atan2(y_basis_local,x_basis_local)
      enddo
    enddo
    !$omp end parallel do
    
    ! setting the Coriolis vector
    !$omp parallel do private(ji,jk,y_basis_local)
    do jk=1,nx
      do ji=1,ny+1
        y_basis_local = sin(grid%dir_geo_v(ji,jk) - 0.5_wp*M_PI)
        grid%fvec_x(ji,jk) = 2._wp*omega*cos(grid%lat_geo_v(ji,jk))*y_basis_local
      enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(ji,jk,y_basis_local)
    do jk=1,nx+1
      do ji=1,ny
        y_basis_local = sin(grid%dir_geo_u(ji,jk) + 0.5_wp*M_PI)
        grid%fvec_y(ji,jk) = 2._wp*omega*cos(grid%lat_geo_u(ji,jk))*y_basis_local
      enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(ji,jk,lat_local,lon_local,r_old,r_new)
    do jk=1,nx+1
      do ji=1,ny+1
        if (ji==ny+1) then
          lat_local = grid%lat_scalar(ji-1) - 0.5*dlat
        else
          lat_local = grid%lat_scalar(ji) + 0.5*dlat
        endif
        if (jk==nx+1) then
          lon_local = grid%lon_scalar(jk-1) + 0.5*dlon
        else
          lon_local = grid%lon_scalar(jk) - 0.5*dlon
        endif
        call find_global_normal(lat_local,lon_local,r_old)
        r_new = matmul(rot,r_old)
        call find_geos(r_new,lat_local,lon_local)
        grid%fvec_z(ji,jk) = 2._wp*omega*sin(lat_local)
      enddo
    enddo
    !$omp end parallel do
    
    ! Physical grid properties
    ! ------------------------
    
    select case (orography_id)
    
      ! no orography
      case(0)
        
        !$omp parallel workshare
        grid%oro = 0._wp
        ! mean surface temperature for an Earth without real orography
        grid%t_const_soil = t_0 + 25._wp*cos(2._wp*grid%lat_geo_scalar)
        !$omp end parallel workshare
      
      ! real orography (and other surface properties)
      case(1)
        
        if (lread_geo) then
          call read_geo(grid)
        else
          
          ! Land fraction
          ! -------------
          
          write(*,*) "Setting the land fraction ..."
          
          nlat_ext = 21600
          nlon_ext = 43200
          
          ! opening the GLCC file
          open(action="read",file="../../grids/phys_sfc_quantities/sfc-fields-usgs-veg30susgs",form="unformatted", &
          access="direct",recl=nlon_ext,newunit=ext_fileunit)
          
          allocate(glcc_raw(nlat_ext,nlon_ext))
          allocate(glcc(nlat_ext,nlon_ext))
          
          !$omp parallel do private(ji)
          do ji=1,nlat_ext
            read(unit=ext_fileunit,rec=ji) glcc_raw(ji,:)
            glcc(ji,:) = ichar(glcc_raw(ji,:))
          end do
          !$omp end parallel do
           
          deallocate(glcc_raw)
          
          delta_lat_ext = M_PI/nlat_ext
          delta_lon_ext = 2._wp*M_PI/nlon_ext
          
          !$omp parallel do private(ji,jk,lat_index_ext,lon_index_ext,left_index_ext,right_index_ext, &
          !$omp lower_index_ext,upper_index_ext,n_points_ext_domain,jm_used,jn_used)
          do jk=1,nx
            do ji=1,ny
              
              ! computing the indices of the GLCC grid point that is the closest to the center of this grid cell
              lat_index_ext = nlat_ext/2 - int(grid%lat_geo_scalar(ji,jk)/delta_lat_ext)
              lon_index_ext = nlon_ext/2 + int(grid%lon_geo_scalar(ji,jk)/delta_lon_ext)
              
              ! computing the subset of the external data to use for the interpolation
              call calc_ext_subset(nlat_ext,nlon_ext,grid%lat_geo_scalar(ji,jk),lat_index_ext,lon_index_ext,n_points_ext_domain, &
                                   lower_index_ext,upper_index_ext,left_index_ext,right_index_ext)
              
              ! looping over all points of the input dataset in the vicinity of the grid cell at hand
              do jm=upper_index_ext,lower_index_ext
                do jn=left_index_ext,right_index_ext
                  
                  ! correcting the indices for boundary cases
                  call correct_ext_data_indices(jm,jn,nlat_ext,nlon_ext,jm_used,jn_used)
                  
                  if (glcc(jm_used,jn_used)/=16) then
                    grid%land_fraction(ji,jk) = grid%land_fraction(ji,jk)+1._wp
                  endif
                  
                enddo
              enddo
              
              grid%land_fraction(ji,jk) = grid%land_fraction(ji,jk)/n_points_ext_domain
              
            enddo
          enddo
          !$omp end parallel do
          
          deallocate(glcc)
          
          write(*,*) "Land fraction set."
          
          ! Lake fraction
          ! -------------
          
          write(*,*) "Setting the lake fraction ..."
          
          nlat_ext = 21600
          nlon_ext = 43200
            
          ! opening the lake depth file
          open(action="read",file="../../grids/phys_sfc_quantities/GlobalLakeDepth.dat",form="unformatted", &
          access="direct",recl=2*nlon_ext,newunit=ext_fileunit)
          
          allocate(lake_depth_gldb_raw(nlat_ext,nlon_ext))
          allocate(lake_depth_gldb(nlat_ext,nlon_ext))
          
          !$omp parallel do private(ji)
          do ji=1,nlat_ext
            read(unit=ext_fileunit,rec=ji) lake_depth_gldb_raw(ji,:)
            lake_depth_gldb(ji,:) = lake_depth_gldb_raw(ji,:)/10._wp
          enddo
          !$omp end parallel do
          
          ! closing the lake depth file
          close(ext_fileunit)
          
          deallocate(lake_depth_gldb_raw)
          
          delta_lat_ext = M_PI/nlat_ext
          delta_lon_ext = 2._wp*M_PI/nlon_ext
          
          ! reading the NCEP NSST land-sea mask
          nlat_ext_sst = 180
          nlon_ext_sst = 360
          allocate(ncep_nsst_lsmask(nlon_ext_sst,nlat_ext_sst))
          
          call nc_check(nf90_open(trim("../../grids/phys_sfc_quantities/lsmask.nc"),NF90_CLOBBER,ncid))
          call nc_check(nf90_inq_varid(ncid,"mask",lsmask_id))
          call nc_check(nf90_get_var(ncid,lsmask_id,ncep_nsst_lsmask))
          call nc_check(nf90_close(ncid))
          
          ! calculating the properties of the NCEP NSST land-sea mask grid
          delta_lat_ext_sst = M_PI/nlat_ext_sst
          delta_lon_ext_sst = 2._wp*M_PI/nlon_ext_sst
          
          !$omp parallel do private(ji,jk,lat_index_ext,lon_index_ext,left_index_ext,right_index_ext, &
          !$omp lower_index_ext,upper_index_ext,n_points_ext_domain,jm_used,jn_used,lon_geo_scalar_used)
          do jk=1,nx
            do ji=1,ny
              
              ! if there is no land in this grid cell, there can also be no lakes in this grid cell
              if (grid%land_fraction(ji,jk)==0._wp) then
                cycle
              endif
              
              ! computing the indices of the GLDB grid point that is the closest to the center of this grid cell
              lat_index_ext = nlat_ext/2 - int(grid%lat_geo_scalar(ji,jk)/delta_lat_ext)
              lon_index_ext = nlon_ext/2 + int(grid%lon_geo_scalar(ji,jk)/delta_lon_ext)
              
              ! computing the subset of the external data to use for the interpolation
              call calc_ext_subset(nlat_ext,nlon_ext,grid%lat_geo_scalar(ji,jk),lat_index_ext,lon_index_ext,n_points_ext_domain, &
                                   lower_index_ext,upper_index_ext,left_index_ext,right_index_ext)
              
              ! looping over all points of the input dataset in the vicinity of the grid cell at hand
              do jm=upper_index_ext,lower_index_ext
                do jn=left_index_ext,right_index_ext
                
                  ! correcting the indices for boundary cases
                  call correct_ext_data_indices(jm,jn,nlat_ext,nlon_ext,jm_used,jn_used)
                  
                  if (lake_depth_gldb(jm_used,jn_used)>0._wp) then
                    
                    lon_geo_scalar_used = grid%lon_geo_scalar(ji,jk)
                    if (lon_geo_scalar_used<0._wp) then
                      lon_geo_scalar_used = lon_geo_scalar_used+2._wp*M_PI
                    endif
                    lon_index_ext = int(lon_geo_scalar_used/delta_lon_ext_sst)
                    
                    ! a lake is only considered a lake if it is not part of the NCEP NSST grid
                    lat_index_ext = nlat_ext_sst/2 - int(grid%lat_geo_scalar(ji,jk)/delta_lat_ext_sst)
                    lon_index_ext = int(lon_geo_scalar_used/delta_lon_ext_sst)
                    lat_index_ext = max(1,lat_index_ext)
                    lat_index_ext = min(nlat_ext_sst,lat_index_ext)
                    lon_index_ext = max(1,lon_index_ext)
                    lon_index_ext = min(nlon_ext_sst,lon_index_ext)
                    
                    if (ncep_nsst_lsmask(lon_index_ext,lat_index_ext)==0) then
                      grid%lake_fraction(ji,jk) = grid%lake_fraction(ji,jk)+1._wp
                    endif
                    
                  endif
                    
                enddo
              enddo
              
              grid%lake_fraction(ji,jk) = grid%lake_fraction(ji,jk)/n_points_ext_domain
              
            enddo
          enddo
          !$omp end parallel do
          
          deallocate(lake_depth_gldb)
          deallocate(ncep_nsst_lsmask)
          
          ! restricting the sum of lake fraction and land fraction to one
          !$omp parallel do private(ji,jk,fractions_sum)
          do jk=1,nx
            do ji=1,ny
              
              if (grid%land_fraction(ji,jk)+grid%lake_fraction(ji,jk)>1._wp) then
                fractions_sum = grid%land_fraction(ji,jk)+grid%lake_fraction(ji,jk)
                grid%land_fraction(ji,jk) = grid%land_fraction(ji,jk)/fractions_sum
                grid%lake_fraction(ji,jk) = grid%lake_fraction(ji,jk)/fractions_sum
              endif
            
            enddo
          enddo
          !$omp end parallel do
          
          write(*,*) "Lake fraction set."
          
          ! Orography
          
          write(*,*) "Setting the orography ..."
          
          nlat_ext = 10801
          nlon_ext = 21601
          
          ! opening the netCDF file
          call nc_check(nf90_open(trim("../../grids/phys_sfc_quantities/" // trim(oro_raw_filename)),NF90_CLOBBER,ncid))
          
          ! reading the variable IDs
          call nc_check(nf90_inq_varid(ncid,"z",etopo_oro_id))
          
          ! allocating memory for reading the grid files
          allocate(etopo_oro(nlon_ext,nlat_ext))
          
          ! reading the arrays
          call nc_check(nf90_get_var(ncid,etopo_oro_id,etopo_oro))
          
          ! closing the netCDF file
          call nc_check(nf90_close(ncid))
      
          ! setting the unfiltered orography
          
          delta_lat_ext = M_PI/nlat_ext
          delta_lon_ext = 2._wp*M_PI/nlon_ext
          
          !$omp parallel do private(ji,jk,lat_index_ext,lon_index_ext,left_index_ext,right_index_ext, &
          !$omp lower_index_ext,upper_index_ext,n_points_ext_domain,jm_used,jn_used)
          do jk=1,nx
            do ji=1,ny
              
              ! if there is primarily sea in this grid cell, the orography remains at zero
              if (grid%land_fraction(ji,jk)+grid%lake_fraction(ji,jk)<0.5_wp) then
                cycle
              endif
              
              ! computing the indices of the GLDB grid point that is the closest to the center of this grid cell
              lat_index_ext = nlat_ext/2 + int(grid%lat_geo_scalar(ji,jk)/delta_lat_ext)
              lon_index_ext = nlon_ext/2 + int(grid%lon_geo_scalar(ji,jk)/delta_lon_ext)
              
              ! computing the subset of the external data to use for the interpolation
              call calc_ext_subset(nlat_ext,nlon_ext,grid%lat_geo_scalar(ji,jk),lat_index_ext,lon_index_ext,n_points_ext_domain, &
                                   lower_index_ext,upper_index_ext,left_index_ext,right_index_ext)
              
              ! looping over all points of the input dataset in the vicinity of the grid cell at hand
              do jm=upper_index_ext,lower_index_ext
                do jn=left_index_ext,right_index_ext
                  
                  ! correcting the indices for boundary cases
                  call correct_ext_data_indices(jm,jn,nlat_ext,nlon_ext,jm_used,jn_used)
                  
                  ! adding the orography value, restrictued to the global minimum of the orography
                  grid%oro(ji,jk) = grid%oro(ji,jk)+max(etopo_oro(jn_used,jm_used),-440)
                  
                enddo
              enddo
              
              grid%oro(ji,jk) = grid%oro(ji,jk)/n_points_ext_domain
              
            enddo
          enddo
          !$omp end parallel do
          
          ! freeing the memory
          deallocate(etopo_oro)
          
          !$omp parallel workshare
          grid%oro_smoothed = grid%oro
          !$omp end parallel workshare
          call smooth_hor_scalar(grid%oro_smoothed)
          if (lsleve) then
            !$omp parallel workshare
            grid%oro = grid%oro_smoothed + 0.2_wp*(grid%oro - grid%oro_smoothed)
            !$omp end parallel workshare
          else
            !$omp parallel workshare
            grid%oro = grid%oro_smoothed
            !$omp end parallel workshare
          endif
          
          write(*,*) "Orography set."
          
          ! Lower boundary soil temperature
          ! -------------------------------
          
          write(*,*) "Setting the lower boundary soil temperature ..."
          
          nlat_ext = 360
          nlon_ext = 720
          
          allocate(ghcn_cams(nlon_ext,nlat_ext,12))
          
          ! reading the GHCN-CAMS data
          call nc_check(nf90_open(trim("../../grids/phys_sfc_quantities/air.mon.ltm.nc"),NF90_CLOBBER,ncid))
          call nc_check(nf90_inq_varid(ncid,"air",ghcn_cams_id))
          call nc_check(nf90_get_var(ncid,ghcn_cams_id,ghcn_cams))
          call nc_check(nf90_close(ncid))
          
          delta_lat_ext = M_PI/nlat_ext
          delta_lon_ext = 2._wp*M_PI/nlon_ext
          
          allocate(invalid_counter(ny,nx))
          
          !$omp parallel workshare
          invalid_counter = 0
          !$omp end parallel workshare
          
          !$omp parallel do private(ji,jk,lat_index_ext,lon_index_ext,left_index_ext,right_index_ext, &
          !$omp lower_index_ext,upper_index_ext,n_points_ext_domain,jm_used,jn_used,lon_geo_scalar_used)
          do jk=1,nx
            do ji=1,ny
              
              ! computing the indices of the GLCC grid point that is the closest to the center of this grid cell
              lat_index_ext = nlat_ext/2 - int(grid%lat_geo_scalar(ji,jk)/delta_lat_ext)
              lon_geo_scalar_used = grid%lon_geo_scalar(ji,jk)
              if (lon_geo_scalar_used<0._wp) then
                lon_geo_scalar_used = lon_geo_scalar_used+2._wp*M_PI
              endif
              lon_index_ext = int(lon_geo_scalar_used/delta_lon_ext)
              
              ! computing the subset of the external data to use for the interpolation
              call calc_ext_subset(nlat_ext,nlon_ext,grid%lat_geo_scalar(ji,jk),lat_index_ext,lon_index_ext,n_points_ext_domain, &
                                   lower_index_ext,upper_index_ext,left_index_ext,right_index_ext)
              
              ! looping over all points of the input dataset in the vicinity of the grid cell at hand
              do jm=upper_index_ext,lower_index_ext
                do jn=left_index_ext,right_index_ext
                  
                  ! correcting the indices for boundary cases
                  call correct_ext_data_indices(jm,jn,nlat_ext,nlon_ext,jm_used,jn_used)
                  
                  ! adding the temperature value at hand to the interpolated value if the temperature value is not invalid
                  if (ghcn_cams(jn_used,jm_used,1)/=-9.96921e36) then
                    grid%t_const_soil(ji,jk) = grid%t_const_soil(ji,jk) &
                                               + sum(ghcn_cams(jn_used,jm_used,:))/12._wp + lapse_rate*grid%oro(ji,jk)
                  else
                    invalid_counter(ji,jk) = invalid_counter(ji,jk)+1
                  endif
                  
                enddo
              enddo
              
              ! computing the average
              if (invalid_counter(ji,jk)<n_points_ext_domain) then
                grid%t_const_soil(ji,jk) = grid%t_const_soil(ji,jk)/(n_points_ext_domain-invalid_counter(ji,jk))
              ! this is the case if all input values were invalid
              else
                ! adding realistic values where no real-world values were found
                grid%t_const_soil(ji,jk) = t_0 + 25._wp*cos(2._wp*grid%lat_geo_scalar(ji,jk))
              endif
              
            enddo
          enddo
          !$omp end parallel do
          
          deallocate(ghcn_cams)
          deallocate(invalid_counter)
          
          write(*,*) "Lower boundary soil temperature set."
          
        endif
      
      ! Schaer wave test orography
      case(2)
        height_mountain = 250._wp
        !$omp parallel do private(jk,x_coord)
        do jk=1,nx
          x_coord = dx*jk - (nx/2 + 1)*dx
          grid%oro(:,jk) = height_mountain*exp(-x_coord**2/5000._wp**2)*cos(M_PI*x_coord/4000._wp)**2
          grid%oro_smoothed(:,jk) = 0.5_wp*height_mountain*exp(-x_coord**2/5000._wp**2)
        enddo
        !$omp end parallel do
      
      ! Schaer advection test orography
      case(3)
        height_mountain = 3000._wp
        !$omp parallel do private(jk,x_coord)
        do jk=1,nx
          x_coord = dx*jk - (nx/2 + 1)*dx
          if (abs(x_coord)<=25e3_wp) then
            grid%oro(:,jk) = height_mountain*cos(0.5_wp*M_PI*x_coord/25000._wp)**2*cos(M_PI*x_coord/8000._wp)**2
          else
            grid%oro(:,jk) = 0._wp
          endif
        enddo
        !$omp end parallel do
    
    endselect
    
    !$omp parallel workshare
    grid%z_w(:,:,n_levels) = grid%oro
    !$omp end parallel workshare
    
    ! Other physical properties of the surface
    ! ----------------------------------------
    
    if (.not. lread_geo) then
      
      density_soil = 1442._wp
      c_p_soil = 830._wp
      c_p_water = 4184._wp
      albedo_water = 0.06_wp
      ! setting the land surface albedo to 0.12 (compare Zdunkowski, Trautmann & Bott:
      ! Radiation in the Atmosphere, 2007, p. 444)
      albedo_soil = 0.12_wp
      albedo_ice = 0.8_wp
      t_conductivity_water = 1.4e-7_wp
      t_conductivity_soil = 7.5e-7_wp
      
      !$omp parallel do private(ji,jk)
      do jk=1,nx
        do ji=1,ny
          
          ! seabreeze land-sea mask
          if (trim(scenario)=="seabreeze") then
            if (jk>=nx/4 .and. jk<=3*nx/4) then
              grid%land_fraction(ji,jk) = 1._wp
            endif
          endif
          
          ! albedo of water
          grid%sfc_albedo(ji,jk) = albedo_water
  
          ! for water the roughness length is set to some sea-typical value, will not be used anyway
          grid%roughness_length(ji,jk) = 0.08_wp
          
          ! will also not be used for water
          grid%sfc_rho_c(ji,jk) = rho_h2o*c_p_water
          
          grid%t_conduc_soil(ji,jk) = t_conductivity_water
          
          ! land is present in this grid cell
          if (grid%land_fraction(ji,jk)>EPSILON_SECURITY) then
            
            ! lakes are included in the soil calculation
            grid%t_conduc_soil(ji,jk) = (grid%land_fraction(ji,jk)*t_conductivity_soil &
                                        + grid%lake_fraction(ji,jk)*t_conductivity_water) &
                                        /(grid%land_fraction(ji,jk)+grid%lake_fraction(ji,jk))
            
            ! lakes are included in the soil calculation
            grid%sfc_rho_c(ji,jk) = (grid%land_fraction(ji,jk)*density_soil*c_p_soil &
                                    + grid%lake_fraction(ji,jk)*rho_h2o*c_p_water) &
                                    /(grid%land_fraction(ji,jk)+grid%lake_fraction(ji,jk))
            
            ! setting the surface albedo of land depending on the latitude
            ! ice
            if (abs(360._wp/(2._wp*M_PI)*grid%lat_geo_scalar(ji,jk))>70._wp) then
              grid%sfc_albedo(ji,jk) = grid%land_fraction(ji,jk)*albedo_ice + (1._wp-grid%land_fraction(ji,jk))*albedo_water
            ! normal soil
            else
              grid%sfc_albedo(ji,jk) = grid%land_fraction(ji,jk)*albedo_soil + (1._wp-grid%land_fraction(ji,jk))*albedo_water
            endif
            
            ! calculating a roughness length depending on the vegetation height
            grid%roughness_length(ji,jk) = vegetation_height_ideal(grid%lat_geo_scalar(ji,jk),grid%z_w(ji,jk,n_levels))/8._wp
            
          endif
          
          ! restricting the roughness length to a minimum
          grid%roughness_length(ji,jk) = max(0.0001_wp,grid%roughness_length(ji,jk))
        
        enddo
      enddo
      !$omp end parallel do
      
    endif
    
    ! giving the user some information on the physical surface quantities
    
    !$omp parallel workshare
    dq_value = minval(grid%land_fraction)
    !$omp end parallel workshare
    write(*,*) "minimum land fraction:",dq_value
    
    !$omp parallel workshare
    dq_value = maxval(grid%land_fraction)
    !$omp end parallel workshare
    write(*,*) "maximum land fraction:",dq_value
    
    !$omp parallel workshare
    dq_value = sum(grid%land_fraction)/(ny*nx)
    !$omp end parallel workshare
    write(*,*) "average land fraction:",dq_value
        
    !$omp parallel workshare
    dq_value = minval(grid%lake_fraction)
    !$omp end parallel workshare
    write(*,*) "minimum lake fraction:",dq_value
    
    !$omp parallel workshare
    dq_value = maxval(grid%lake_fraction)
    !$omp end parallel workshare
    write(*,*) "maximum lake fraction:",dq_value
    
    !$omp parallel workshare
    dq_value = sum(grid%lake_fraction)/(ny*nx)
    !$omp end parallel workshare
    write(*,*) "average lake fraction:",dq_value
    
    !$omp parallel workshare
    dq_value = minval(grid%oro)
    !$omp end parallel workshare
    write(*,*) "minimum orography:",dq_value,"m"
    
    !$omp parallel workshare
    dq_value = maxval(grid%oro)
    !$omp end parallel workshare
    write(*,*) "maximum orography:",dq_value,"m"
    
    !$omp parallel workshare
    dq_value = minval(grid%t_const_soil)
    !$omp end parallel workshare
    write(*,*) "minimum background soil temperature:",dq_value,"K"
    
    !$omp parallel workshare
    dq_value = maxval(grid%t_const_soil)
    !$omp end parallel workshare
    write(*,*) "maximum background soil temperature:",dq_value,"K"
    
    !$omp parallel workshare
    dq_value = minval(grid%sfc_rho_c)
    !$omp end parallel workshare
    write(*,*) "minimum volumetric heat capacity of the soil:",dq_value,"J/(m**3K)"
    
    !$omp parallel workshare
    dq_value = maxval(grid%sfc_rho_c)
    !$omp end parallel workshare
    write(*,*) "maximum volumetric heat capacity of the soil:",dq_value,"J/(m**3K)"
    
    !$omp parallel workshare
    dq_value = minval(grid%t_conduc_soil)
    !$omp end parallel workshare
    write(*,*) "minimum temperature conductivity of the soil:",dq_value,"m**2/s"
    
    !$omp parallel workshare
    dq_value = maxval(grid%t_conduc_soil)
    !$omp end parallel workshare
    write(*,*) "maximum temperature conductivity of the soil:",dq_value,"m**2/s"
    
    ! Vertical grid
    ! -------------
    
    ! calculating the vertical positions of the scalar points
    ! the heights are defined according to k = A_k + B_k*surface with A_0 = toa, A_{NO_OF_LEVELS} = 0, B_0 = 0, B_{NO_OF_LEVELS} = 1
    !$omp parallel do private(ji,jk,jl,z_rel,sigma_z,A,B,vertical_vector_pre,max_oro,toa_oro)
    do jk=1,nx
      do ji=1,ny
        ! filling up vertical_vector_pre
        do jl=1,n_levels
          z_rel = 1._wp-(jl-1._wp)/n_layers
          sigma_z = z_rel**stretching_parameter
          A = sigma_z*toa ! the height without orography
          ! including orography
          if (jl>n_flat_layers+1) then
            if (lsleve) then
              toa_oro = vertical_vector_pre(n_flat_layers+1)
              if (orography_id==2) then
                vertical_vector_pre(jl) = A + sinh((toa_oro-A)/5.e3_wp)/sinh(toa_oro/5.e3_wp)*grid%oro_smoothed(ji,jk) &
                                          + sinh((toa_oro-A)/2.e3_wp)/sinh(toa_oro/2.e3_wp) &
                                          *(grid%z_w(ji,jk,n_levels) - grid%oro_smoothed(ji,jk))
              else
                vertical_vector_pre(jl) = A + sinh((toa_oro-A)/(0.5_wp*toa_oro))/sinh(toa_oro/(0.5_wp*toa_oro)) &
                                          *grid%oro_smoothed(ji,jk) &
                                          + sinh((toa_oro-A)/(0.25_wp*toa_oro))/sinh(toa_oro/(0.25_wp*toa_oro)) &
                                          *(grid%z_w(ji,jk,n_levels) - grid%oro_smoothed(ji,jk))
              endif
            else
              B = (jl-(n_flat_layers+1._wp))/n_oro_layers
              vertical_vector_pre(jl) = A + B*grid%z_w(ji,jk,n_levels)
            endif
          else
            vertical_vector_pre(jl) = A
          endif
        
          ! check
          if (jl>1) then
            if (vertical_vector_pre(jl)>=vertical_vector_pre(jl-1)) then
              write(*,*) "Problem in vertical grid generation. You might have to change SLEVE parameters."
              write(*,*) "Aborting."
              call exit(1)
            endif
          endif
        
        enddo
          
        ! check
        if (ji==1 .and. jk==1) then
          max_oro = maxval(grid%z_w(:,:,n_levels))
          if (max_oro>=vertical_vector_pre(n_flat_layers+1)) then
            write(*,*) "Maximum of orography larger or equal to the height of the lowest flat level."
            write(*,*) "Aborting."
            call exit(1)
          endif
        endif
        
        ! placing the scalar points in the middle between the preliminary values of the adjacent levels
        do jl=1,n_layers
          grid%z_scalar(ji,jk,jl) = 0.5_wp*(vertical_vector_pre(jl) + vertical_vector_pre(jl+1))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! setting the height of the u-vector points
    ! inner domain
    !$omp parallel do private(jk)
    do jk=2,nx
      grid%z_u(:,jk,:) = 0.5_wp*(grid%z_scalar(:,jk-1,:) + grid%z_scalar(:,jk,:))
    enddo
    !$omp end parallel do
    ! boundaries
    if (lperiodic) then
      !$omp parallel workshare
      grid%z_u(:,1,:) = 0.5_wp*(grid%z_scalar(:,1,:) + grid%z_scalar(:,nx,:))
      grid%z_u(:,nx+1,:) = grid%z_u(:,1,:)
      !$omp end parallel workshare
    else
      !$omp parallel workshare
      grid%z_u(:,1,:) = grid%z_scalar(:,1,:) + 0.5_wp*(grid%z_scalar(:,1,:) - grid%z_scalar(:,2,:))
      grid%z_u(:,nx+1,:) = grid%z_scalar(:,nx,:) + 0.5_wp*(grid%z_scalar(:,nx,:) - grid%z_scalar(:,nx-1,:))
      !$omp end parallel workshare
    endif
    
    ! setting the height of the v-vector points
    ! inner domain
    !$omp parallel do private(ji)
    do ji=2,ny
      grid%z_v(ji,:,:) = 0.5_wp*(grid%z_scalar(ji-1,:,:) + grid%z_scalar(ji,:,:))
    enddo
    !$omp end parallel do
    ! boundaries
    if (lperiodic) then
      !$omp parallel workshare
      grid%z_v(1,:,:) = 0.5_wp*(grid%z_scalar(1,:,:) + grid%z_scalar(ny,:,:))
      grid%z_v(ny+1,:,:) = grid%z_v(1,:,:)
      !$omp end parallel workshare
    else
      !$omp parallel workshare
      grid%z_v(1,:,:) = grid%z_scalar(1,:,:) + 0.5_wp*(grid%z_scalar(1,:,:) - grid%z_scalar(2,:,:))
      grid%z_v(ny+1,:,:) = grid%z_scalar(ny,:,:) + 0.5_wp*(grid%z_scalar(ny,:,:) - grid%z_scalar(ny-1,:,:))
      !$omp end parallel workshare
    endif
    
    ! setting dx
    !$omp parallel do private(ji)
    do ji=1,ny
      grid%dx(ji,:,:) = dx*cos(grid%lat_scalar(ji))*(r_e + grid%z_u(ji,:,:))/r_e
    enddo
    !$omp end parallel do
    ! setting dy
    !$omp parallel workshare
    grid%dy = dy*(r_e + grid%z_v)/r_e
    !$omp end parallel workshare
    
    ! plane geometry grid distances
    if (lplane) then
      !$omp parallel workshare
      grid%dx = dx
      grid%dy = dy
      !$omp end parallel workshare
    endif
    
    ! setting dx of the dual grid
    !$omp parallel do private(ji)
    do ji=1,ny+1
      if (ji==ny+1) then
        grid%dx_dual(ji,:,:) = dx*cos(grid%lat_scalar(ji-1) - 0.5_wp*dlat)*(r_e + grid%z_v(ji,:,:))/r_e
      else
        grid%dx_dual(ji,:,:) = dx*cos(grid%lat_scalar(ji) + 0.5_wp*dlat)*(r_e + grid%z_v(ji,:,:))/r_e
      endif
    enddo
    !$omp end parallel do
    ! setting dy of the dual grid
    !$omp parallel workshare
    grid%dy_dual = dy*(r_e + grid%z_u)/r_e
    !$omp end parallel workshare
    
    ! plane geometry dual grid distances
    if (lplane) then
      !$omp parallel workshare
      grid%dx_dual = dx
      grid%dy_dual = dy
      !$omp end parallel workshare
    endif
    
    ! calculating the coordinate slopes
    call grad_hor_cov(grid%z_scalar,grid%slope_x,grid%slope_y,grid)
    
    ! setting the z coordinates of the vertical vector points
    !$omp parallel do private(jl)
    do jl=1,n_layers
      if (jl==1) then
        grid%z_w(:,:,jl) = toa
      else
        grid%z_w(:,:,jl) = 0.5_wp*(grid%z_scalar(:,:,jl-1) + grid%z_scalar(:,:,jl))
      endif
    enddo
    !$omp end parallel do
    
    ! calculating n_damping_layers
    n_damping_levels = 0
    do jl=1,n_levels
      !$omp parallel workshare
      max_z = maxval(grid%z_w(:,:,jl))
      !$omp end parallel workshare
      if (max_z>klemp_begin_rel*toa) then
        n_damping_levels = n_damping_levels + 1
      endif
    enddo
    
    ! setting the layer thicknesses
    !$omp parallel do private(jl)
    do jl=1,n_layers
      grid%layer_thickness(:,:,jl) = grid%z_w(:,:,jl) - grid%z_w(:,:,jl+1)
    enddo
    !$omp end parallel do
    
    ! setting the vertical distances between the scalar data points
    !$omp parallel do private(jl)
    do jl=1,n_levels
      if (jl==1) then
        grid%dz(:,:,jl) = 2._wp*(toa - grid%z_scalar(:,:,jl))
      elseif (jl==n_levels) then
        grid%dz(:,:,jl) = 2._wp*(grid%z_scalar(:,:,jl-1) - grid%z_w(:,:,jl))
      else
        grid%dz(:,:,jl) = grid%z_scalar(:,:,jl-1) - grid%z_scalar(:,:,jl)
      endif
    enddo
    !$omp end parallel do
    
    ! setting the horizontal areas at the surface
    !$omp parallel do private(ji,jk)
    do jk=1,nx
      do ji=1,ny
        grid%area_z(ji,jk,n_levels) = patch_area(grid%lat_scalar(ji),dlon,dlat)*(r_e + grid%z_w(ji,jk,n_levels))**2/r_e**2
      enddo
    enddo
    !$omp end parallel do

    ! setting the horizontal areas at the higher points (above the surface)
    !$omp parallel do private(jl)
    do jl=1,n_layers
      grid%area_z(:,:,jl) = grid%area_z(:,:,n_levels)*(r_e + grid%z_w(:,:,jl))**2 &
      /(r_e + grid%z_w(:,:,n_levels))**2
    enddo
    !$omp end parallel do
    
    ! plane geometry horizontal areas
    if (lplane) then
      !$omp parallel workshare
      grid%area_z(:,:,:) = dx*dy
      !$omp end parallel workshare
    endif
    
    ! some DQ
    !$omp parallel workshare
    dq_value = sum(grid%land_fraction*grid%area_z(:,:,n_levels))/sum(grid%area_z(:,:,n_levels))
    !$omp end parallel workshare
    write(*,*) "share of the surface covered by land:",dq_value
    
    ! the mean velocity area can be set now
    grid%mean_velocity_area = 2._wp*sum(grid%area_z(:,:,n_levels))/size(grid%area_z(:,:,n_levels))
    
    ! setting the vertical areas in x-direction
    !$omp parallel do private(ji,jk,jl,lower_z,upper_z,lower_length)
    do jl=1,n_layers
      do jk=1,nx+1
        do ji=1,ny
          ! left boundary
          if (jk==1) then
            lower_z = grid%z_w(ji,1,jl+1)+0.5_wp*(grid%z_w(ji,1,jl+1)-grid%z_w(ji,2,jl+1))
            upper_z = grid%z_w(ji,1,jl)+0.5_wp*(grid%z_w(ji,1,jl)-grid%z_w(ji,2,jl))
          ! right boundary
          elseif (jk==nx+1) then
            lower_z = grid%z_w(ji,nx,jl+1)+0.5_wp*(grid%z_w(ji,nx,jl+1)-grid%z_w(ji,nx-1,jl+1))
            upper_z = grid%z_w(ji,nx,jl)+0.5_wp*(grid%z_w(ji,nx,jl)-grid%z_w(ji,nx-1,jl))
          ! inner domain
          else
            lower_z = 0.5_wp*(grid%z_w(ji,jk-1,jl+1) + grid%z_w(ji,jk,jl+1))
            upper_z = 0.5_wp*(grid%z_w(ji,jk-1,jl) + grid%z_w(ji,jk,jl))
          endif
          lower_length = dy*(r_e+lower_z)/r_e
          ! plane geometry
          if (lplane) then
            lower_length = dy
          endif
          grid%area_x(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! setting the vertical areas in y-direction
    !$omp parallel do private(ji,jk,jl,lower_z,upper_z,lower_length)
    do jl=1,n_layers
      do jk=1,nx
        do ji=1,ny+1
          ! upper boundary
          if (ji==1) then
            lower_z = grid%z_w(1,jk,jl+1)+0.5_wp*(grid%z_w(1,jk,jl+1)-grid%z_w(2,jk,jl+1))
            upper_z = grid%z_w(1,jk,jl)+0.5_wp*(grid%z_w(1,jk,jl)-grid%z_w(2,jk,jl))
          ! lower boundary
          elseif (ji==ny+1) then
            lower_z = grid%z_w(ny,jk,jl+1)+0.5_wp*(grid%z_w(ny,jk,jl+1)-grid%z_w(ny-1,jk,jl+1))
            upper_z = grid%z_w(ny,jk,jl)+0.5_wp*(grid%z_w(ny,jk,jl)-grid%z_w(ny-1,jk,jl))
          ! inner domain
          else
            lower_z = 0.5_wp*(grid%z_w(ji-1,jk,jl+1) + grid%z_w(ji,jk,jl+1))
            upper_z = 0.5_wp*(grid%z_w(ji-1,jk,jl) + grid%z_w(ji,jk,jl))
          endif
          if (ji==ny+1) then
            lower_length = dx*cos(grid%lat_scalar(ji-1)-0.5_wp*dlat)*(r_e+lower_z)/r_e
          else
            lower_length = dx*cos(grid%lat_scalar(ji)+0.5_wp*dlat)*(r_e+lower_z)/r_e
          endif
          ! plane geometry
          if (lplane) then
            lower_length = dx
          endif
          grid%area_y(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! setting the horizontal dual areas
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do jk=1,nx+1
        do ji=1,ny+1
        
          ! setting the vertical position of the areas
          if (jk==1) then
            grid%z_area_dual_z(ji,jk,jl) = grid%z_v(ji,1,jl) + 0.5_wp*(grid%z_v(ji,1,jl)-grid%z_v(ji,2,jl))
          elseif (jk==nx+1) then
            grid%z_area_dual_z(ji,jk,jl) = grid%z_v(ji,nx,jl) + 0.5_wp*(grid%z_v(ji,nx,jl)-grid%z_v(ji,nx-1,jl))
          else
            grid%z_area_dual_z(ji,jk,jl) = 0.5_wp*(grid%z_v(ji,jk-1,jl)+grid%z_v(ji,jk,jl))
          endif
          
          ! setting the area itself
          if (ji==ny+1) then
            grid%area_dual_z(ji,jk,jl) = patch_area(grid%lat_scalar(ji-1) - 0.5_wp*dlat,dlon,dlat) &
            *(r_e + grid%z_area_dual_z(ji,jk,jl))**2/r_e**2
          else
            grid%area_dual_z(ji,jk,jl) = patch_area(grid%lat_scalar(ji) + 0.5_wp*dlat,dlon,dlat) &
            *(r_e + grid%z_area_dual_z(ji,jk,jl))**2/r_e**2
          endif
          ! plane geometry
          if (lplane) then
            grid%area_dual_z(ji,jk,jl) = dx*dy
          endif
          
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! setting the vertical dual areas in x-direction
    !$omp parallel do private(ji,jk,jl,lower_z,lower_length,upper_z)
    do jl=1,n_levels
      do jk=1,nx
        do ji=1,ny+1
          if (jl==n_levels) then
            if (ji==1) then
              lower_z = grid%z_w(1,jk,jl) + 0.5_wp*(grid%z_w(1,jk,jl)-grid%z_w(2,jk,jl))
            elseif (ji==ny+1) then
              lower_z = grid%z_w(ny,jk,jl) + 0.5_wp*(grid%z_w(ny,jk,jl)-grid%z_w(ny-1,jk,jl))
            else
              lower_z = 0.5_wp*(grid%z_w(ji-1,jk,jl) + grid%z_w(ji,jk,jl))
            endif
            lower_length = grid%dy(ji,jk,n_layers)*(r_e + lower_z)/(r_e + grid%z_v(ji,jk,n_layers))
          else
            lower_z = grid%z_v(ji,jk,jl)
            lower_length = grid%dy(ji,jk,jl)
          endif
          if (jl==1) then
            if (ji==1) then
              upper_z = grid%z_w(1,jk,jl) + 0.5_wp*(grid%z_w(1,jk,jl)-grid%z_w(2,jk,jl))
            elseif (ji==ny+1) then
              upper_z = grid%z_w(ny,jk,jl) + 0.5_wp*(grid%z_w(ny,jk,jl)-grid%z_w(ny-1,jk,jl))
            else
              upper_z = 0.5_wp*(grid%z_w(ji-1,jk,jl) + grid%z_w(ji,jk,jl))
            endif
          else
            upper_z = grid%z_v(ji,jk,jl-1)
          endif
          ! plane geometry
          if (lplane) then
            lower_length=dy
          endif
          grid%area_dual_x(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! setting the vertical dual areas in y-direction
    !$omp parallel do private(ji,jk,jl,lower_z,lower_length,upper_z)
    do jl=1,n_levels
      do jk=1,nx+1
        do ji=1,ny
          if (jl==n_levels) then
            if (jk==1) then
              lower_z = grid%z_w(ji,1,jl) + 0.5_wp*(grid%z_w(ji,1,jl)-grid%z_w(ji,2,jl))
            elseif (jk==nx+1) then
              lower_z = grid%z_w(ji,nx,jl) + 0.5_wp*(grid%z_w(ji,nx,jl)-grid%z_w(ji,nx-1,jl))
            else
              lower_z = 0.5_wp*(grid%z_w(ji,jk-1,jl) + grid%z_w(ji,jk,jl))
            endif
            lower_length = grid%dx(ji,jk,n_layers)*(r_e + lower_z)/(r_e + grid%z_u(ji,jk,n_layers))
          else
            lower_z = grid%z_u(ji,jk,jl)
            lower_length = grid%dx(ji,jk,jl)
          endif
          if (jl==1) then
            if (jk==1) then
              upper_z = grid%z_w(ji,1,jl) + 0.5_wp*(grid%z_w(ji,1,jl)-grid%z_w(ji,2,jl))
            elseif (jk==nx+1) then
              upper_z = grid%z_w(ji,nx,jl) + 0.5_wp*(grid%z_w(ji,nx,jl)-grid%z_w(ji,nx-1,jl))
            else
              upper_z = 0.5_wp*(grid%z_w(ji,jk-1,jl) + grid%z_w(ji,jk,jl))
            endif
          else
            upper_z = grid%z_u(ji,jk,jl-1)
          endif
          ! plane geometry
          if (lplane) then
            lower_length = dx
          endif
          grid%area_dual_y(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! setting the volume of the grid boxes
    !$omp parallel do private(jl)
    do jl=1,n_layers
      grid%volume(:,:,jl) = 1._wp/3._wp*((r_e + grid%z_w(:,:,jl))**3 - (r_e + grid%z_w(:,:,jl+1))**3) &
      /(r_e + grid%z_w(:,:,jl+1))**2*grid%area_z(:,:,jl+1)
      ! plane geometry
      if (lplane) then
        grid%volume(:,:,jl) = grid%area_z(:,:,jl+1)*(grid%z_w(:,:,jl) - grid%z_w(:,:,jl+1))
      endif
    enddo
    !$omp end parallel do
    
    ! Derived quantities
    ! ------------------
    
    ! setting the inner product weights
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do jk=1,nx
        do ji=1,ny
          grid%inner_product_weights(1,ji,jk,jl) = grid%area_x(ji,jk+1,jl)*grid%dx(ji,jk+1,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(2,ji,jk,jl) = grid%area_y(ji,jk,jl)*grid%dy(ji,jk,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(3,ji,jk,jl) = grid%area_x(ji,jk,jl)*grid%dx(ji,jk,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(4,ji,jk,jl) = grid%area_y(ji+1,jk,jl)*grid%dy(ji+1,jk,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(5,ji,jk,jl) = grid%area_z(ji,jk,jl)*grid%dz(ji,jk,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(6,ji,jk,jl) = grid%area_z(ji,jk,jl+1)*grid%dz(ji,jk,jl+1)/(2._wp*grid%volume(ji,jk,jl))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! setting the TRSK weights
    ! u
    do ji=1,ny
      base_area = patch_area(grid%lat_scalar(ji),dlon,dlat)
      grid%trsk_weights_u(1,ji) = (0.5_wp - patch_area(grid%lat_scalar(ji)+0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area)
      if (lplane) then
        grid%trsk_weights_u(1,ji) = grid%trsk_weights_u(1,ji)*dx/dx
      else
        grid%trsk_weights_u(1,ji) = grid%trsk_weights_u(1,ji)*dx*cos(grid%lat_scalar(ji)+0.5_wp*dlat)/(dx*cos(grid%lat_scalar(ji)))
      endif
      grid%trsk_weights_u(2,ji) = -(0.5_wp - patch_area(grid%lat_scalar(ji)+0.25_wp*dlat,dlon,0.5_wp*dlat)/base_area)
      if (lplane) then
        grid%trsk_weights_u(2,ji) = grid%trsk_weights_u(2,ji)*dy/dx
      else
        grid%trsk_weights_u(2,ji) = grid%trsk_weights_u(2,ji)*dy/(dx*cos(grid%lat_scalar(ji)))
      endif
      grid%trsk_weights_u(3,ji) = -(0.5_wp - (patch_area(grid%lat_scalar(ji)+0.25_wp*dlat,dlon,0.5_wp*dlat) &
      + patch_area(grid%lat_scalar(ji)-0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat))/base_area)
      if (lplane) then
        grid%trsk_weights_u(3,ji) = grid%trsk_weights_u(3,ji)*dx/dx
      else
        grid%trsk_weights_u(3,ji) = grid%trsk_weights_u(3,ji)*dx*cos(grid%lat_scalar(ji)-0.5_wp*dlat)/(dx*cos(grid%lat_scalar(ji)))
      endif
    enddo
    ! v
    do ji=1,ny+1
      if (ji==ny+1) then
        lat_lower_center = grid%lat_scalar(ny)-dlat
      else
        lat_lower_center = grid%lat_scalar(ji)
      endif
      base_area = patch_area(lat_lower_center,dlon,dlat)
      grid%trsk_weights_v(1,ji) = -(0.5_wp - patch_area(lat_lower_center+0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area)
      if (ji==1) then
        lat_upper_center = grid%lat_scalar(1)+dlat
      else
        lat_upper_center = grid%lat_scalar(ji-1)
      endif
      base_area = patch_area(lat_upper_center,dlon,dlat)
      grid%trsk_weights_v(2,ji) = -(0.5_wp - patch_area(lat_upper_center-0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area)
    enddo
    
    ! Soil grid
    ! ---------
    
    sigma_soil = 0.352_wp
    grid%z_t_const = -10._wp

    ! the surface is always at zero
    grid%z_soil_interface(1) = 0._wp
    do jl=2,nsoillays+1
      grid%z_soil_interface(jl) = grid%z_soil_interface(jl-1) + sigma_soil**(nsoillays-jl)
    enddo
    
    rescale_factor = grid%z_t_const/grid%z_soil_interface(nsoillays+1)
    
    do jl=2,nsoillays+1
      grid%z_soil_interface(jl) = rescale_factor*grid%z_soil_interface(jl)
    enddo
    do jl=1,nsoillays
      grid%z_soil_center(jl) = 0.5_wp*(grid%z_soil_interface(jl) + grid%z_soil_interface(jl+1))
    enddo
    
    write(*,*) "Thickness of the uppermost soil layer: ", -grid%z_soil_interface(2), "m."

    ! writing some grid properties to a file if required by the user
    if (lwrite_grid) then
      call write_grid(grid)
    endif
    
  end subroutine grid_setup
  
  subroutine bg_setup(grid)
    
    ! This subroutine sets up the background state.
    
    type(t_grid), intent(inout) :: grid ! the model grid
    
    ! local variables
    real(wp) :: pressure ! pressure at the respective gridpoint
    real(wp) :: b        ! parameter needed for the hydrostatic initialization routine
    real(wp) :: c        ! parameter needed for the hydrostatic initialization routine
    integer  :: ji       ! horizontal index
    integer  :: jk       ! horizontal index
    integer  :: jl       ! layer index
    
    ! integrating the hydrostatic background state according to the given temperature profile and pressure in the lowest layer
    !$omp parallel do private(ji,jk,jl,b,c,pressure)
    do jk=1,nx
      do ji=1,ny
        ! integrating from bottom to top
        do jl=n_layers,1,-1
          ! setting the geopotential
          grid%gravity_potential(ji,jk,jl) = geopot(grid%z_scalar(ji,jk,jl))
          ! lowest layer
          if (jl==n_layers) then
            pressure = bg_pres(grid%z_scalar(ji,jk,jl))
            grid%exner_bg(ji,jk,jl) = (pressure/p_0)**(r_d/c_d_p)
          ! other layers
          else
            ! solving a quadratic equation for the Exner pressure
            b = -0.5_wp*grid%exner_bg(ji,jk,jl+1)/bg_temp(grid%z_scalar(ji,jk,jl+1)) &
            *(bg_temp(grid%z_scalar(ji,jk,jl)) - bg_temp(grid%z_scalar(ji,jk,jl+1)) + 2.0_wp/ &
            c_d_p*(grid%gravity_potential(ji,jk,jl) - grid%gravity_potential(ji,jk,jl+1)))
            c = grid%exner_bg(ji,jk,jl+1)**2*bg_temp(grid%z_scalar(ji,jk,jl))/bg_temp(grid%z_scalar(ji,jk,jl+1))
            grid%exner_bg(ji,jk,jl) = b+sqrt(b**2+c)
          endif
          grid%theta_v_bg(ji,jk,jl) = bg_temp(grid%z_scalar(ji,jk,jl))/grid%exner_bg(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! calculating the gradient of the background Exner pressure (only needs to be done once)
    call grad_vert(grid%exner_bg,grid%exner_bg_grad_w,grid)
    call grad_hor(grid%exner_bg,grid%exner_bg_grad_u,grid%exner_bg_grad_v,grid%exner_bg_grad_w,grid)
    
    ! calculating the vertical acceleration due to gravity
    call grad_vert(grid%gravity_potential,grid%gravity_m_v,grid)
    
  end subroutine bg_setup
  
  subroutine smooth_hor_scalar(array)
    
    ! This subroutine smoothes a scalar field on one layer.
    
    real(wp), intent(inout) :: array(ny,nx) ! The field to smooth.
    
    ! local variables
    real(wp) :: original_array(ny,nx) ! the unsmoothed input array
    integer  :: ji                    ! horizontal index
    integer  :: jk                    ! horizontal index
    
    ! copying the original array
    original_array = array
    
    ! inner domin
    !$omp parallel do private(ji,jk)
    do jk=2,nx-1
      do ji=2,ny-1
        array(ji,jk) = sum(original_array(ji-1:ji+1,jk-1:jk+1))/9._wp
      enddo
    enddo
    !$omp end parallel do
    
    ! corners
    array(1,1) = sum(original_array(1:2,1:2))/4._wp
    array(1,nx) = sum(original_array(1:2,nx-1:nx))/4._wp
    array(ny,1) = sum(original_array(ny-1:ny,1:2))/4._wp
    array(ny,nx) = sum(original_array(ny-1:ny,nx-1:nx))/4._wp
    
    ! boundaries
    !$omp parallel do private(jk)
    do jk=2,nx-1
      array(1,jk) = sum(original_array(1:2,jk-1:jk+1))/6._wp
    enddo
    !$omp end parallel do
    !$omp parallel do private(jk)
    do jk=2,nx-1
      array(ny,jk) = sum(original_array(ny-1:ny,jk-1:jk+1))/6._wp
    enddo
    !$omp end parallel do
    !$omp parallel do private(ji)
    do ji=2,ny-1
      array(ji,1) = sum(original_array(ji-1:ji+1,1:2))/6._wp
    enddo
    !$omp end parallel do
    !$omp parallel do private(ji)
    do ji=2,ny-1
      array(ji,nx) = sum(original_array(ji-1:ji+1,nx-1:nx))/6._wp
    enddo
    !$omp end parallel do
  
  end subroutine smooth_hor_scalar
  
  function find_min_index(input_vector)
    
    ! This function finds the index where a vector assumes its minimum.
    
    real(wp) :: input_vector(:) ! the vector to operate with
    integer  :: find_min_index  ! the result
    
    ! local variables
    integer  :: ji
    real(wp) :: current_min
    
    find_min_index = 1
    current_min = input_vector(1)
    
    do ji=2,size(input_vector)
      if (input_vector(ji)<current_min) then
        current_min = input_vector(ji)
        find_min_index = ji
      endif
    enddo
    
  end function find_min_index
  
  function patch_area(center_lat,dx_as_angle,dy_as_angle)
  
    ! This function calculates the area of a quadrilateral grid cell.
  
    real(wp) :: center_lat  ! latitude at the center of the patch
    real(wp) :: dx_as_angle ! delta x as angle
    real(wp) :: dy_as_angle ! delta y as angle
    real(wp) :: patch_area  ! the result
  
    ! computing the result
    patch_area = r_e**2*dx_as_angle*(sin(center_lat + 0.5_wp*dy_as_angle) - sin(center_lat - 0.5_wp*dy_as_angle))
    
    ! plane geometry
    if (lplane) then
      patch_area = r_e**2*dx_as_angle*dy_as_angle
    endif
  
  end function patch_area
  
  function vertical_face_area(lower_z,upper_z,lower_length)
    
    ! This function calculates the area of a vertical face.
    
    real(wp) :: lower_z            ! geometric height of the lower boundary of the face
    real(wp) :: upper_z            ! geometric height of the upper boundary of the face
    real(wp) :: lower_length       ! length of the lower boundary of the face
    real(wp) :: vertical_face_area ! the result
    
    vertical_face_area = 0.5_wp*lower_length*(r_e + upper_z + r_e + lower_z) &
    /(r_e + lower_z)*(upper_z - lower_z)
    
    ! plane geometry
    if (lplane) then
      vertical_face_area = lower_length*(upper_z - lower_z)
    endif
    
  end function vertical_face_area
  
  function vegetation_height_ideal(latitude,oro)
    
    ! This function calculates a latitude- and height-dependant idealized vegetation height.

    real(wp) :: latitude                ! latitude of this point
    real(wp) :: oro                     ! height of the terrain at this point
    real(wp) :: vegetation_height_ideal ! result
    
    ! local variables
    real(wp) :: vegetation_height_equator ! the vegetation height at the equator

    ! setting the vegetation height at the equator
    vegetation_height_equator = 20._wp

    vegetation_height_ideal = vegetation_height_equator*cos(latitude)*exp(-oro/1500._wp)

  end function vegetation_height_ideal
  
  subroutine find_global_normal(lat,lon,r)
    
    ! This subroutine calculates the Cartesian normal vector of a point given its geographical coordinates.
    
    real(wp), intent(in)  :: lat  ! input latitude
    real(wp), intent(in)  :: lon  ! input longitude
    real(wp), intent(out) :: r(3) ! positional vector (result)

    r(1) = cos(lat)*cos(lon)
    r(2) = cos(lat)*sin(lon)
    r(3) = sin(lat)

  end subroutine find_global_normal
  
  subroutine find_geos(r,lat_out,lon_out)
    
    ! This subroutine calculates the geographical coordinates of a point given its Cartesian coordinates.
    
    real(wp), intent(in)  :: r(3)    ! positional vector
    real(wp), intent(out) :: lat_out ! output latitude
    real(wp), intent(out) :: lon_out ! output longitude
    
    lat_out = asin(r(3)/sqrt(r(1)**2+r(2)**2+r(3)**2))
    lon_out = atan2(r(2),r(1))

  end subroutine find_geos

  subroutine calc_local_i(lon,result_vec)

    ! This subroutine calculates the local eastward basis vector.
    
    real(wp), intent(in)  :: lon           ! geographical longitude
    real(wp), intent(out) :: result_vec(3) ! result
    
    result_vec(1) = -sin(lon)
    result_vec(2) = cos(lon)
    result_vec(3) = 0._wp
  
  end subroutine calc_local_i

  subroutine calc_local_j(lat,lon,result_vec)

    ! This subroutine calculates the local northward basis vector.
    
    real(wp), intent(in)  :: lat           ! geographical latitude
    real(wp), intent(in)  :: lon           ! geographical longitude
    real(wp), intent(out) :: result_vec(3) ! result
    
    result_vec(1) = -sin(lat)*cos(lon)
    result_vec(2) = -sin(lat)*sin(lon)
    result_vec(3) = cos(lat)

  end subroutine calc_local_j
  
  subroutine calc_ext_subset(nlat_ext,nlon_ext,lat_c_value,lat_index_ext,lon_index_ext,n_points_ext_domain, &
                             lower_index_ext,upper_index_ext,left_index_ext,right_index_ext)
    
    ! This subroutine calculates which subset of the external dataset to use to interpolate to a grid cell.
    
    integer,  intent(in)    :: nlat_ext            ! number of points on the latitude axis of the external data grid
    integer,  intent(in)    :: nlon_ext            ! number of points on the longitude axis of the external data grid
    real(wp), intent(in)    :: lat_c_value         ! latitude of the center of the cell to which to interpolate
    integer,  intent(inout) :: lat_index_ext       ! closest latitude index of the external data
    integer,  intent(inout) :: lon_index_ext       ! closest longitude index of the external data
    integer,  intent(out)   :: lower_index_ext     ! index of the lower latitude boundary of the computed subset of the external data
    integer,  intent(out)   :: upper_index_ext     ! index of the upper latitude boundary of the computed subset of the external data
    integer,  intent(out)   :: left_index_ext      ! longitude index of the western boundary of the computed subset of the external data
    integer,  intent(out)   :: right_index_ext     ! longitude index of the eastern boundary of the computed subset of the external data
    integer,  intent(inout) :: n_points_ext_domain ! number of points of the computed subset of the external data
    
    ! local variables
    integer  :: lat_index_span_ext ! latitude index width of the subset of the external data
    integer  :: lon_index_span_ext ! longitude index width of the subset of the external data
    real(wp) :: delta_lat_ext      ! latitude resolution of the external grid
    real(wp) :: delta_lon_ext      ! longitude resolution of the external grid
    
    ! computing helper variables
    delta_lat_ext = M_PI/nlat_ext
    delta_lon_ext = 2._wp*M_PI/nlon_ext
    lat_index_span_ext = int(eff_hor_res/(r_e*delta_lat_ext))
    
    lat_index_ext = max(1,lat_index_ext)
    lat_index_ext = min(nlat_ext,lat_index_ext)
    lon_index_ext = max(1,lon_index_ext)
    lon_index_ext = min(nlon_ext,lon_index_ext)
    
    lon_index_span_ext = int(min(eff_hor_res/(r_e*delta_lon_ext*max(cos(lat_c_value),EPSILON_SECURITY)),0._wp+nlon_ext))
    lon_index_span_ext = min(lon_index_span_ext,nlon_ext)
    
    n_points_ext_domain = (lat_index_span_ext+1)*(lon_index_span_ext+1)
    
    lower_index_ext = lat_index_ext + lat_index_span_ext/2
    upper_index_ext = lat_index_ext - lat_index_span_ext/2
    left_index_ext = lon_index_ext - lon_index_span_ext/2
    right_index_ext = lon_index_ext + lon_index_span_ext/2
    
    ! updating n_points_ext_domain
    n_points_ext_domain = (lower_index_ext-upper_index_ext+1)*(right_index_ext-left_index_ext+1)
    
  end subroutine calc_ext_subset
  
  subroutine correct_ext_data_indices(jm,jn,nlat_ext,nlon_ext,jm_used,jn_used)
    
    ! This subroutine calculates which indices of an external dataset to actually use.
    
    integer, intent(in)  :: jm       ! latitude index
    integer, intent(in)  :: jn       ! longitude index
    integer, intent(in)  :: nlat_ext ! maximum latitude index
    integer, intent(in)  :: nlon_ext ! maximum longitude index
    integer, intent(out) :: jm_used  ! corrected latitude index
    integer, intent(out) :: jn_used  ! corrected longitude index
    
    jm_used = jm
    if (jm_used<1) then
      jm_used = 1
    endif
    if (jm_used>nlat_ext) then
      jm_used = nlat_ext
    endif
    jn_used = jn
    if (jn_used<1) then
      jn_used = jn_used + nlon_ext
    endif
    if (jn_used>nlon_ext) then
      jn_used = jn_used - nlon_ext
    endif
            
  end subroutine correct_ext_data_indices

end module mo_grid_generator






