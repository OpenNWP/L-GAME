! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

! This file contains the calculation of the grid properties as well as some type definitions.

module grid_generator

  use netcdf
  use definitions,        only: wp,t_grid
  use run_nml,            only: nlins,ncols,nlays,dy,dx,toa,nlays_oro,sigma,scenario,lat_center, &
                                lon_center,lplane
  use constants,          only: re,density_water,T_0,M_PI,p_0,omega,gravity,p_0_standard, &
                                lapse_rate,surface_temp,tropo_height,inv_height,t_grad_inv
  use surface_nml,        only: nsoillays,orography_id
  use gradient_operators, only: grad,grad_hor_cov
  use dictionary,         only: specific_gas_constants,spec_heat_capacities_p_gas
  use io_nml,             only: lwrite_grid,lread_oro,lread_land_sea,lset_oro, &
                                oro_raw_filename
  use read_write_grid,    only: write_grid,read_oro,read_land_sea
  use set_initial_state,  only: bg_temp,bg_pres,geopot
  use bc_nml,             only: lperiodic
  use set_initial_state,  only: nc_check

  implicit none
  
  private
  
  public :: grid_setup
  public :: bg_setup
  
  contains
  
  subroutine grid_setup(grid)
  
    type(t_grid), intent(inout) :: grid ! the model grid
    ! local variables
    real(wp) :: lat_left_upper               ! latitude coordinate of upper left corner
    real(wp) :: lon_left_upper               ! longitude coordinate of upper left corner
    real(wp) :: dlat                         ! mesh size in y direction as angle
    real(wp) :: dlon                         ! mesh size in x direction as angle
    real(wp) :: max_oro                      ! variable for orography check
    real(wp) :: A                            ! variable for calculating the vertical grid
    real(wp) :: B                            ! variable for calculating the vertical grid
    real(wp) :: sigma_z                      ! variable for calculating the vertical grid
    real(wp) :: z_rel                        ! variable for calculating the vertical grid
    real(wp) :: vertical_vector_pre(nlays+1) ! variable for calculating the vertical grid
    real(wp) :: base_area                    ! variable for calculating the vertical grid
    real(wp) :: lower_z,upper_z,lower_length ! variables needed for area calculations
    real(wp) :: height_mountain              ! height of Gaussian mountain (needed for test case)
    real(wp) :: sigma_mountain               ! standard deviation of Gaussian mountain (needed for test case)
    real(wp) :: x_coord                      ! help variable needed for the Schär test
    real(wp) :: rescale_factor               ! soil grid rescaling factor
    real(wp) :: sigma_soil                   ! sigma of the soil grid
    real(wp) :: density_soil                 ! typical density of soil
    real(wp) :: c_p_soil                     ! typical c_p of soil
    real(wp) :: c_p_water                    ! typical c_p of water
    real(wp) :: lat_lower_center             ! variable for calculating the TRSK weights
    real(wp) :: lat_upper_center             ! variable for calculating the TRSK weights
    real(wp) :: rot_y(3,3)                   ! rotation matrix around the global y-axis
    real(wp) :: rot_z(3,3)                   ! rotation matrix around the global z-axis
    real(wp) :: rot(3,3)                     ! complete rotation matrix
    real(wp) :: r_old(3)                     ! positional vector before rotation
    real(wp) :: r_new(3)                     ! positional vector after rotation
    real(wp) :: basis_old(3)                 ! old local basis vector
    real(wp) :: basis_new(3)                 ! new local basis vector
    real(wp) :: local_i(3)                   ! local i-vector
    real(wp) :: local_j(3)                   ! local j-vector
    real(wp) :: x_basis_local,y_basis_local  ! local Cartesian components of the local basis vector
    real(wp) :: lat_local,lon_local          ! helper variables
    integer  :: ji,jk,jl                     ! loop indices
    
    ! Horizontal grid properties
    ! --------------------------
    
    ! setting the center and direction of the grid
    grid%lat_center = lat_center
    grid%lon_center = lon_center
    
    ! setting the latitude and longitude coordinates of the scalar grid points
    ! setting the dy of the model grid
    dlat = dy/re
    dlon = dx/re
    lat_left_upper = (nlins-1._wp)/2._wp*dlat
    lon_left_upper = -(ncols-1._wp)/2._wp*dlon
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji)
    do ji=1,nlins
      grid%lat_scalar(ji) = lat_left_upper - dlat*(ji-1._wp)
      ! this will be modified later
      grid%lat_geo_scalar(ji,:) = grid%lat_scalar(ji)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=1,ncols
      grid%lon_scalar(jk) = lon_left_upper + dlon*(jk-1._wp)
      ! this will be modified later
      grid%lon_geo_scalar(:,jk) = grid%lon_scalar(jk)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! this will be modified later
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins
      do jk=1,ncols+1
        grid%lat_geo_u(ji,jk) = grid%lat_scalar(ji)
        if (jk==ncols+1) then
          grid%lon_geo_u(ji,jk) = grid%lon_scalar(jk-1) + 0.5_wp*dlon
        else
          grid%lon_geo_u(ji,jk) = grid%lon_scalar(jk) - 0.5_wp*dlon
        endif
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! this will be modified later
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins+1
      do jk=1,ncols
        if (ji==nlins+1) then
          grid%lat_geo_v(ji,jk) = grid%lat_scalar(ji-1) - 0.5_wp*dlat
        else
          grid%lat_geo_v(ji,jk) = grid%lat_scalar(ji) + 0.5_wp*dlat
        endif
        grid%lon_geo_v(ji,jk) =  grid%lon_scalar(jk)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! this will be modified later
    !$OMP PARALLEL
    !$OMP WORKSHARE
    grid%dir_geo_u(:,:) = 0._wp
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    
    ! this will be modified later
    !$OMP PARALLEL
    !$OMP WORKSHARE
    grid%dir_geo_v(:,:) = 0.5_wp*M_PI
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    
    ! this will be modified later
    !$OMP PARALLEL
    !$OMP WORKSHARE
    grid%dir_geo_u_scalar(:,:) = 0._wp
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    
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
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,r_old,r_new,basis_old,basis_new,local_i,local_j,x_basis_local,y_basis_local)
    do ji=1,nlins
      do jk=1,ncols
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
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! u-vector points, including directions
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,r_old,r_new,basis_old,basis_new,local_i,local_j,x_basis_local,y_basis_local)
    do ji=1,nlins
      do jk=1,ncols+1
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
    !$OMP END DO
    !$OMP END PARALLEL
	
    ! v-vector points, including directions
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,r_old,r_new,basis_old,basis_new,local_i,local_j,x_basis_local,y_basis_local)
    do ji=1,nlins+1
      do jk=1,ncols
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
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the Coriolis vector
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,y_basis_local)
    do ji=1,nlins+1
      do jk=1,ncols
        y_basis_local = sin(grid%dir_geo_v(ji,jk) - 0.5_wp*M_PI)
        grid%fvec_x(ji,jk) = 2._wp*omega*cos(grid%lat_geo_v(ji,jk))*y_basis_local
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,y_basis_local)
    do ji=1,nlins
      do jk=1,ncols+1
        y_basis_local = sin(grid%dir_geo_u(ji,jk) + 0.5_wp*M_PI)
        grid%fvec_y(ji,jk) = 2._wp*omega*cos(grid%lat_geo_u(ji,jk))*y_basis_local
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,lat_local,lon_local,r_old,r_new)
    do ji=1,nlins+1
      do jk=1,ncols+1
        if (ji==nlins+1) then
          lat_local = grid%lat_scalar(ji-1) - 0.5*dlat
        else
          lat_local = grid%lat_scalar(ji) + 0.5*dlat
        endif
        if (jk==ncols+1) then
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
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! Physical grid properties
    ! ------------------------
    
    ! reading the land-sea mask if configured to do so
    if (lread_land_sea) then
      call read_land_sea(grid)
    endif
    
    select case (orography_id)
    
      ! no orography
      case(0)
        !$OMP PARALLEL
        !$OMP WORKSHARE
        grid%z_w(:,:,nlays+1) = 0._wp
        !$OMP END WORKSHARE
        !$OMP END PARALLEL
      
      ! real orography
      case(1)
        if (lread_oro) then
          call read_oro(grid)
        elseif (lset_oro) then
          call set_orography(grid)
          call smooth_hor_scalar(grid%z_w(:,:,nlays+1))
        endif
      
      ! Schaer orography
      case(2)
        height_mountain = 250._wp
        sigma_mountain = 5000._wp/sqrt(2._wp)
        !$OMP PARALLEL
        !$OMP DO PRIVATE(jk,x_coord)
        do jk=1,ncols
          x_coord = dx*jk - (ncols/2 + 1)*dx
          grid%z_w(:,jk,nlays+1) = height_mountain*exp(-x_coord**2/(2._wp*sigma_mountain**2))*cos(M_PI*x_coord/4000._wp)**2
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    
    endselect
    
    ! idealized soil properties are being set here
    density_soil = 1442._wp
    c_p_soil = 830._wp
    c_p_water = 4184._wp
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins
      do jk=1,ncols
		    
       grid%t_const_soil(ji,jk) = T_0 + 25._wp*cos(2._wp*grid%lat_geo_scalar(ji,jk))
            
        ! albedo of water
        grid%sfc_albedo(ji,jk) = 0.06_wp

        ! for water, the roughness_length is set to some sea-typical value, will not be used anyway
        grid%roughness_length(ji,jk) = 0.08_wp
		
        ! will also not be used for water
        grid%sfc_rho_c(ji,jk) = density_water*c_p_water
             
        ! will only be used for land points
        grid%t_conduc_soil(ji,jk) = 7.5e-7_wp
		
		! land
        if (grid%is_land(ji,jk)==1) then
        
          grid%sfc_rho_c(ji,jk) = density_soil*c_p_soil
          
          ! setting the surface albedo of land depending on the latitude
          ! ice
          if (abs(360._wp/(2._wp*M_PI)*grid%lat_geo_scalar(ji,jk)) > 70._wp) then
            grid%sfc_albedo(ji,jk) = 0.8_wp
          ! normal soil
          else
            grid%sfc_albedo(ji,jk) = 0.12_wp
          endif
          
          ! calculating a roughness length depending on the vegetation height
          grid%roughness_length(ji,jk) = vegetation_height_ideal(grid%lat_geo_scalar(ji,jk),grid%z_w(ji,jk,nlays+1))/8._wp
          
        endif
		
	    ! restricting the roughness length to a minimum
        grid%roughness_length(ji,jk) = max(0.0001_wp,grid%roughness_length(ji,jk))
      
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! Vertical grid
    ! -------------
    
    ! calculating the vertical positions of the scalar points
    ! the heights are defined according to k = A_k + B_k*surface with A_0 = toa, A_{NO_OF_LEVELS} = 0, B_0 = 0, B_{NO_OF_LEVELS} = 1
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,z_rel,sigma_z,A,B,vertical_vector_pre,max_oro)
    do ji=1,nlins
      do jk=1,ncols
        ! filling up vertical_vector_pre
        do jl=1,nlays+1
          z_rel = 1._wp-(jl-1._wp)/nlays ! z/toa
          sigma_z = z_rel**sigma
          A = sigma_z*toa ! the height without orography
          ! B corrects for orography
          if (jl>=nlays-nlays_oro+1._wp) then
            B = (jl-(nlays-nlays_oro+1._wp))/nlays_oro
          else
            B = 0._wp
          endif
          vertical_vector_pre(jl)=A+B*grid%z_w(ji,jk,nlays+1)
        enddo
        
        ! doing a check
        if (ji==1 .and. jk==1) then
          max_oro = maxval(grid%z_w(:,:,nlays+1))
          if (max_oro >= vertical_vector_pre(nlays-nlays_oro+1)) then
            write(*,*) "Maximum of orography larger or equal to the height of the lowest flat level."
            write(*,*) "Aborting."
            call exit(1)
          endif
        endif
        
        ! placing the scalar points in the middle between the preliminary values of the adjacent levels
        do jl=1,nlays
          grid%z_scalar(ji,jk,jl) = 0.5_wp*(vertical_vector_pre(jl) + vertical_vector_pre(jl+1))
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the height of the u-vector points
    ! inner domain
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=2,ncols
      grid%z_u(:,jk,:) = 0.5_wp*(grid%z_scalar(:,jk-1,:) + grid%z_scalar(:,jk,:))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! boundaries
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%z_u(:,1,:) = 0.5_wp*(grid%z_scalar(:,1,:) + grid%z_scalar(:,ncols,:))
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%z_u(:,ncols+1,:) = grid%z_u(:,1,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    else
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%z_u(:,1,:) = 2._wp*grid%z_scalar(:,1,:) - grid%z_scalar(:,2,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%z_u(:,ncols+1,:) = 2._wp*grid%z_scalar(:,ncols,:) - grid%z_scalar(:,ncols-1,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! setting the height of the v-vector points
    ! inner domain
    do ji=2,nlins
      grid%z_v(ji,:,:) = 0.5_wp*(grid%z_scalar(ji-1,:,:) + grid%z_scalar(ji,:,:))
    enddo
    ! boundaries
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%z_v(1,:,:) = 0.5_wp*(grid%z_scalar(1,:,:) + grid%z_scalar(nlins,:,:))
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%z_v(nlins+1,:,:) = grid%z_v(nlins,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    else
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%z_v(1,:,:) = 2._wp*grid%z_scalar(1,:,:) - grid%z_scalar(2,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%z_v(nlins+1,:,:) = 2._wp*grid%z_scalar(nlins,:,:) - grid%z_scalar(nlins-1,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! setting dx
    !$OMP PARALLEL
    !$OMP WORKSHARE
    grid%dx = dx*(re + grid%z_u)/re
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! setting dy
    !$OMP PARALLEL
    !$OMP WORKSHARE
    grid%dy = dy*(re + grid%z_v)/re
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    
    ! plane geometry grid distances
    if (lplane) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%dx = dx
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%dy = dy
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! setting dx of the dual grid
    !$OMP PARALLEL
    !$OMP WORKSHARE
    grid%dx_dual = dx*(re + grid%z_v)/re
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! setting dy of the dual grid
    !$OMP PARALLEL
    !$OMP WORKSHARE
    grid%dy_dual = dy*(re + grid%z_u)/re
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    
    ! plane geometry dual grid distances
    if (lplane) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%dx_dual = dx
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%dy_dual = dy
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! calculating the coordinate slopes
    call grad_hor_cov(grid%z_scalar,grid%slope_x,grid%slope_y,grid)
    
    ! setting the z coordinates of the vertical vector points
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl)
    do jl=1,nlays
      if (jl==1) then
        grid%z_w(:,:,jl) = toa
      else
        grid%z_w(:,:,jl) = 0.5_wp*(grid%z_scalar(:,:,jl-1) + grid%z_scalar(:,:,jl))
      endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the vertical distances between the scalar data points
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl)
    do jl=1,nlays+1
      if (jl==1) then
        grid%dz(:,:,jl) = 2._wp*(toa - grid%z_scalar(:,:,jl))
      elseif (jl==nlays+1) then
        grid%dz(:,:,jl) = 2._wp*(grid%z_scalar(:,:,jl-1) - grid%z_w(:,:,jl))
      else
        grid%dz(:,:,jl) = grid%z_scalar(:,:,jl-1) - grid%z_scalar(:,:,jl)
      endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the horizontal areas at the surface
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins
      do jk=1,ncols
        grid%area_z(ji,jk,nlays+1) = patch_area(grid%lat_scalar(ji),dlon,dlat)*(re + grid%z_w(ji,jk,nlays+1))**2/re**2
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! setting the horizontal areas at the higher points (above the surface)
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl)
    do jl=1,nlays
      grid%area_z(:,:,jl) = grid%area_z(:,:,nlays+1)*(re + grid%z_w(:,:,jl))**2 &
      /(re + grid%z_w(:,:,nlays+1))**2
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! plane geometry horizontal areas
    if (lplane) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%area_z(:,:,:) = dx*dy
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! the mean velocity area can be set now
    grid%mean_velocity_area = 2._wp*sum(grid%area_z(:,:,nlays+1))/size(grid%area_z(:,:,nlays+1))
    
    ! setting the vertical areas in x-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,lower_z,upper_z,lower_length)
    do ji=1,nlins
      do jk=1,ncols+1
        do jl=1,nlays
          ! left boundary
          if (jk==1) then
            lower_z = grid%z_w(ji,1,jl+1)+0.5_wp*(grid%z_w(ji,1,jl+1)-grid%z_w(ji,2,jl+1))
            upper_z = grid%z_w(ji,1,jl)+0.5_wp*(grid%z_w(ji,1,jl)-grid%z_w(ji,2,jl))
          ! right boundary
          elseif (jk==ncols+1) then
            lower_z = grid%z_w(ji,ncols,jl+1)+0.5_wp*(grid%z_w(ji,ncols,jl+1)-grid%z_w(ji,ncols-1,jl+1))
            upper_z = grid%z_w(ji,ncols,jl)+0.5_wp*(grid%z_w(ji,ncols,jl)-grid%z_w(ji,ncols-1,jl))
          ! inner domain
          else
            lower_z = 0.5_wp*(grid%z_w(ji,jk-1,jl+1) + grid%z_w(ji,jk,jl+1))
            upper_z = 0.5_wp*(grid%z_w(ji,jk-1,jl) + grid%z_w(ji,jk,jl))
          endif
          lower_length = dy*(re+lower_z)/re
          grid%area_x(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the vertical areas in y-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,lower_z,upper_z,lower_length)
    do ji=1,nlins+1
      do jk=1,ncols
        do jl=1,nlays
          ! upper boundary
          if (ji==1) then
            lower_z = grid%z_w(1,jk,jl+1)+0.5_wp*(grid%z_w(1,jk,jl+1)-grid%z_w(2,jk,jl+1))
            upper_z = grid%z_w(1,jk,jl)+0.5_wp*(grid%z_w(1,jk,jl)-grid%z_w(2,jk,jl))
          ! lower boundary
          elseif (ji==nlins+1) then
            lower_z = grid%z_w(nlins,jk,jl+1)+0.5_wp*(grid%z_w(nlins,jk,jl+1)-grid%z_w(nlins-1,jk,jl+1))
            upper_z = grid%z_w(nlins,jk,jl)+0.5_wp*(grid%z_w(nlins,jk,jl)-grid%z_w(nlins-1,jk,jl))
          ! inner domain
          else
            lower_z = 0.5_wp*(grid%z_w(ji-1,jk,jl+1) + grid%z_w(ji,jk,jl+1))
            upper_z = 0.5_wp*(grid%z_w(ji-1,jk,jl) + grid%z_w(ji,jk,jl))
          endif
          if (ji==nlins+1) then
            lower_length = dx*cos(grid%lat_scalar(ji-1)-0.5_wp*dlat)*(re+lower_z)/re
          else
            lower_length = dx*cos(grid%lat_scalar(ji)+0.5_wp*dlat)*(re+lower_z)/re
          endif
          grid%area_y(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the horizontal dual areas
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins+1
      do jk=1,ncols+1
        do jl=1,nlays
          ! setting the vertical position of the areas
          if (jk==1) then
            grid%z_area_dual_z(ji,jk,jl) = grid%z_v(ji,1,jl) + 0.5_wp*(grid%z_v(ji,1,jl)-grid%z_v(ji,2,jl))
          elseif (jk==ncols+1) then
            grid%z_area_dual_z(ji,jk,jl) = grid%z_v(ji,ncols,jl) + 0.5_wp*(grid%z_v(ji,ncols,jl)-grid%z_v(ji,ncols-1,jl))
          else
            grid%z_area_dual_z(ji,jk,jl) = 0.5_wp*(grid%z_v(ji,jk-1,jl)+grid%z_v(ji,jk,jl))
          endif
          ! setting the area itself
          if (ji==nlins+1) then
            grid%area_dual_z(ji,jk,jl) = patch_area(grid%lat_scalar(ji-1) - 0.5_wp*dlat,dlon,dlat) &
            *(re + grid%z_area_dual_z(ji,jk,jl))**2/re**2
          else
            grid%area_dual_z(ji,jk,jl) = patch_area(grid%lat_scalar(ji) + 0.5_wp*dlat,dlon,dlat) &
            *(re + grid%z_area_dual_z(ji,jk,jl))**2/re**2
          endif
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the vertical dual areas in x-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,lower_z,lower_length,upper_z)
    do ji=1,nlins+1
      do jk=1,ncols
        do jl=1,nlays+1
          if (jl==nlays+1) then
            if (ji==1) then
              lower_z = grid%z_w(1,jk,jl) + 0.5_wp*(grid%z_w(1,jk,jl)-grid%z_w(2,jk,jl))
            elseif (ji==nlins+1) then
              lower_z = grid%z_w(nlins,jk,jl) + 0.5_wp*(grid%z_w(nlins,jk,jl)-grid%z_w(nlins-1,jk,jl))
            else
              lower_z = 0.5_wp*(grid%z_w(ji-1,jk,jl) + grid%z_w(ji,jk,jl))
            endif
            lower_length = grid%dy(ji,jk,nlays)*(re + lower_z)/(re + grid%z_v(ji,jk,nlays))
          else
            lower_z = grid%z_v(ji,jk,jl)
            lower_length = grid%dy(ji,jk,jl)
          endif
          if (jl==1) then
            if (ji==1) then
              upper_z = grid%z_w(1,jk,jl) + 0.5_wp*(grid%z_w(1,jk,jl)-grid%z_w(2,jk,jl))
            elseif (ji==nlins+1) then
              upper_z = grid%z_w(nlins,jk,jl) + 0.5_wp*(grid%z_w(nlins,jk,jl)-grid%z_w(nlins-1,jk,jl))
            else
              upper_z = 0.5_wp*(grid%z_w(ji-1,jk,jl) + grid%z_w(ji,jk,jl))
            endif
          else
            upper_z = grid%z_v(ji,jk,jl-1)
          endif
          grid%area_dual_x(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the vertical dual areas in y-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,lower_z,lower_length,upper_z)
    do ji=1,nlins
      do jk=1,ncols+1
        do jl=1,nlays+1
          if (jl==nlays+1) then
            if (jk==1) then
              lower_z = grid%z_w(ji,1,jl) + 0.5_wp*(grid%z_w(ji,2,jl)-grid%z_w(ji,1,jl))
            elseif (jk==ncols+1) then
              lower_z = grid%z_w(ji,ncols,jl) + 0.5_wp*(grid%z_w(ji,ncols,jl)-grid%z_w(ji,ncols-1,jl))
            else
              lower_z = 0.5_wp*(grid%z_w(ji,jk-1,jl) + grid%z_w(ji,jk,jl))
            endif
            lower_length = grid%dx(ji,jk,nlays)*(re + lower_z)/(re + grid%z_u(ji,jk,nlays))
          else
            lower_z = grid%z_u(ji,jk,jl)
            lower_length = grid%dx(ji,jk,jl)
          endif
          if (jl==1) then
            if (jk==1) then
              upper_z = grid%z_w(ji,1,jl) + 0.5_wp*(grid%z_w(ji,2,jl)-grid%z_w(ji,1,jl))
            elseif (jk==ncols+1) then
              upper_z = grid%z_w(ji,ncols,jl) + 0.5_wp*(grid%z_w(ji,ncols,jl)-grid%z_w(ji,ncols-1,jl))
            else
              upper_z = 0.5_wp*(grid%z_w(ji,jk-1,jl) + grid%z_w(ji,jk,jl))
            endif
          else
            upper_z = grid%z_u(ji,jk,jl-1)
          endif
          grid%area_dual_y(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! setting the volume of the grid boxes
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji)
    do jl=1,nlays
      grid%volume(:,:,jl) = 1._wp/3._wp*((re + grid%z_w(:,:,jl))**3 - (re + grid%z_w(:,:,jl+1))**3) &
      /(re + grid%z_w(:,:,jl+1))**2*grid%area_z(:,:,jl+1)
      ! plane geometry
      if (lplane) then
        grid%volume(:,:,jl) = grid%area_z(:,:,jl+1)*(grid%z_w(:,:,jl) - grid%z_w(:,:,jl+1))
      endif
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! Derived quantities
    ! ------------------
    
    ! setting the inner product weights
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          grid%inner_product_weights(ji,jk,jl,1) = grid%area_x(ji,jk+1,jl)*grid%dx(ji,jk+1,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,2) = grid%area_y(ji,jk,jl)*grid%dy(ji,jk,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,3) = grid%area_x(ji,jk,jl)*grid%dx(ji,jk,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,4) = grid%area_y(ji+1,jk,jl)*grid%dy(ji+1,jk,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,5) = grid%area_z(ji,jk,jl)*grid%dz(ji,jk,jl)/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,6) = grid%area_z(ji,jk,jl+1)*grid%dz(ji,jk,jl+1)/(2._wp*grid%volume(ji,jk,jl))
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the TRSK weights
    ! u
    do ji=1,nlins
      base_area = patch_area(grid%lat_scalar(ji),dlon,dlat)
      grid%trsk_weights_u(ji,1) = (0.5_wp - patch_area(grid%lat_scalar(ji)+0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area) &
      *dx*cos(grid%lat_scalar(ji)+0.5_wp*dlat)/(dx*cos(grid%lat_scalar(ji)))
      grid%trsk_weights_u(ji,2) = -(0.5_wp - patch_area(grid%lat_scalar(ji)+0.25_wp*dlat,dlon,0.5_wp*dlat)/base_area) &
      *dy/(dx*cos(grid%lat_scalar(ji)))
      grid%trsk_weights_u(ji,3) = -(0.5_wp - (patch_area(grid%lat_scalar(ji)+0.25_wp*dlat,dlon,0.5_wp*dlat) &
      + patch_area(grid%lat_scalar(ji)-0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat))/base_area) &
      *dx*cos(grid%lat_scalar(ji)-0.5_wp*dlat)/(dx*cos(grid%lat_scalar(ji)))
      grid%trsk_weights_u(ji,4) = grid%trsk_weights_u(ji,3)
      grid%trsk_weights_u(ji,5) = -grid%trsk_weights_u(ji,2)
      grid%trsk_weights_u(ji,6) = grid%trsk_weights_u(ji,1)
    enddo
    ! v
    do ji=1,nlins+1
      if (ji==nlins+1) then
        lat_lower_center = grid%lat_scalar(nlins)-dlat
      else
        lat_lower_center = grid%lat_scalar(ji)
      endif
      base_area = patch_area(lat_lower_center,dlon,dlat)
      grid%trsk_weights_v(ji,1) = -(0.5_wp - patch_area(lat_lower_center+0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area)
      grid%trsk_weights_v(ji,2) = grid%trsk_weights_v(ji,1)
      if (ji==1) then
        lat_upper_center = grid%lat_scalar(1)+dlat
      else
        lat_upper_center = grid%lat_scalar(ji-1)
      endif
      base_area = patch_area(lat_upper_center,dlon,dlat)
      grid%trsk_weights_v(ji,3) = -(0.5_wp - patch_area(lat_upper_center-0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area)
      grid%trsk_weights_v(ji,4) = grid%trsk_weights_v(ji,3)
    enddo
    
    ! Soil grid
    ! ---------
    
    sigma_soil = 0.36_wp
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
    real(wp) :: temperature ! temperature at the respective grid point
    real(wp) :: pressure    ! pressure at the respective grid point
    real(wp) :: b,c         ! abbreviations needed for the hydrostatic initialization routine
    integer  :: ji,jk,jl    ! index variables
    
    ! integrating the hydrostatic background state according to the given temperature profile and pressure in the lowest layer
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,b,c,temperature,pressure)
    do ji=1,nlins
      do jk=1,ncols
        ! integrating from bottom to top
        do jl=nlays,1,-1
          temperature = bg_temp(grid%z_scalar(ji,jk,jl))
          ! lowest layer
          if (jl==nlays) then
            pressure = bg_pres(grid%z_scalar(ji,jk,jl))
            grid%exner_bg(ji,jk,jl) = (pressure/p_0)**(specific_gas_constants(0)/spec_heat_capacities_p_gas(0))
            grid%theta_bg(ji,jk,jl) = temperature/grid%exner_bg(ji,jk,jl)
          ! other layers
          else
            ! solving a quadratic equation for the Exner pressure
            b = -0.5_wp*grid%exner_bg(ji,jk,jl+1)/bg_temp(grid%z_scalar(ji,jk,jl+1)) &
            *(temperature - bg_temp(grid%z_scalar(ji,jk,jl+1)) + 2.0_wp/ &
            spec_heat_capacities_p_gas(0)*(geopot(grid%z_scalar(ji,jk,jl)) - geopot(grid%z_scalar(ji,jk,jl+1))))
            c = grid%exner_bg(ji,jk,jl+1)**2*temperature/bg_temp(grid%z_scalar(ji,jk,jl+1))
            grid%exner_bg(ji,jk,jl) = b+sqrt(b**2+c)
            grid%theta_bg(ji,jk,jl) = temperature/grid%exner_bg(ji,jk,jl)
          endif
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! calculating the gradient of the background Exner pressure (only needs to be done once)
    call grad(grid%exner_bg,grid%exner_bg_grad_u,grid%exner_bg_grad_v,grid%exner_bg_grad_w,grid)
  
  end subroutine bg_setup
  
  subroutine set_orography(grid)
  
    ! This subroutine interpolates the real orography from ETOPO1.
    
    type(t_grid), intent(inout) :: grid
    
    ! local variables
    integer               :: ncid                   ! ID of the NetCDF file
    character(len=64)     :: filename               ! input filename
    integer               :: lat_in_id              ! variable ID of the latitude vector
    integer               :: lon_in_id              ! variable ID of the longitude vector
    integer               :: z_in_id                ! variable ID of the orography
    integer               :: no_of_lat_points       ! number of latitude points of the input dataset
    integer               :: no_of_lon_points       ! number of longitude points of the input dataset
    real(wp), allocatable :: latitude_input(:)      ! latitudes of the input dataset
    real(wp), allocatable :: longitude_input(:)     ! longitudes of the input dataset
    real(wp), allocatable :: lat_distance_vector(:) ! latitudes distance vector
    real(wp), allocatable :: lon_distance_vector(:) ! longitudes distance vector
    integer, allocatable  :: z_input(:,:)           ! the input dataset
    integer               :: lat_index,lon_index    ! minimum distance indices
    integer               :: ji,jk,jl               ! loop indices
    
    no_of_lat_points = 10801
    no_of_lon_points = 21601
    
    ! the filename of the grid file including the relative path
    filename = "../../grids/phys_sfc_properties/" // trim(oro_raw_filename)
    
    ! opening the NetCDF file
    call nc_check(nf90_open(trim(filename),NF90_CLOBBER,ncid))
    
    ! reading the variable IDs
    call nc_check(nf90_inq_varid(ncid,"y",lat_in_id))
    call nc_check(nf90_inq_varid(ncid,"x",lon_in_id))
    call nc_check(nf90_inq_varid(ncid,"z",z_in_id))
    
    ! allocating memory for reading the grid files
    allocate(latitude_input(no_of_lat_points))
    allocate(longitude_input(no_of_lon_points))
    allocate(z_input(no_of_lon_points,no_of_lat_points))
    
    ! reading the arrays
    call nc_check(nf90_get_var(ncid,lat_in_id,latitude_input))
    call nc_check(nf90_get_var(ncid,lon_in_id,longitude_input))
    call nc_check(nf90_get_var(ncid,z_in_id,z_input))
    
    ! closing the NetCDF file
    call nc_check(nf90_close(ncid))
    
    ! allocating memory for the distance vectors
    allocate(lat_distance_vector(no_of_lat_points))
    allocate(lon_distance_vector(no_of_lon_points))

    ! setting the unfiltered orography
    do ji=1,nlins
      do jk=1,ncols
        ! default
        grid%z_w(ji,jk,nlays+1) = 0._wp
    
        do jl=1,no_of_lat_points
          lat_distance_vector(jl) = abs(2._wp*M_PI*latitude_input(jl)/360._wp - grid%lat_geo_scalar(ji,jk))
        enddo
    
        do jl=1,no_of_lon_points
          lon_distance_vector(jl) = abs(2._wp*M_PI*longitude_input(jl)/360._wp - grid%lon_geo_scalar(ji,jk))
        enddo
    
        lat_index = find_min_index(lat_distance_vector)
        lon_index = find_min_index(lon_distance_vector)
          
        grid%z_w(ji,jk,nlays+1) = real(z_input(lon_index,lat_index),wp)

        ! check
        if (grid%z_w(ji,jk,nlays+1)<-382._wp .or. grid%z_w(ji,jk,nlays+1)>8850._wp) then
          write(*,*) "Warning: orography value out of usual range."
        endif
      
      enddo
    enddo

    ! freeing the memory
    deallocate(lat_distance_vector)
    deallocate(lon_distance_vector)
    deallocate(z_input)
    deallocate(latitude_input)
    deallocate(longitude_input)
  
  end subroutine set_orography
  
  subroutine smooth_hor_scalar(array)
  
    ! This subroutine smoothes a scalar field on one layer.
    
    real(wp), intent(inout) :: array(nlins,ncols)
    
    ! local variables
    real(wp) :: original_array(nlins,ncols) ! the unsmoothed input array
    integer  :: ji,jk                       ! loop indices
    
    ! copying the original array
    original_array = array
    
    ! inner domin
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=2,nlins-1
      do jk=2,ncols-1
        array(ji,jk) = sum(original_array(ji-1:ji+1,jk-1:jk+1))/9._wp
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! corners
    array(1,1) = sum(original_array(1:2,1:2))/4._wp
    array(1,ncols) = sum(original_array(1:2,ncols-1:ncols))/4._wp
    array(nlins,1) = sum(original_array(nlins-1:nlins,1:2))/4._wp
    array(nlins,ncols) = sum(original_array(nlins-1:nlins,ncols-1:ncols))/4._wp
  
    ! boundaries
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=2,ncols-1
      array(1,jk) = sum(original_array(1:2,jk-1:jk+1))/6._wp
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=2,ncols-1
      array(nlins,jk) = sum(original_array(nlins-1:nlins,jk-1:jk+1))/6._wp
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji)
    do ji=2,nlins-1
      array(ji,1) = sum(original_array(ji-1:ji+1,1:2))/6._wp
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji)
    do ji=2,nlins-1
      array(ji,ncols) = sum(original_array(ji-1:ji+1,ncols-1:ncols))/6._wp
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine smooth_hor_scalar
  
  function find_min_index(input_vector)
  
    ! This function finds the index where a vector assumes its minimum.
    
    real(wp) :: input_vector(:)
    integer  :: find_min_index ! the result
    
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
  
    ! input
    real(wp) :: center_lat  ! latitude at the center of the patch
    real(wp) :: dx_as_angle ! delta x as angle
    real(wp) :: dy_as_angle ! delta y as angle
    ! output
    real(wp) :: patch_area  ! the result
  
    ! computing the result
    patch_area = re**2*dx_as_angle*(sin(center_lat + 0.5_wp*dy_as_angle) - sin(center_lat - 0.5_wp*dy_as_angle))
    
    ! plane geometry
    if (lplane) then
      patch_area = re**2*dx_as_angle*dy_as_angle
    endif
  
  end function patch_area
  
  function vertical_face_area(lower_z,upper_z,lower_length)
  
    ! This function calculates the area of a vertical face.
    
    ! input
    real(wp) :: lower_z            ! geometric height of the lower boundary of the face
    real(wp) :: upper_z            ! geometric height of the upper boundary of the face
    real(wp) :: lower_length       ! length of the lower boundary of the face
    ! output
    real(wp) :: vertical_face_area ! the result
    
    vertical_face_area = 0.5_wp*lower_length*(re + upper_z + re + lower_z) &
    /(re + lower_z)*(upper_z - lower_z)
    
    ! plane geometry
    if (lplane) then
      vertical_face_area = lower_length*(upper_z - lower_z)
    endif
  
  end function vertical_face_area
  
  function vegetation_height_ideal(latitude,oro)
  
    ! This function calculates a latitude- and height-dependant idealized vegetation height.

    ! input arguments
    real(wp) :: latitude                ! latitude of this point
    real(wp) :: oro                     ! height of the terrain at this point
    ! output
    real(wp) :: vegetation_height_ideal ! the result
    
    ! local variables
    real(wp) :: vegetation_height_equator ! the vegetation height at the equator

    ! setting the vegetation height at the equator
    vegetation_height_equator = 20._wp

    vegetation_height_ideal = vegetation_height_equator*cos(latitude)*exp(-oro/1500._wp)

  end function vegetation_height_ideal
  
  subroutine find_global_normal(lat,lon,r)

    ! This subroutine calculates the Cartesian normal vector of a point given its geographical coordinates.
    
    ! input arguments and output
    real(wp), intent(in)  :: lat
    real(wp), intent(in)  :: lon
    real(wp), intent(out) :: r(3) ! positional vector (result)

    r(1) = cos(lat)*cos(lon)
    r(2) = cos(lat)*sin(lon)
    r(3) = sin(lat)

  end subroutine find_global_normal
  
  subroutine find_geos(r,lat_out,lon_out)

    ! This subroutine calculates the geographical coordinates of a point given its Cartesian coordinates.

    ! input arguments and output
    real(wp), intent(in)  :: r(3)    ! positional vector
    real(wp), intent(out) :: lat_out
    real(wp), intent(out) :: lon_out
    
    lat_out = asin(r(3)/sqrt(r(1)**2+r(2)**2+r(3)**2))
    lon_out = atan2(r(2),r(1))

  end subroutine find_geos

  subroutine calc_local_i(lon,result_vec)

    ! This function calculates the local eastward basis vector.
    
    ! input arguments and output
    real(wp), intent(in)  :: lon
    real(wp), intent(out) :: result_vec(3)
    
    result_vec(1) = -sin(lon)
    result_vec(2) = cos(lon)
    result_vec(3) = 0._wp
  
  end subroutine calc_local_i

  subroutine calc_local_j(lat,lon,result_vec)

    ! This subroutine calculates the local northward basis vector.
    
    ! input arguments and output
    real(wp), intent(in)  :: lat
    real(wp), intent(in)  :: lon
    real(wp), intent(out) :: result_vec(3)
    
    result_vec(1) = -sin(lat)*cos(lon)
    result_vec(2) = -sin(lat)*sin(lon)
    result_vec(3) = cos(lat)

  end subroutine calc_local_j

end module grid_generator






