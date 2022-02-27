! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

! This file contains the calculation of the grid properties as well as some type definitions.

module grid_generator

  use definitions,        only: wp,t_grid
  use run_nml,            only: nlins,ncols,nlays,dy,dx,toa,nlays_oro,sigma,scenario,lat_center_deg, &
                                lon_center_deg,x_dir_deg
  use constants,          only: re,density_water,T_0,M_PI,p_0,omega,gravity,p_0_standard, &
                                lapse_rate,surface_temp,tropo_height,inv_height,t_grad_inv
  use surface_nml,        only: nsoillays,orography_id
  use gradient_operators, only: grad,grad_hor_cov
  use dictionary,         only: specific_gas_constants,spec_heat_capacities_p_gas
  use io_nml,             only: lwrite_grid,lread_grid
  use read_write_grid,    only: write_grid,read_grid
  use set_initial_state,  only: bg_temp,bg_pres,geopot
  use bc_nml,             only: lperiodic

  implicit none
  
  private
  
  public :: grid_setup
  public :: bg_setup
  
  interface
    ! interface to the geos95 library
    real(C_DOUBLE) function calculate_distance_h(latitude_a,longitude_a,latitude_b,longitude_b,radius) &
    bind(c, name = "calculate_distance_h")
      use, intrinsic::iso_c_binding
      implicit none
      real(C_DOUBLE), value :: latitude_a
      real(C_DOUBLE), value :: longitude_a
      real(C_DOUBLE), value :: latitude_b
      real(C_DOUBLE), value :: longitude_b
      real(C_DOUBLE), value :: radius
    end function calculate_distance_h
  end interface
  
  contains
  
  subroutine grid_setup(grid)
  
    type(t_grid), intent(inout) :: grid ! the model grid
    ! local variables
    real(wp) :: lat_left_upper                 ! latitude coordinate of upper left corner
    real(wp) :: lon_left_upper                 ! longitude coordinate of upper left corner
    real(wp) :: dlat                           ! mesh size in y direction as angle
    real(wp) :: dlon                           ! mesh size in x direction as angle
    integer  :: ji,jk,jl                       ! loop indices
    real(wp) :: max_oro                        ! variable for orography check
    real(wp) :: A                              ! variable for calculating the vertical grid
    real(wp) :: B                              ! variable for calculating the vertical grid
    real(wp) :: sigma_z                        ! variable for calculating the vertical grid
    real(wp) :: rel                          ! variable for calculating the vertical grid
    real(wp) :: vertical_vector_pre(nlays+1) ! variable for calculating the vertical grid
    real(wp) :: base_area                      ! variable for calculating the vertical grid
    real(wp) :: lower_z,upper_z,lower_length   ! variables needed for area calculations
    real(wp) :: height_mountain                ! height of Gaussian mountain (needed for test case)
    real(wp) :: sigma_mountain                 ! standard deviation of Gaussian mountain (needed for test case)
    real(wp) :: x_coord                        ! help variable needed for the Schär test
    real(wp) :: rescale_factor                 ! soil grid rescaling factor
    real(wp) :: sigma_soil                     ! sigma of the soil grid
    real(wp) :: density_soil                   ! typical density of soil
    real(wp) :: c_p_soil                       ! typical c_p of soil
    real(wp) :: c_p_water                      ! typical c_p of water
    
    ! setting the center and direction of the grid
    grid%lat_center = 2._wp*M_PI*lat_center_deg/360._wp
    grid%lon_center = 2._wp*M_PI*lon_center_deg/360._wp
    grid%x_dir_center = 2._wp*M_PI*x_dir_deg/360._wp
    
    ! setting the latitude and longitude coordinates of the scalar grid points
    ! setting the dy of the model grid
    dlat = dy/re
    dlon = dx/re
    lat_left_upper = (nlins-1._wp)/2._wp*dlat
    lon_left_upper = -(ncols-1._wp)/2._wp*dlon
    do ji=1,nlins
      grid%lat_scalar(ji) = lat_left_upper - dlat*(ji-1._wp)
    enddo
    do jk=1,ncols
      grid%lon_scalar(jk) = lon_left_upper + dlon*(jk-1._wp)
    enddo
    
    ! setting the Coriolis vector at the grid points
    do ji=1,nlins+1
      do jk=1,ncols
        grid%fvec_x(ji,jk) = 0._wp
      enddo
    enddo
    do ji=1,nlins
      do jk=1,ncols+1
        grid%fvec_y(ji,jk) = omega
      enddo
    enddo
    do ji=1,nlins+1
      do jk=1,ncols+1
        grid%fvec_z(ji,jk) = omega
      enddo
    enddo
    
    ! setting the physical surface properties, including orography
    select case (orography_id)
    
      ! no orography
      case(0)
        grid%geo_w(:,:,nlays+1) = 0._wp
      
      ! real orography
      case(1)
        ! reading the grid
        if (lread_grid) then
          call read_grid(grid)
        else
        
          ! real orography is not yet implemented
          grid%geo_w(:,:,nlays+1) = 0._wp
        
          ! idealized soil properties are being set here if no grid file shall be read in
          density_soil = 1442._wp
          c_p_soil = 830._wp
          c_p_water = 4184._wp
    
          do ji=1,nlins
            do jk=1,ncols
		
              grid%t_const_soil(ji,jk) = T_0 + 25._wp*cos(2._wp*grid%lat_geo_scalar(ji,jk))
            
              ! albedo of water
              grid%sfc_albedo(ji,jk) = 0.06_wp

              ! for water, the roughness_length is set to some sea-typical value, will not be used anyway
              grid%roughness_length(ji,jk) = 0.08_wp
		
              ! will also not be used for water
              grid%sfc_rho_c(ji,jk) = density_water*c_p_water
		
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
                grid%roughness_length(ji,jk) = vegetation_height_ideal(grid%lat_geo_scalar(ji,jk),grid%geo_w(ji,jk,nlays+1))/8._wp
          
              endif
		
		      ! restricting the roughness length to a minimum
              grid%roughness_length(ji,jk) = max(0.0001_wp,grid%roughness_length(ji,jk))
      
            enddo
          enddo
        endif
      
      ! Schaer orography
      case(2)
        height_mountain = 250._wp
        sigma_mountain = 5000._wp/sqrt(2._wp)
        do ji=1,nlins
          do jk=1,ncols
            x_coord = calculate_distance_h(grid%lat_scalar(ji),grid%lon_scalar(jk),0._wp,0._wp,re)
            grid%geo_w(ji,jk,nlays+1) = height_mountain*exp(-x_coord**2/(2._wp*sigma_mountain**2))*cos(M_PI*x_coord/4000._wp)**2
          enddo
        enddo
    
    endselect
  
    ! calculating the vertical positions of the scalar points
    ! the heights are defined according to k = A_k + B_k*surface with A_0 = toa, A_{NO_OF_LEVELS} = 0, B_0 = 0, B_{NO_OF_LEVELS} = 1
    do ji=1,nlins
      do jk=1,ncols
        ! filling up vertical_vector_pre
        do jl=1,nlays+1
          rel = 1._wp-(jl-1._wp)/nlays ! z/toa
          sigma_z = rel**sigma
          A = sigma_z*toa ! the height without orography
          ! B corrects for orography
          if (jl >= nlays-nlays_oro+1._wp) then
            B = (jl-(nlays-nlays_oro+1._wp))/nlays_oro
          else
            B = 0
          endif
          vertical_vector_pre(jl)=A+B*grid%geo_w(ji,jk,nlays+1)
        enddo
        
        ! doing a check
        if (ji == 1 .and. jk == 1) then
          max_oro = maxval(grid%geo_w(:,:,nlays+1))
          if (max_oro >= vertical_vector_pre(max(nlays-nlays_oro,1))) then
            write(*,*) "Maximum of orography larger or equal to the height of the lowest flat level."
            write(*,*) "Aborting."
            call exit(1)
          endif
        endif
        
        ! placing the scalar points in the middle between the preliminary values of the adjacent levels
        do jl=1,nlays
          grid%geo_scal(ji,jk,jl) = 0.5_wp*(vertical_vector_pre(jl) + vertical_vector_pre(jl+1))
        enddo
      enddo
    enddo
    
    ! setting the height of the u-vector points
    ! inner domain
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=2,ncols
      grid%geo_u(:,jk,:) = 0.5_wp*(grid%geo_scal(:,jk-1,:) + grid%geo_scal(:,jk,:))
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! boundaries
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%geo_u(:,1,:) = 0.5_wp*(grid%geo_scal(:,1,:) + grid%geo_scal(:,ncols,:))
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%geo_u(:,ncols+1,:) = grid%geo_u(:,1,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    else
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%geo_u(:,1,:) = 2._wp*grid%geo_scal(:,1,:) - grid%geo_scal(:,2,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%geo_u(:,ncols+1,:) = 2._wp*grid%geo_scal(:,ncols,:) - grid%geo_scal(:,ncols-1,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! setting the height of the v-vector points
    ! inner domain
    do ji=2,nlins
      grid%geo_v(ji,:,:) = 0.5_wp*(grid%geo_scal(ji-1,:,:) + grid%geo_scal(ji,:,:))
    enddo
    ! boundaries
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%geo_v(1,:,:) = 0.5_wp*(grid%geo_scal(1,:,:) + grid%geo_scal(nlins,:,:))
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%geo_v(nlins+1,:,:) = grid%geo_v(nlins,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    else
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%geo_v(1,:,:) = 2._wp*grid%geo_scal(1,:,:) - grid%geo_scal(2,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      grid%geo_v(nlins+1,:,:) = 2._wp*grid%geo_scal(nlins,:,:) - grid%geo_scal(nlins-1,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! setting dx
    !$OMP PARALLEL
    !$OMP WORKSHARE
    grid%dx = dx*(re + grid%geo_u)/re
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! setting dy
    !$OMP PARALLEL
    !$OMP WORKSHARE
    grid%dy = dy*(re + grid%geo_v)/re
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    
    ! calculating the coordinate slopes
    call grad_hor_cov(grid%geo_scal,grid%slope_x,grid%slope_y,grid)
    
    ! setting the z coordinates of the vertical vector points
    do jl=1,nlays
      if (jl==1) then
        grid%geo_w(:,:,jl) = toa
      else
        grid%geo_w(:,:,jl) = 0.5_wp*(grid%geo_scal(:,:,jl-1) + grid%geo_scal(:,:,jl))
      endif
    enddo
    
    ! setting the vertical distances between the scalar data points
    do jl=1,nlays+1
      if (jl==1) then
        grid%dz(:,:,jl) = 2._wp*(toa - grid%geo_scal(:,:,jl))
      elseif (jl==nlays+1) then
        grid%dz(:,:,jl) = 2._wp*(grid%geo_scal(:,:,jl-1) - grid%geo_w(:,:,jl))
      else
        grid%dz(:,:,jl) = grid%geo_scal(:,:,jl-1) - grid%geo_scal(:,:,jl)
      endif
    enddo
    
    ! setting the horizontal areas at the surface
    do ji=1,nlins
      do jk=1,ncols
        grid%area_z(ji,jk,nlays+1) = patch_area(grid%lat_scalar(ji),dlon,dlat)*(re + grid%geo_w(ji,jk,nlays+1))**2/re**2
      enddo
    enddo

    ! setting the horizontal areas at the higher points (above the surface)
    do jl=1,nlays
      grid%area_z(:,:,jl) = grid%area_z(:,:,nlays+1)*(re + grid%geo_w(:,:,jl))**2 &
      /(re + grid%geo_w(:,:,nlays+1))**2
    enddo
    
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
            lower_z = grid%geo_w(ji,1,jl+1)+0.5_wp*(grid%geo_w(ji,1,jl+1)-grid%geo_w(ji,2,jl+1))
            upper_z = grid%geo_w(ji,1,jl)+0.5_wp*(grid%geo_w(ji,1,jl)-grid%geo_w(ji,2,jl))
          ! right boundary
          elseif (jk==ncols+1) then
            lower_z = grid%geo_w(ji,ncols,jl+1)+0.5_wp*(grid%geo_w(ji,ncols,jl+1)-grid%geo_w(ji,ncols-1,jl+1))
            upper_z = grid%geo_w(ji,ncols,jl)+0.5_wp*(grid%geo_w(ji,ncols,jl)-grid%geo_w(ji,ncols-1,jl))
          ! inner domain
          else
            lower_z = 0.5_wp*(grid%geo_w(ji,jk-1,jl+1) + grid%geo_w(ji,jk,jl+1))
            upper_z = 0.5_wp*(grid%geo_w(ji,jk-1,jl) + grid%geo_w(ji,jk,jl))
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
            lower_z = grid%geo_w(1,jk,jl+1)+0.5_wp*(grid%geo_w(1,jk,jl+1)-grid%geo_w(2,jk,jl+1))
            upper_z = grid%geo_w(1,jk,jl)+0.5_wp*(grid%geo_w(1,jk,jl)-grid%geo_w(2,jk,jl))
          ! lower boundary
          elseif (ji==nlins+1) then
            lower_z = grid%geo_w(nlins,jk,jl+1)+0.5_wp*(grid%geo_w(nlins,jk,jl+1)-grid%geo_w(nlins-1,jk,jl+1))
            upper_z = grid%geo_w(nlins,jk,jl)+0.5_wp*(grid%geo_w(nlins,jk,jl)-grid%geo_w(nlins-1,jk,jl))
          ! inner domain
          else
            lower_z = 0.5_wp*(grid%geo_w(ji-1,jk,jl+1) + grid%geo_w(ji,jk,jl+1))
            upper_z = 0.5_wp*(grid%geo_w(ji-1,jk,jl) + grid%geo_w(ji,jk,jl))
          endif
          lower_length = dx*cos(grid%lat_scalar(ji)+0.5_wp*dlat)*(re+lower_z)/re
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
          grid%geo_area_dual_z(ji,jk,jl) = 0.25_wp*(grid%geo_scal(ji,jk,jl)+grid%geo_scal(ji+1,jk,jl) &
          +grid%geo_scal(ji+1,jk+1,jl)+grid%geo_scal(ji,jk+1,jl))
          grid%area_dual_z(ji,jk,jl) = patch_area(grid%lat_scalar(ji) + 0.5_wp*dlat,dlon,dlat) &
          *(re + grid%geo_area_dual_z(ji,jk,jl))**2/re**2
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
            lower_z = 0.5_wp*(grid%geo_w(ji,jk+1,jl) + grid%geo_w(ji+1,jk+1,jl))
            lower_length = grid%dy(ji,jk+1,jl-1)*(re + lower_z)/ &
            (re + 0.5_wp*(grid%geo_scal(ji,jk+1,jl-1) + grid%geo_scal(ji+1,jk+1,jl-1)))
          else
            lower_z = 0.5_wp*(grid%geo_scal(ji,jk+1,jl) + grid%geo_scal(ji+1,jk+1,jl))
            lower_length = grid%dy(ji,jk+1,jl)
          endif
          if (jl==1) then
            upper_z = 0.5_wp*(grid%geo_w(ji,jk+1,jl) + grid%geo_w(ji+1,jk+1,jl))
          else
            upper_z = 0.5_wp*(grid%geo_scal(ji,jk+1,jl-1) + grid%geo_scal(ji+1,jk+1,jl-1))
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
            lower_z = 0.5_wp*(grid%geo_w(ji+1,jk,jl) + grid%geo_w(ji+1,jk+1,jl))
            lower_length = grid%dx(ji+1,jk,jl-1)*(re + lower_z)/ &
            (re + 0.5_wp*(grid%geo_scal(ji+1,jk,jl-1) + grid%geo_scal(ji+1,jk+1,jl-1)))
          else
            lower_z = 0.5_wp*(grid%geo_scal(ji+1,jk,jl) + grid%geo_scal(ji+1,jk+1,jl))
            lower_length = grid%dx(ji+1,jk,jl)
          endif
          if (jl==1) then
            upper_z = 0.5_wp*(grid%geo_w(ji+1,jk,jl) + grid%geo_w(ji+1,jk+1,jl))
          else
            upper_z = 0.5_wp*(grid%geo_scal(ji+1,jk,jl-1) + grid%geo_scal(ji+1,jk+1,jl-1))
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
      grid%volume(:,:,jl) = 1._wp/3._wp*((re + grid%geo_w(:,:,jl))**3 - (re + grid%geo_w(:,:,jl+1))**3) &
      /(re + grid%geo_w(:,:,jl+1))**2*grid%area_z(:,:,jl+1)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
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
        *dx*cos(0.5_wp*(grid%lat_scalar(ji)+grid%lat_scalar(ji+1)))/(dx*cos(grid%lat_scalar(ji)))
        grid%trsk_weights_u(ji,2) = -(0.5_wp - patch_area(grid%lat_scalar(ji)+0.25_wp*dlat,dlon,0.5_wp*dlat)/base_area) &
        *dy/(dx*cos(grid%lat_scalar(ji)))
        grid%trsk_weights_u(ji,3) = -(0.5_wp - (patch_area(grid%lat_scalar(ji)+0.25_wp*dlat,dlon,0.5_wp*dlat) &
        +patch_area(grid%lat_scalar(ji)-0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat))/base_area) &
        *dx*cos(0.5_wp*(grid%lat_scalar(ji)+grid%lat_scalar(ji)))/(dx*cos(grid%lat_scalar(ji)))
        grid%trsk_weights_u(ji,4) = grid%trsk_weights_u(ji,3)
        grid%trsk_weights_u(ji,5) = -grid%trsk_weights_u(ji,2)
        grid%trsk_weights_u(ji,6) = grid%trsk_weights_u(ji,1)
    enddo
    ! v
    do ji=1,nlins+1
        base_area = patch_area(grid%lat_scalar(ji),dlon,dlat)
        grid%trsk_weights_v(ji,1) = -(0.5_wp - patch_area(grid%lat_scalar(ji)+0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area)
        grid%trsk_weights_v(ji,2) = grid%trsk_weights_v(ji,1)
        base_area = patch_area(grid%lat_scalar(ji+1),dlon,dlat)
        grid%trsk_weights_v(ji,3) = -(0.5_wp - patch_area(grid%lat_scalar(ji+1)-0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area)
        grid%trsk_weights_v(ji,4) = grid%trsk_weights_v(ji,3)
    enddo
    
    ! soil grid
    sigma_soil = 0.36_wp
    grid%t_const = -10._wp

    ! the surface is always at zero
    grid%soil_interface(1) = 0._wp
    do jl=2,nsoillays+1
      grid%soil_interface(jl) = grid%soil_interface(jl-1) + sigma_soil**(nsoillays-jl)
    enddo
    
    rescale_factor = grid%t_const/grid%soil_interface(nsoillays+1)
    
    do jl=2,nsoillays+1
      grid%soil_interface(jl) = rescale_factor*grid%soil_interface(jl)
    enddo
    do jl=1,nsoillays
      grid%soil_center(jl) = 0.5_wp*(grid%soil_interface(jl) + grid%soil_interface(jl+1))
    enddo
    
    write(*,*) "Thickness of the uppermost soil layer: ", -grid%soil_interface(2), "m."
	
    ! writing the costly grid properties to a file if required by the user
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
    do ji=1,nlins
      do jk=1,ncols
        ! integrating from bottom to top
        do jl=nlays,1,-1
          temperature = bg_temp(grid%geo_scal(ji,jk,jl))
          ! lowest layer
          if (jl==nlays) then
            pressure = bg_pres(grid%geo_scal(ji,jk,jl))
            grid%exner_bg(ji,jk,jl) = (pressure/p_0)**(specific_gas_constants(0)/spec_heat_capacities_p_gas(0))
            grid%theta_bg(ji,jk,jl) = temperature/grid%exner_bg(ji,jk,jl)
          ! other layers
          else
            ! solving a quadratic equation for the Exner pressure
            b = -0.5_wp*grid%exner_bg(ji,jk,jl+1)/bg_temp(grid%geo_scal(ji,jk,jl+1)) &
            *(temperature - bg_temp(grid%geo_scal(ji,jk,jl+1)) + 2.0_wp/ &
            spec_heat_capacities_p_gas(0)*(geopot(grid%geo_scal(ji,jk,jl)) - geopot(grid%geo_scal(ji,jk,jl+1))))
            c = grid%exner_bg(ji,jk,jl+1)**2*temperature/bg_temp(grid%geo_scal(ji,jk,jl+1))
            grid%exner_bg(ji,jk,jl) = b+sqrt(b**2+c)
            grid%theta_bg(ji,jk,jl) = temperature/grid%exner_bg(ji,jk,jl)
          endif
        enddo
      enddo
    enddo
    
    ! calculating the gradient of the background Exner pressure (only needs to be done once)
    call grad(grid%exner_bg,grid%exner_bg_grad_u,grid%exner_bg_grad_v,grid%exner_bg_grad_w,grid)
  
  end subroutine bg_setup
  
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

  
end module grid_generator






