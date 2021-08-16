! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

! This file contains the calculation of the grid properties as well as some type definitions.

module grid_generator

  use definitions,        only: wp,t_grid
  use run_nml,            only: nlins,ncols,nlays,dy,dx,toa,nlays_oro,sigma,re,omega,p_0,gravity, &
                                lapse_rate,surface_temp,tropo_height,inv_height,t_grad_inv,p_0_standard, &
                                scenario
  use gradient_operators, only: grad_hor_cov_extended,grad
  use thermodynamics,     only: gas_constant_diagnostics,spec_heat_cap_diagnostics_p

  implicit none
  
  private
  
  public :: grid_setup
  public :: bg_setup
  public :: bg_temp
  public :: bg_pres
  public :: geopot
  
  interface
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
    real(wp) :: lat_left_lower  ! latitude coordinate of lower left corner
    real(wp) :: lon_left_lower  ! longitude coordinate of lower left corner
    real(wp) :: dlat            ! mesh size in y direction as angle
    real(wp) :: dlon            ! mesh size in x direction as angle
    integer  :: ji,jk,jl        ! loop indices
    real(wp) :: max_oro         ! variable for orography check
    real(wp) :: A               ! variable for calculating the vertical grid
    real(wp) :: B               ! variable for calculating the vertical grid
    real(wp) :: sigma_z         ! variable for calculating the vertical grid
    real(wp) :: z_rel           ! variable for calculating the vertical grid
    real(wp) :: z_vertical_vector_pre(nlays+1)
                                ! variable for calculating the vertical grid
    real(wp) :: base_area       ! variable for calculating the vertical grid
    real(wp) :: lower_z,upper_z,lower_length
                                ! variables needed for area calculations
    real(wp) :: height_mountain ! height of Gaussian mountain (needed for test case)
    real(wp) :: sigma_mountain  ! standard deviation of Gaussian mountain (needed for test case)
    
    ! setting the latitude and longitude coordinates of the scalar grid points
    ! setting the dy of the model grid
    dlat = dy/re
    dlon = dx/re
    lat_left_lower = -(nlins+2-1)/2*dlat
    lon_left_lower = -(ncols+2-1)/2*dlon
    do ji=1,nlins+2
      grid%lat_scalar(ji) = lat_left_lower + dlat*(ji - 1)
    enddo
    do ji=1,ncols+2
      grid%lon_scalar(ji) = lon_left_lower + dlon*(ji - 1)
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
    
    ! setting up the orography of the grid
    select case (trim(scenario))
    
      case("standard")
        do ji=1,nlins+2
          do jk=1,ncols+2
            grid%z_geo_w(ji,jk,nlays+1) = 0._wp
          enddo
        enddo

      case("resting_mountain","schaer")
        height_mountain = 100._wp
        sigma_mountain = 7000._wp
        do ji=1,nlins+2
          do jk=1,ncols+2
            grid%z_geo_w(ji,jk,nlays+1) = height_mountain*exp(-calculate_distance_h(grid%lat_scalar(ji),grid%lon_scalar(jk), &
            0._wp,0._wp,re)**2/(2._wp*sigma_mountain**2))
          enddo
        enddo
    
    endselect
  
    ! calculating the vertical positions of the scalar points
    ! the heights are defined according to z_k = A_k + B_k*z_surface with A_0 = toa, A_{NO_OF_LEVELS} = 0, B_0 = 0, B_{NO_OF_LEVELS} = 1
    do ji=1,nlins+2
      do jk=1,ncols+2
        ! filling up z_vertical_vector_pre
        do jl=1,nlays+1
          z_rel = 1._wp - (jl - 1._wp)/(nlays+1) ! z/toa
          sigma_z = z_rel**sigma
          A = sigma_z*toa; ! the height without orography
          ! B corrects for orography
          if (jl >= nlays - nlays_oro + 1._wp) then
            B = (jl - (nlays - nlays_oro + 1._wp))/nlays_oro
          else
            B = 0
          endif
          z_vertical_vector_pre(jl) = A + B*grid%z_geo_w(ji,jk,nlays+1)
        enddo
        
        ! doing a check
        if (ji == 1 .and. jk == 1) then
          max_oro = maxval(grid%z_geo_w(:,:,nlays+1))
          if (max_oro >= z_vertical_vector_pre(nlays - nlays_oro)) then
            write(*,*) "Maximum of orography larger or equal to the height of the lowest flat level."
            write(*,*) "Aborting."
            call exit(1)
          endif
        endif
        
        ! placing the scalar points in the middle between the preliminary values of the adjacent levels
        do jl=1,nlays
          grid%z_geo_scal(ji,jk,jl) = 0.5_wp*(z_vertical_vector_pre(jl) + z_vertical_vector_pre(jl + 1))
        enddo
      enddo
    enddo
    
    ! setting dx
    do ji=1,nlins+2
      do jk=1,ncols+1
        do jl=1,nlays
          grid%z_geo_u(ji,jk,jl) = 0.5_wp*(grid%z_geo_scal(ji,jk,jl) + grid%z_geo_scal(ji,jk+1,jl))
          grid%dx(ji,jk,jl) = dx*(re + grid%z_geo_u(ji,jk,jl))/re
        enddo
      enddo
    enddo
    
    ! setting dy
    do ji=1,nlins+1
      do jk=1,ncols+2
        do jl=1,nlays
          grid%z_geo_v(ji,jk,jl) = 0.5_wp*(grid%z_geo_scal(ji,jk,jl) + grid%z_geo_scal(ji+1,jk,jl))
          grid%dy(ji,jk,jl) = dy*(re + grid%z_geo_v(ji,jk,jl))/re
        enddo
      enddo
    enddo
    
    ! calculating the coordinate slopes
    call grad_hor_cov_extended(grid%z_geo_scal, grid%slope_x, grid%slope_y, grid)
    
    ! setting the z coordinates of the vertical vector points
    do ji=1,nlins+2
      do jk=1,ncols+2
        do jl=1,nlays
          if (jl == 1) then
            grid%z_geo_w(ji,jk,jl) = toa
          else
            grid%z_geo_w(ji,jk,jl) = 0.5_wp*(grid%z_geo_scal(ji,jk,jl-1) + grid%z_geo_scal(ji,jk,jl))
          endif
        enddo
      enddo
    enddo
    
    ! setting the vertical distances between the scalar data points
    do ji=1,nlins+2
      do jk=1,ncols+2
        do jl=1,nlays+1
          if (jl == 1) then
            grid%dz(ji,jk,jl) = 2._wp*(toa - grid%z_geo_scal(ji,jk,jl))
          elseif (jl == nlays+1) then
            grid%dz(ji,jk,jl) = 2._wp*(grid%z_geo_scal(ji,jk,jl-1) - grid%z_geo_w(ji,jk,jl))
          else
            grid%dz(ji,jk,jl) = grid%z_geo_scal(ji,jk,jl-1) - grid%z_geo_scal(ji,jk,jl)
          endif
        enddo
      enddo
    enddo
    
    ! setting the horizontal areas at the surface
    do ji=1,nlins
      do jk=1,ncols
        grid%area_z(ji,jk,nlays+1) = patch_area(grid%lat_scalar(ji+1),dlon,dlat)*(re + grid%z_geo_w(ji+1,jk+1,nlays+1))**2 &
          /re**2
      enddo
    enddo

    ! setting the horizontal areas at the higher points (above the surface)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          grid%area_z(ji,jk,jl) = grid%area_z(ji,jk,nlays+1)*(re + grid%z_geo_w(ji+1,jk+1,jl))**2 &
          /(re + grid%z_geo_w(ji+1,jk+1,nlays+1))**2
        enddo
      enddo
    enddo
    
    ! setting the vertical areas in x-direction
    do ji=1,nlins
      do jk=1,ncols+1
        do jl=1,nlays
          lower_z = 0.5_wp*(grid%z_geo_w(ji+1,jk,jl+1) + grid%z_geo_w(ji+1,jk+1,jl+1))
          upper_z = 0.5_wp*(grid%z_geo_w(ji+1,jk,jl  ) + grid%z_geo_w(ji+1,jk+1,jl  ))
          lower_length = dy*(re+lower_z)/re
          grid%area_x(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    
    ! setting the vertical areas in y-direction
    do ji=1,nlins+1
      do jk=1,ncols
        do jl=1,nlays
          lower_z = 0.5_wp*(grid%z_geo_w(ji,jk+1,jl+1) + grid%z_geo_w(ji+1,jk+1,jl+1))
          upper_z = 0.5_wp*(grid%z_geo_w(ji,jk+1,jl  ) + grid%z_geo_w(ji+1,jk+1,jl  ))
          lower_length = dx*cos(0.5_wp*(grid%lat_scalar(ji)+grid%lat_scalar(ji+1)))*(re+lower_z)/re
          grid%area_y(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    
    ! setting the horizontal dual areas
    do ji=1,nlins+1
      do jk=1,ncols+1
        do jl=1,nlays
          grid%z_geo_area_dual_z(ji,jk,jl) = 0.25_wp*(grid%z_geo_scal(ji,jk,jl)+grid%z_geo_scal(ji+1,jk,jl) &
          +grid%z_geo_scal(ji+1,jk+1,jl)+grid%z_geo_scal(ji,jk+1,jl))
          grid%area_dual_z(ji,jk,jl) = patch_area(0.5_wp*(grid%lat_scalar(ji) + grid%lat_scalar(ji+1)),dlon,dlat) &
          *(re + grid%z_geo_area_dual_z(ji,jk,jl))**2/re**2
        enddo
      enddo
    enddo
    
    ! setting the vertical dual areas in x-direction
    do ji=1,nlins+1
      do jk=1,ncols
        do jl=1,nlays+1
          if (jl==nlays+1) then
            lower_z = 0.5_wp*(grid%z_geo_w(ji,jk+1,jl) + grid%z_geo_w(ji+1,jk+1,jl))
            lower_length = grid%dy(ji,jk+1,jl-1)*(re + lower_z)/ &
            (re + 0.5_wp*(grid%z_geo_scal(ji,jk+1,jl-1) + grid%z_geo_scal(ji+1,jk+1,jl-1)))
          else
            lower_z = 0.5_wp*(grid%z_geo_scal(ji,jk+1,jl) + grid%z_geo_scal(ji+1,jk+1,jl))
            lower_length = grid%dy(ji,jk+1,jl)
          endif
          if (jl==1) then
            upper_z = 0.5_wp*(grid%z_geo_w(ji,jk+1,jl) + grid%z_geo_w(ji+1,jk+1,jl))
          else
            upper_z = 0.5_wp*(grid%z_geo_scal(ji,jk+1,jl-1) + grid%z_geo_scal(ji+1,jk+1,jl-1))
          endif
          grid%area_dual_x(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo
    
    ! setting the vertical dual areas in y-direction
    do ji=1,nlins
      do jk=1,ncols+1
        do jl=1,nlays+1
          if (jl==nlays+1) then
            lower_z = 0.5_wp*(grid%z_geo_w(ji+1,jk,jl) + grid%z_geo_w(ji+1,jk+1,jl))
            lower_length = grid%dx(ji+1,jk,jl-1)*(re + lower_z)/ &
            (re + 0.5_wp*(grid%z_geo_scal(ji+1,jk,jl-1) + grid%z_geo_scal(ji+1,jk+1,jl-1)))
          else
            lower_z = 0.5_wp*(grid%z_geo_scal(ji+1,jk,jl) + grid%z_geo_scal(ji+1,jk+1,jl))
            lower_length = grid%dx(ji+1,jk,jl)
          endif
          if (jl==1) then
            upper_z = 0.5_wp*(grid%z_geo_w(ji+1,jk,jl) + grid%z_geo_w(ji+1,jk+1,jl))
          else
            upper_z = 0.5_wp*(grid%z_geo_scal(ji+1,jk,jl-1) + grid%z_geo_scal(ji+1,jk+1,jl-1))
          endif
          grid%area_dual_y(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
        enddo
      enddo
    enddo

    ! setting the volume of the grid boxes
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          grid%volume(ji,jk,jl) = 1._wp/3._wp*((re + grid%z_geo_w(ji+1,jk+1,jl))**3 - (re + grid%z_geo_w(ji+1,jk+1,jl+1))**3) &
          /(re + grid%z_geo_w(ji+1,jk+1,jl+1))**2*grid%area_z(ji,jk,jl+1)
        enddo
      enddo
    enddo
    
    ! setting the inner product weights
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          grid%inner_product_weights(ji,jk,jl,1) = grid%area_x(ji  ,jk+1,jl  )*grid%dx(ji+1,jk+1,jl  )/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,2) = grid%area_y(ji+1,jk  ,jl  )*grid%dy(ji+1,jk+1,jl  )/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,3) = grid%area_x(ji  ,jk  ,jl  )*grid%dx(ji+1,jk  ,jl  )/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,4) = grid%area_y(ji  ,jk  ,jl  )*grid%dy(ji  ,jk+1,jl  )/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,5) = grid%area_z(ji  ,jk  ,jl  )*grid%dz(ji+1,jk+1,jl  )/(2._wp*grid%volume(ji,jk,jl))
          grid%inner_product_weights(ji,jk,jl,6) = grid%area_z(ji  ,jk  ,jl+1)*grid%dz(ji+1,jk+1,jl+1)/(2._wp*grid%volume(ji,jk,jl))
        enddo
      enddo
    enddo
    
    ! setting the TRSK weights
    ! u
    do ji=1,nlins
      do jk=1,ncols-1
        base_area = patch_area(grid%lat_scalar(ji+1),dlon,dlat)
        grid%trsk_weights_u(ji,jk,1) = (0.5_wp - patch_area(grid%lat_scalar(ji+1)+0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area) &
        *dx*cos(0.5_wp*(grid%lat_scalar(ji+1)+grid%lat_scalar(ji+2)))/(dx*cos(grid%lat_scalar(ji+1)))
        grid%trsk_weights_u(ji,jk,2) = -(0.5_wp - patch_area(grid%lat_scalar(ji+1)+0.25_wp*dlat,dlon,0.5_wp*dlat)/base_area) &
        *dy/(dx*cos(grid%lat_scalar(ji+1)))
        grid%trsk_weights_u(ji,jk,3) = -(0.5_wp - (patch_area(grid%lat_scalar(ji+1)+0.25_wp*dlat,dlon,0.5_wp*dlat) &
        +patch_area(grid%lat_scalar(ji+1)-0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat))/base_area) &
        *dx*cos(0.5_wp*(grid%lat_scalar(ji  )+grid%lat_scalar(ji+1)))/(dx*cos(grid%lat_scalar(ji+1)))
        grid%trsk_weights_u(ji,jk,4) = grid%trsk_weights_u(ji,jk,3)
        grid%trsk_weights_u(ji,jk,5) = -grid%trsk_weights_u(ji,jk,2)
        grid%trsk_weights_u(ji,jk,6) = grid%trsk_weights_u(ji,jk,1)
      enddo
    enddo
    ! v
    do ji=1,nlins-1
      do jk=1,ncols
        base_area = patch_area(grid%lat_scalar(ji+1),dlon,dlat)
        grid%trsk_weights_v(ji,jk,1) = 0.5_wp - patch_area(grid%lat_scalar(ji+1)+0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area
        grid%trsk_weights_v(ji,jk,2) = -grid%trsk_weights_v(ji,jk,1)
        base_area = patch_area(grid%lat_scalar(ji+2),dlon,dlat)
        grid%trsk_weights_v(ji,jk,3) = -(0.5_wp - patch_area(grid%lat_scalar(ji+2)-0.25_wp*dlat,0.5_wp*dlon,0.5_wp*dlat)/base_area)
        grid%trsk_weights_v(ji,jk,4) = -grid%trsk_weights_v(ji,jk,3)
      enddo
    enddo
  
  end subroutine grid_setup
  
  subroutine bg_setup(grid)
  
    ! This subroutine sets up the background state.
    type(t_grid), intent(inout) :: grid     ! the model grid
    
    ! local variables
    integer                     :: ji,jk,jl ! index variables
    real(wp)                    :: temperature, pressure
                                            ! temperature and pressure at the respective grid point
    real(wp)                    :: b,c      ! abbreviations needed for the hydrostatic initialization routine
    
    ! integrating the hydrostatic background state according to the given temperature profile and pressure in the lowest layer
    do ji=1,nlins+2
      do jk=1,ncols+2
        ! integrating from bottom to top
        do jl=nlays,1,-1
          temperature = bg_temp(grid%z_geo_scal(ji,jk,jl))
          ! lowest layer
          if (jl == nlays) then
            pressure    = bg_pres(grid%z_geo_scal(ji,jk,jl))
            grid%exner_bg(ji,jk,jl) = (pressure/p_0)**(gas_constant_diagnostics(1)/spec_heat_cap_diagnostics_p(1))
            grid%theta_bg(ji,jk,jl) = temperature/grid%exner_bg(ji,jk,jl)
          ! other layers
          else
            ! solving a quadratic equation for the Exner pressure
            b = -0.5_wp*grid%exner_bg(ji,jk,jl+1)/bg_temp(grid%z_geo_scal(ji,jk,jl+1)) &
            *(temperature - bg_temp(grid%z_geo_scal(ji,jk,jl+1)) + 2.0_wp/ &
            spec_heat_cap_diagnostics_p(1)*(geopot(grid%z_geo_scal(ji,jk,jl)) - geopot(grid%z_geo_scal(ji,jk,jl+1))))
            c = grid%exner_bg(ji,jk,jl+1)**2*temperature/bg_temp(grid%z_geo_scal(ji,jk,jl+1))
            grid%exner_bg(ji,jk,jl) = b+sqrt(b**2+c)
            grid%theta_bg(ji,jk,jl) = temperature/grid%exner_bg(ji,jk,jl)
          endif
        enddo
      enddo
    enddo
    
    ! calculating the gradient of the background Exner pressure (only needs to be done once)
    call grad(grid%exner_bg(2:nlins+1,2:ncols+1,:),grid%exner_bg_grad_u,grid%exner_bg_grad_v,grid%exner_bg_grad_w,grid)
  
  end subroutine bg_setup
  
  function patch_area(center_lat,dx_as_angle,dy_as_angle)
  
    ! calculates the area of a quadrilateral grid cell
  
    ! input
    real(wp) :: center_lat  ! latitude at the center of the patch
    real(wp) :: dx_as_angle ! delta x as angle
    real(wp) :: dy_as_angle ! delta y as angle
    ! output
    real(wp) :: patch_area  ! the result
  
    patch_area = re**2*dx_as_angle*(sin(center_lat + 0.5_wp*dy_as_angle) - sin(center_lat - 0.5_wp*dy_as_angle))
  
  end function patch_area
  
  function vertical_face_area(lower_z,upper_z,lower_length)
  
    ! calculates the surface of a vertical face
    ! input
    real(wp) :: lower_z            ! geometric height of the lower boundary of the face
    real(wp) :: upper_z            ! geometric height of the upper boundary of the face
    real(wp) :: lower_length       ! length of the lower boundary of the face
    ! output
    real(wp) :: vertical_face_area ! the result
    
    vertical_face_area = 0.5_wp*lower_length*(re + upper_z + re + lower_z) &
    /(re + lower_z)*(upper_z - lower_z)
  
  end function vertical_face_area
  
  function bg_temp(z_height)
  
    ! This function returns the temperature of the background state.
    
    real(wp), intent(in) :: z_height
    ! output
    real(wp)             :: bg_temp

    ! troposphere
    if (z_height < tropo_height) then  
      bg_temp = surface_temp - lapse_rate*z_height
    ! constant temperature layer
    elseif (z_height < inv_height) then
      bg_temp = surface_temp - lapse_rate*tropo_height
    ! inversion
    else
      bg_temp = surface_temp - lapse_rate*tropo_height + t_grad_inv*(z_height - inv_height)
    endif
  
  end function bg_temp

  function bg_pres(z_height)
  
    ! This function returns the pressure of the background state (only used in the lowest layer during the initialization).
    
    real(wp), intent(in) :: z_height
    ! output
    real(wp)             :: bg_pres

    if (z_height < inv_height) then  
      bg_pres = p_0_standard*(1 - lapse_rate*z_height/surface_temp)**(gravity/(gas_constant_diagnostics(1)*lapse_rate));
    elseif (z_height < tropo_height) then
      bg_pres = p_0_standard*(1 - lapse_rate*tropo_height/surface_temp)**(gravity/(gas_constant_diagnostics(1)*lapse_rate)) &
      *exp(-gravity*(z_height - tropo_height)/(gas_constant_diagnostics(1)*(surface_temp - lapse_rate*tropo_height)));
    else
      write(*,*) "Argument of bg_pres is above the inversion height. This is unrealistic in the lowest layer."
      write(*,*) "Aborting."
      call exit(1)
    endif
  
  end function bg_pres
  
  function geopot(z_height)
  
    ! This function returns the geopotential as a function of the geometrical height.
  
    real(wp), intent(in) :: z_height
    ! output
    real(wp)             :: geopot
    
    geopot = -gravity*re**2/(re+z_height)+gravity*re
  
  end function geopot
  
end module grid_generator






