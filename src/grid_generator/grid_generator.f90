! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file contains the calculation of the grid properties as well as some type definitions.

module grid_generator

	use definitions,        only: wp,t_grid
	use run_nml,            only: nlins,ncols,nlays,dy,dx,toa,nlays_oro,sigma,re
	use gradient_operators, only: grad_hor_cov

	implicit none
	
	private
	
	public :: grid_setup
	
	contains
	
	subroutine grid_setup(grid)
	
		type(t_grid), intent(inout) :: grid
		! local variables
		real(wp) :: lat_left_lower ! latitude coordinate of lower left corner
		real(wp) :: lon_left_lower ! longitude coordinate of lower left corner
		real(wp) :: dlat           ! mesh size in y direction as angle
		real(wp) :: dlon           ! mesh size in x direction as angle
		integer  :: ji, jk, jl     ! loop indices
		real(wp) :: max_oro        ! variable for orography check
		real(wp) :: A              ! variable for calculating the vertical grid
		real(wp) :: B              ! variable for calculating the vertical grid
		real(wp) :: sigma_z        ! variable for calculating the vertical grid
		real(wp) :: z_rel          ! variable for calculating the vertical grid
		real(wp) :: z_vertical_vector_pre(nlays+1)
		                           ! variable for calculating the vertical grid
		real(wp) :: lower_z,upper_z,lower_length
		                           ! variables needed for area calculations
		
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
		
		! setting up the orography of the grid
		do ji=1,nlins+2
			do jk=1,ncols+2
				grid%z_geo_w(ji,jk,nlays+1) = 0._wp
			enddo
		enddo
	
		! calculating the vertical positions of the scalar points
		! the heights are defined according to z_k = A_k + B_k*z_surface with A_0 = toa, A_{NO_OF_LEVELS} = 0, B_0 = 0, B_{NO_OF_LEVELS} = 1
		do ji=1,nlins+2
			do jk=1,ncols+2
				! filling up z_vertical_vector_pre
				do jl=1,nlays+1
					z_rel = 1 - jl/nlays; ! z/toa
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
					grid%z_geo_scal(ji,jk,jl) = 0.5_wp*(z_vertical_vector_pre(jl) + z_vertical_vector_pre(jl + 1));
				enddo
			enddo
		enddo
		
		! setting dy
		do ji=1,nlins+1
			do jk=1,ncols+2
				do jl=1,nlays
					grid%dy(ji,jk,jl) = dy*(re + 0.5_wp*(grid%z_geo_scal(ji,jk,jl) + grid%z_geo_scal(ji+1,jk,jl)))/re
				enddo
			enddo
		enddo
		
		! setting dx
		do ji=1,nlins+2
			do jk=1,ncols+1
				do jl=1,nlays
					grid%dx(ji,jk,jl) = dx*(re + 0.5_wp*(grid%z_geo_scal(ji,jk,jl) + grid%z_geo_scal(ji,jk+1,jl)))/re
				enddo
			enddo
		enddo
		
		! calculating the coordinate slopes
		call grad_hor_cov(grid%z_geo_scal, grid%slope_x, grid%slope_y, grid)
		
		! setting the z coordinates of the vertical vector points
		do ji=1,nlins
			do jk=1,ncols
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
				grid%area_z(ji,jk,nlays+1) = patch_area(grid%lat_scalar(ji+1))*(re + grid%z_geo_w(ji+1,jk+1,nlays+1))**2 &
					/re**2
			enddo
		enddo

		! setting the horizontal areas at the higher points (above the surface)
		do ji=1,nlins
			do jk=1,ncols
				do jl=1,nlays+1
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
					lower_length = dx*cos(0.5_wp*(grid%lat_scalar(jk)+grid%lat_scalar(jk+1)))*(re+lower_z)/re
					grid%area_y(ji,jk,jl) = vertical_face_area(lower_z,upper_z,lower_length)
				enddo
			enddo
		enddo
		
		! setting the horizontal dual areas
		do ji=1,nlins+1
			do jk=1,ncols+1
				do jl=1,nlays
					grid%area_dual_z(ji,jk,jl) = patch_area(0.5_wp*(grid%lat_scalar(ji) + grid%lat_scalar(ji+1))) &
					*(re + 0.25_wp*(grid%z_geo_scal(ji,jk,jl)+grid%z_geo_scal(ji+1,jk,jl) &
					+grid%z_geo_scal(ji+1,jk+1,jl)+grid%z_geo_scal(ji,jk+1,jl)))**2/re**2
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
	
	end subroutine grid_setup
	
	function patch_area(center_lat)
	
		! calculates the surface of a quadrilateral grid cell
	
		! input
		! latitude at the center of the patch
		real(wp) :: center_lat
		! output
		real(wp) :: patch_area  ! the result
		! local variables
		real(wp) :: dy_as_angle ! mesh size in y-direction as angle
		real(wp) :: dx_as_angle ! mesh size in x-direction as angle
	
		dy_as_angle = dy/re
		dx_as_angle = dx/re
	
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
	
end module grid_generator






