! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module averaging

	use definitions, only: t_grid,wp
	use run_nml,     only: nlays,nlays_oro,nlins,ncols

	! This module contains averaging operators.
	
	implicit none
	
	private
	
	public :: vertical_contravariant_corr
	public :: hor_cov_to_con
	
	contains
	
	function vertical_contravariant_corr(vector_field_x,vector_field_y,ji,jk,jl,grid)
	
		! calculates (the vertical contravariant component - the vertical covariant component)
		! of a vector field out of the horizontal contravariant components
		
		real(wp),     intent(in) :: vector_field_x(:,:,:) ! x-component of vector field to work with
		real(wp),     intent(in) :: vector_field_y(:,:,:) ! y-component of vector field to work with
		integer,      intent(in) :: ji,jk,jl              ! spatial indices
		type(t_grid), intent(in) :: grid                  ! model grid
		
		real(wp) :: vertical_contravariant_corr
		vertical_contravariant_corr = 0._wp
		
		if (jl >= nlays - nlays_oro + 1) then
		    ! highest level following orography
			if (jl == nlays - nlays_oro + 1) then
				vertical_contravariant_corr = vertical_contravariant_corr - &
				0.5_wp*vector_field_x(ji+1,jk+1,jl  )*grid%slope_x(ji+1,jk+1,jl  )*grid%inner_product_weights(ji,jk,jl  ,1) - &
				0.5_wp*vector_field_y(ji+1,jk+1,jl  )*grid%slope_y(ji+1,jk+1,jl  )*grid%inner_product_weights(ji,jk,jl  ,2) - &
				0.5_wp*vector_field_x(ji+1,jk  ,jl  )*grid%slope_x(ji+1,jk  ,jl  )*grid%inner_product_weights(ji,jk,jl  ,3) - &
				0.5_wp*vector_field_y(ji+1,jk+1,jl  )*grid%slope_y(ji+1,jk+1,jl  )*grid%inner_product_weights(ji,jk,jl  ,4)
			! level at the surface
			elseif (jl == nlays + 1) then
				vertical_contravariant_corr = vertical_contravariant_corr - &
				vector_field_x(ji+1,jk+1,jl-1)*grid%slope_x(ji+1,jk+1,jl-1)*grid%inner_product_weights(ji,jk,jl-1,1) - &
				vector_field_y(ji+1,jk+1,jl-1)*grid%slope_y(ji+1,jk+1,jl-1)*grid%inner_product_weights(ji,jk,jl-1,2) - &
				vector_field_x(ji+1,jk  ,jl-1)*grid%slope_x(ji+1,jk  ,jl-1)*grid%inner_product_weights(ji,jk,jl-1,3) - &
				vector_field_y(ji+1,jk+1,jl-1)*grid%slope_y(ji+1,jk+1,jl-1)*grid%inner_product_weights(ji,jk,jl-1,4)
			! levels in between
			else
				vertical_contravariant_corr = vertical_contravariant_corr - &
				0.5_wp*vector_field_x(ji+1,jk+1,jl-1)*grid%slope_x(ji+1,jk+1,jl-1)*grid%inner_product_weights(ji,jk,jl-1,1) - &
				0.5_wp*vector_field_y(ji+1,jk+1,jl-1)*grid%slope_y(ji+1,jk+1,jl-1)*grid%inner_product_weights(ji,jk,jl-1,2) - &
				0.5_wp*vector_field_x(ji+1,jk  ,jl-1)*grid%slope_x(ji+1,jk  ,jl-1)*grid%inner_product_weights(ji,jk,jl-1,3) - &
				0.5_wp*vector_field_y(ji+1,jk+1,jl-1)*grid%slope_y(ji+1,jk+1,jl-1)*grid%inner_product_weights(ji,jk,jl-1,4)
				vertical_contravariant_corr = vertical_contravariant_corr - &
				0.5_wp*vector_field_x(ji+1,jk+1,jl  )*grid%slope_x(ji+1,jk+1,jl  )*grid%inner_product_weights(ji,jk,jl  ,1) - &
				0.5_wp*vector_field_y(ji+1,jk+1,jl  )*grid%slope_y(ji+1,jk+1,jl  )*grid%inner_product_weights(ji,jk,jl  ,2) - &
				0.5_wp*vector_field_x(ji+1,jk  ,jl  )*grid%slope_x(ji+1,jk  ,jl  )*grid%inner_product_weights(ji,jk,jl  ,3) - &
				0.5_wp*vector_field_y(ji+1,jk+1,jl  )*grid%slope_y(ji+1,jk+1,jl  )*grid%inner_product_weights(ji,jk,jl  ,4)
			endif
		endif
	
	end function vertical_contravariant_corr
	
	subroutine hor_cov_to_con(result_field_x,result_field_y,result_field_z,grid)
	
		! This subroutine computes the terrain correction of the gradient.
		real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
		real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
		real(wp),     intent(in)    :: result_field_z(:,:,:) ! z-component of resulting vector field
		type(t_grid), intent(in)    :: grid                  ! the grid properties
	
		! local variables
		integer                     :: ji,jk,jl              ! loop indices
		
		! correction to the x-component
		do ji=1,nlins
			do jk=1,ncols-1
				do jl=1,nlays
					result_field_x(ji,jk,jl) = result_field_x(ji,jk,jl) - &
					grid%slope_x(ji+1,jk+1,jl)*remap_ver2hor_x(result_field_z,grid,ji,jk,jl)
				enddo
			enddo
		enddo
		
		! correction to the y-component
		do ji=1,nlins-1
			do jk=1,ncols
				do jl=1,nlays
					result_field_y(ji,jk,jl) = result_field_y(ji,jk,jl) - &
					grid%slope_y(ji+1,jk+1,jl)*remap_ver2hor_y(result_field_z,grid,ji,jk,jl)
				enddo
			enddo
		enddo
	
	end subroutine hor_cov_to_con
	
	function remap_ver2hor_x(vertical_cov,grid,ji,jk,jl)
	
		! This function remaps a vertical covariant component of
		! a vector field to a position of a vector component in x-direction.
	
		real(wp),     intent(in)    :: vertical_cov(:,:,:)   ! z-component of vector field to work with
		type(t_grid), intent(in)    :: grid                  ! the grid properties
		integer,      intent(in)    :: ji,jk,jl              ! positional indices
	
		real(wp)                    :: remap_ver2hor_x       ! the result
		
	    remap_ver2hor_x = grid%inner_product_weights(ji,jk  ,jl,5)*vertical_cov(ji,jk  ,jl)
	    remap_ver2hor_x = grid%inner_product_weights(ji,jk+1,jl,5)*vertical_cov(ji,jk+1,jl)
		! layer below
		if (jl < nlays) then
			remap_ver2hor_x = grid%inner_product_weights(ji,jk  ,jl,6)*vertical_cov(ji,jk  ,jl+1)
			remap_ver2hor_x = grid%inner_product_weights(ji,jk+1,jl,6)*vertical_cov(ji,jk+1,jl+1)
		endif
		! horizontal average
    	remap_ver2hor_x = 0.5_wp*remap_ver2hor_x
	
	end function remap_ver2hor_x
	
	function remap_ver2hor_y(vertical_cov,grid,ji,jk,jl)
	
		! This function remaps a vertical covariant component of
		! a vector field to a position of a vector component in y-direction.
	
		real(wp),     intent(in)    :: vertical_cov(:,:,:)   ! z-component of vector field to work with
		type(t_grid), intent(in)    :: grid                  ! the grid properties
		integer,      intent(in)    :: ji,jk,jl              ! positional indices
	
		real(wp)                    :: remap_ver2hor_y       ! the result
		
	    remap_ver2hor_y = grid%inner_product_weights(ji  ,jk,jl,5)*vertical_cov(ji  ,jk,jl)
	    remap_ver2hor_y = grid%inner_product_weights(ji+1,jk,jl,5)*vertical_cov(ji+1,jk,jl)
		! layer below
		if (jl < nlays) then
			remap_ver2hor_y = grid%inner_product_weights(ji  ,jk,jl,6)*vertical_cov(ji  ,jk,jl+1)
			remap_ver2hor_y = grid%inner_product_weights(ji+1,jk,jl,6)*vertical_cov(ji+1,jk,jl+1)
		endif
		! horizontal average
    	remap_ver2hor_y = 0.5_wp*remap_ver2hor_y
	
	end function remap_ver2hor_y

end module averaging







