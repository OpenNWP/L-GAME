! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module averaging

	use definitions, only: t_grid,wp
	use run_nml,     only: nlays, nlays_oro

	! This module contains averaging operators.
	
	implicit none
	
	private
	
	public :: vertical_contravariant_corr
	
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
				0.5_wp*vector_field_x(ji  ,jk+1,jl  )*grid%slope_x(ji  ,jk+1,jl  )*grid%inner_product_weights(ji,jk,jl  ,1) - &
				0.5_wp*vector_field_y(ji  ,jk  ,jl  )*grid%slope_y(ji  ,jk  ,jl  )*grid%inner_product_weights(ji,jk,jl  ,2) - &
				0.5_wp*vector_field_x(ji  ,jk  ,jl  )*grid%slope_x(ji  ,jk  ,jl  )*grid%inner_product_weights(ji,jk,jl  ,3) - &
				0.5_wp*vector_field_y(ji+1,jk  ,jl  )*grid%slope_y(ji+1,jk  ,jl  )*grid%inner_product_weights(ji,jk,jl  ,4)
			! level at the surface
			elseif (jl == nlays + 1) then
				vertical_contravariant_corr = vertical_contravariant_corr - &
				vector_field_x(ji  ,jk+1,jl-1)*grid%slope_x(ji  ,jk+1,jl-1)*grid%inner_product_weights(ji,jk,jl-1,1) - &
				vector_field_y(ji  ,jk  ,jl-1)*grid%slope_y(ji  ,jk  ,jl-1)*grid%inner_product_weights(ji,jk,jl-1,2) - &
				vector_field_x(ji  ,jk  ,jl-1)*grid%slope_x(ji  ,jk  ,jl-1)*grid%inner_product_weights(ji,jk,jl-1,3) - &
				vector_field_y(ji+1,jk  ,jl-1)*grid%slope_y(ji+1,jk  ,jl-1)*grid%inner_product_weights(ji,jk,jl-1,4)
			! levels in between
			else
				vertical_contravariant_corr = vertical_contravariant_corr - &
				0.5_wp*vector_field_x(ji  ,jk+1,jl-1)*grid%slope_x(ji  ,jk+1,jl-1)*grid%inner_product_weights(ji,jk,jl-1,1) - &
				0.5_wp*vector_field_y(ji  ,jk  ,jl-1)*grid%slope_y(ji  ,jk  ,jl-1)*grid%inner_product_weights(ji,jk,jl-1,2) - &
				0.5_wp*vector_field_x(ji  ,jk  ,jl-1)*grid%slope_x(ji  ,jk  ,jl-1)*grid%inner_product_weights(ji,jk,jl-1,3) - &
				0.5_wp*vector_field_y(ji+1,jk  ,jl-1)*grid%slope_y(ji+1,jk  ,jl-1)*grid%inner_product_weights(ji,jk,jl-1,4)
				vertical_contravariant_corr = vertical_contravariant_corr - &
				0.5_wp*vector_field_x(ji  ,jk+1,jl  )*grid%slope_x(ji  ,jk+1,jl  )*grid%inner_product_weights(ji,jk,jl  ,1) - &
				0.5_wp*vector_field_y(ji  ,jk  ,jl  )*grid%slope_y(ji  ,jk  ,jl  )*grid%inner_product_weights(ji,jk,jl  ,2) - &
				0.5_wp*vector_field_x(ji  ,jk  ,jl  )*grid%slope_x(ji  ,jk  ,jl  )*grid%inner_product_weights(ji,jk,jl  ,3) - &
				0.5_wp*vector_field_y(ji+1,jk  ,jl  )*grid%slope_y(ji+1,jk  ,jl  )*grid%inner_product_weights(ji,jk,jl  ,4)
			endif
		endif
	
	end function vertical_contravariant_corr

end module averaging







