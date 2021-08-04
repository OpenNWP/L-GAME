! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module multiplications

	! This module is a collection of various multiplications of vector and/or scalar fields.
	
	use definitions, only: t_grid,wp
	use run_nml,     only: nlins,ncols,nlays
	
	implicit none
	
	private
	
	public :: scalar_times_vector_scalar_h
	public :: scalar_times_vector
	
	contains

	subroutine scalar_times_vector_scalar_h(scalar_field,in_vector_x,in_vector_y, &
											result_vector_x,result_vector_y)
	
		! Multiplication of a scalar with a vector field at horizontal points.
		
		real(wp),     intent(in)    :: scalar_field(:,:,:)
		real(wp),     intent(in)    :: in_vector_x(:,:,:)
		real(wp),     intent(in)    :: in_vector_y(:,:,:)
		real(wp),     intent(inout) :: result_vector_x(:,:,:)
		real(wp),     intent(inout) :: result_vector_y(:,:,:)
	
		! local variables
		integer                     :: ji,jk ! loop indices
		
		result_vector_x(:,:,:) = 0._wp
		result_vector_y(:,:,:) = 0._wp
		
		do jk=1,ncols+1
			result_vector_x(:,jk,:) = 0.5_wp*(scalar_field(:,jk,:) + scalar_field(:,jk,:))*in_vector_x(:,jk,:)
		enddo
		
		do ji=1,nlins+1
			result_vector_y(ji,:,:) = 0.5_wp*(scalar_field(ji,:,:) + scalar_field(ji+1,:,:))*in_vector_y(ji,:,:)
		enddo
	
	end subroutine scalar_times_vector_scalar_h
	
	subroutine scalar_times_vector_scalar_v(scalar_field,in_vector_z,result_vector_z,grid)
	
		! Multiplication of a scalar with a vector field at vertical points.
		
		real(wp),     intent(in)    :: scalar_field(:,:,:)
		real(wp),     intent(in)    :: in_vector_z(:,:,:)
		real(wp),     intent(inout) :: result_vector_z(:,:,:)
		type(t_grid), intent(in)    :: grid
	
		! local variables
		integer                     :: jl ! loop index
		
		do jl=2,nlays
			result_vector_z(:,:,jl) = 0.5_wp*(scalar_field(:,:,jl-1) + scalar_field(:,:,jl))*in_vector_z(:,:,jl)
		enddo
		
		! linear extrapolation to the TOA
		result_vector_z(:,:,1) = (scalar_field(:,:,1) +    &
		(scalar_field(:,:,2) - scalar_field(:,:,3))        &
		/(grid%z_geo_scal(:,:,2) - grid%z_geo_scal(:,:,3)) &
		*(grid%z_geo_w(:,:,1) - grid%z_geo_scal(:,:,1)))   &
		*in_vector_z(:,:,1)
		
		! linear extrapolation to the surface
		result_vector_z(:,:,nlays+1) = (scalar_field(:,:,nlays) +    &
		(scalar_field(:,:,nlays-1) - scalar_field(:,:,nlays))        &
		/(grid%z_geo_scal(:,:,nlays-1) - grid%z_geo_scal(:,:,nlays)) &
		*(grid%z_geo_w(:,:,nlays+1) - grid%z_geo_scal(:,:,nlays)))   &
		*in_vector_z(:,:,nlays+1)
	
	end subroutine scalar_times_vector_scalar_v
	
	subroutine scalar_times_vector(scalar_field,in_vector_x,in_vector_y,in_vector_z, &
								   result_vector_x,result_vector_y,result_vector_z,grid)
	
		! Multiplication of a scalar with a vector field.
		
		real(wp),     intent(in)    :: scalar_field(:,:,:)
		real(wp),     intent(in)    :: in_vector_x(:,:,:)
		real(wp),     intent(in)    :: in_vector_y(:,:,:)
		real(wp),     intent(in)    :: in_vector_z(:,:,:)
		real(wp),     intent(inout) :: result_vector_x(:,:,:)
		real(wp),     intent(inout) :: result_vector_y(:,:,:)
		real(wp),     intent(inout) :: result_vector_z(:,:,:)
		type(t_grid), intent(in)    :: grid
		
		call scalar_times_vector_scalar_h(scalar_field,in_vector_x,in_vector_y,result_vector_x,result_vector_y)
		call scalar_times_vector_scalar_v(scalar_field,in_vector_z,result_vector_z,grid)
	
	end subroutine scalar_times_vector

end module multiplications









