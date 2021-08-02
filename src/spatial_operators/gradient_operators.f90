! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file is a collection of gradient operators.

module gradient_operators

	use definitions, only: t_grid,wp
	use run_nml,     only: nlins,ncols,nlays,toa
		
	implicit none
	
	private
	
	public :: grad_hor_cov
	public :: grad
	
	contains
	
	subroutine grad_hor_cov(scalar_field,result_field_x,result_field_y,grid)

		! This subroutine computes the horizontal covariant gradient of a scalar field.
		real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
		real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
		real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
		type(t_grid), intent(in)    :: grid                  ! the grid properties
		! local variables
		integer                     :: ji,jk                  ! loop variables

		! calculating the x component of the gradient
		do ji=1,nlins
			do jk=1,ncols-1
				result_field_x(ji,jk+1,:) = (scalar_field(ji,jk+1,:) - scalar_field(ji,jk,:))/grid%dx(ji,jk,:)
			enddo
		enddo

		! calculating the y component of the gradient
		do ji=1,nlins-1
			do jk=1,ncols
				result_field_y(ji+1,jk,:) = (scalar_field(ji+1,jk,:) - scalar_field(ji,jk,:))/grid%dy(ji,jk,:)
			enddo
		enddo

	end subroutine grad_hor_cov
	
	subroutine grad_vert_cov(scalar_field,result_field,grid)

		! This subroutine computes the vertical covariant gradient of a scalar field.
		real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
		real(wp),     intent(inout) :: result_field(:,:,:)   ! z-component of resulting vector field
		type(t_grid), intent(in)    :: grid                  ! the grid properties
		! local variables
		integer                     :: ji,jk,jl              ! loop variables

		! calculating the vertical gradient in the inner levels
		do ji=1,nlins
			do jk=1,ncols-1
				do jl=2,nlays
					result_field(ji,jk,jl) = (scalar_field(ji,jk,jl-1) - scalar_field(ji,jk,jl))/grid%dz(ji,jk,jl)
				enddo
			enddo
		enddo

		! linear extrapolation to the TOA
		do ji=1,nlins
			do jk=1,ncols
				result_field(ji,jk,1) = result_field(ji,jk,2) + &
				(result_field(ji,jk,2) - result_field(ji,jk,3))/(grid%z_geo_w(ji,jk,2) - grid%z_geo_w(ji,jk,3)) &
				*(toa - grid%z_geo_w(ji,jk,2))
			enddo
		enddo
		
		! linear extrapolation to the surface
		do ji=1,nlins
			do jk=1,ncols
				result_field(ji,jk,nlays+1) = result_field(ji,jk,nlays) + &
				(result_field(ji,jk,nlays-1) - result_field(ji,jk,nlays))/(grid%z_geo_w(ji,jk,nlays-1) - grid%z_geo_w(ji,jk,nlays)) &
				*(grid%z_geo_w(ji,jk,nlays+1) - grid%z_geo_w(ji,jk,nlays))
			enddo
		enddo

	end subroutine grad_vert_cov
	
	subroutine grad_cov(scalar_field,result_field_x,result_field_y,result_field_z,grid)
	
		! This subroutine computes the covariant gradient of a scalar field.
		real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
		real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
		real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
		real(wp),     intent(inout) :: result_field_z(:,:,:) ! z-component of resulting vector field
		type(t_grid), intent(in)    :: grid                  ! the grid properties
		
		call grad_hor_cov(scalar_field,result_field_x,result_field_y,grid)
		call grad_vert_cov(scalar_field,result_field_z,grid)
	
	end subroutine grad_cov
	
	subroutine grad(scalar_field,result_field_x,result_field_y,result_field_z,grid)
	
		! This subroutine computes the covariant gradient of a scalar field.
		real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
		real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
		real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
		real(wp),     intent(inout) :: result_field_z(:,:,:) ! z-component of resulting vector field
		type(t_grid), intent(in)    :: grid                  ! the grid properties
		
		call grad_cov(scalar_field,result_field_x,result_field_y,result_field_z,grid)
	
	end subroutine grad

end module gradient_operators








