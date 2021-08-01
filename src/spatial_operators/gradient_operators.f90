! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file is a collection of gradient operators.

module gradient_operators

	use definitions, only: t_grid,wp
	use run_nml,     only: nlins,ncols,nlevs
		
	implicit none
	
	private
	
	contains
	
	subroutine grad_hor_cov(scalar_field, result_field_x, result_field_y, grid)

		! This subroutine computes the gradient of a scalar field.
		real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
		real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
		real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
		type(t_grid), intent(in)    :: grid                  ! the grid properties
		! local variables
		integer                     :: ji,jk                  ! loop variables

		! calculating the x component of the gradient
		do ji = 1,nlins
			do jk = 1,ncols-1
				result_field_x(ji,jk+1,:) = (scalar_field(ji,jk+1,:) - scalar_field(ji,jk,:))/grid%dx(ji,jk,:)
			enddo
		enddo

		! calculating the y component of the gradient
		do ji = 1,nlins-1
			do jk = 1,ncols
				result_field_y(ji+1,jk,:) = (scalar_field(ji+1,jk,:) - scalar_field(ji,jk,:))/grid%dy(ji,jk,:)
			enddo
		enddo

	end subroutine grad_hor_cov

end module gradient_operators
