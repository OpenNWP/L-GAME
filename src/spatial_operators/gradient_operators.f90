! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file is a collection of gradient operators.

module gradient_operators

	use grid_generator, only: t_vector, t_grid
		
	implicit none
	
	private
	
	contains
	
	subroutine grad_hor(scalar_field, result_field, grid)

		! This subroutine computes the gradient of a scalar field.
		type(t_scalar), intent(in)      :: scalar_field(:,:,:) ! scalar field of which to calculate the gradient
		type(t_vector_h), intent(inout) :: result_field        ! resulting vector field
		type(t_grid),   intent(in)      :: grid                ! the grid properties

		! calculating the x component of the gradient
		do ji = 1,nlins
			do jk = 1,ncols-1
				result_field%x(ji,jk+1,:) = (scalar_field(ji,jk+1,:) - scalar_field(ji,jk,:))/grid=>dx(ji,:)
			enddo
		enddo

		! calculating the y component of the gradient
		do ji = 1,nlins-1
			do jk = 1,ncols
				result_field%y(ji+1,jk:) = (scalar_field(ji,jk,:) - scalar_field(ji+1,jk,:))/grid=>dy(:)
			enddo
		enddo

	end subroutine grad_hor

end module gradient_operators
