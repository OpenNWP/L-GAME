! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file is a collection of gradient operators.

module divergence_operators

	use grid_generator, only: t_grid
	
	implicit none
	
	private
	
	contains

	subroutine divv_h(vector_field, result_field, grid)

		! This subroutine computes the gradient of a scalar field.
		real(wp), intent(in)     :: vector_field_x(:,:,:) ! scalar field of which to calculate the gradient
		real(wp), intent(inout)  :: vector_field_y(:,:,:) ! x component of the resulting vector field
		real(wp), intent(inout)  :: result_field  (:,:,:) ! y component of the resulting vector field
		type(t_grid), intent(in) :: grid                  ! the grid properties

		! performing the actual calculation
		do ji = 1,ncols
			do jk = 1,nlins
				result_field(ji,jk,:) = (scalar_field(ji,jk,:) - scalar_field(ji,jk+1,:))/dy(jk,:)
			enddo
		enddo

	end subroutine divv_h

end module divergence_operators
