! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file is a collection of gradient operators.

module divergence_operators

	use io, only: wp
	use grid_generator, only: t_grid, t_vector_h
	
	implicit none
	
	private
	
	contains

	subroutine divv_h(vector_field, result_field, grid)

		! This subroutine computes the gradient of a scalar field.
		type(t_vector_h), intent(in) :: vector_field(:,:,:) ! scalar field of which to calculate the gradient
		real(wp), intent(inout)      :: result_field  (:,:,:) ! y component of the resulting vector field
		type(t_grid), intent(in)     :: grid                  ! the grid properties

		! performing the actual calculation
		do ji = 1,ncols
			do jk = 1,nlins
				result_field(ji,jk,:) = (vector_field%x(ji,jk,:) - vector_field%x(ji,jk+1,:))/dy(jk,:)
			enddo
		enddo

	end subroutine divv_h

end module divergence_operators
