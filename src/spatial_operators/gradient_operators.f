! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file is a collection of gradient operators.

module gradient_operators

	subroutine grad_hor(scalar_field, result_field_x, result_field_y, grid)

		! This subroutine computes the gradient of a scalar field.
		real(wp), intent(in)    :: scalar_field(:,:,:)      ! scalar field of which to calculate the gradient
		real(wp), intent(inout) :: result_field_x(:,:,:)    ! x component of the resulting vector field
		real(wp), intent(inout) :: result_field_y(:,:,:)    ! y component of the resulting vector field
		real(wp), intent(in)    :: grid                     ! the grid properties

		! calculating the x component of the gradient
		do ji = 1,ncols-1
			do jk = 1,nlins-1
				result_field_x(ji,jk:) = (scalar_field(ji+1,jk,:) - scalar_field(ji,jk,:))/dx(ji,:)
			enddo
		enddo

		! calculating the y component of the gradient
		do ji = 1,ncols-1
			do jk = 1,nlins-1
				result_field_y(ji,jk:) = (scalar_field(ji,jk,:) - scalar_field(ji,jk+1,:))/dy(jk,:)
			enddo
		enddo

	end subroutine grad_hor

end module gradient_operators
