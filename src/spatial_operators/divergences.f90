! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file is a collection of gradient operators.

module divergence_operators

	use definitions, only: wp,t_grid
	use run_nml,     only: nlins,ncols,nlays
	
	implicit none
	
	private
	
	contains

	subroutine divv_h(vector_field_x, vector_field_y, result_field, grid)

		! This subroutine computes the gradient of a scalar field.
		real(wp),     intent(in)    :: vector_field_x(:,:,:) ! x-component of horizontal vector field of which to calculate the divergence
		real(wp),     intent(in)    :: vector_field_y(:,:,:) ! y-component of horizontal vector field of which to calculate the divergence
		real(wp),     intent(inout) :: result_field(:,:,:)   ! resulting scalar field
		type(t_grid), intent(in)    :: grid                  ! the grid properties
		! local variables
		integer                     :: ji,jk                  ! loop variables

		! performing the actual calculation
		do ji = 1,nlins
			do jk = 1,ncols
				result_field(ji,jk,:) = (&
				vector_field_x(ji  ,jk+1,:)*grid%area_x(ji  ,jk+1,:) - &
				vector_field_y(ji  ,jk  ,:)*grid%area_y(ji  ,jk  ,:) - &
				vector_field_x(ji  ,jk  ,:)*grid%area_x(ji  ,jk  ,:) + &
				vector_field_y(ji+1,jk  ,:)*grid%area_y(ji+1,jk  ,:)) &
				/grid%volume(ji,jk,:)
			enddo
		enddo

	end subroutine divv_h

end module divergence_operators
