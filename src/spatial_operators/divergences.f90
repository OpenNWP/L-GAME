! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file is a collection of gradient operators.

module divergence_operators

	use definitions, only: wp,t_grid
	use run_nml,     only: nlins,ncols,nlevs
	
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
				result_field(ji,jk,:) = 0._wp
			enddo
		enddo

	end subroutine divv_h

end module divergence_operators
