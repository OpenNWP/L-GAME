! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module vertical_slice_solvers

	! This module contains the implicit vertical routines (implicit part of the HEVI scheme).

	use run_nml,     only: nlins,ncols,wp,nlays
	use definitions, only: t_grid

	implicit none
	
	private
	
	public :: three_band_solver_ver
	
	contains
	
	subroutine three_band_solver_ver(grid)

		type(t_grid), intent(in) :: grid  ! model grid

		! local variables
		integer                  :: ji,jk ! loop variables

		do ji=1,nlins
			do jk=1,ncols
		
				
				
			enddo
		enddo

	end subroutine three_band_solver_ver
	
	subroutine thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector)

		! This subroutine solves a system of linear equations with a three-band matrix.
		
		real(wp), intent(in)    :: c_vector(:)
		real(wp), intent(in)    :: d_vector(:)
		real(wp), intent(in)    :: e_vector(:)
		real(wp), intent(in)    :: r_vector(:)
		real(wp), intent(inout) :: solution_vector(:)
		
		! local variables
		real(wp) :: e_prime_vector(nlays-1) ! help vector for solving the matrix equation
		real(wp) :: r_prime_vector(nlays)   ! help vector for solving the matrix equation
		integer  :: jl                      ! loop index
		
		! downward sweep (matrix)
		if (d_vector(1) /= 0._wp) then
			e_prime_vector(1) = e_vector(1)/d_vector(1)
		else
			e_prime_vector(1) = 0._wp
		endif
		do jl=2,nlays-1
			if (d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1) /= 0) then
				e_prime_vector(jl) = e_vector(jl)/(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
			else
				e_prime_vector(jl) = 0._wp
			endif
		enddo
		! downward sweep (right-hand side)
		if (d_vector(1) /= 0) then
			r_prime_vector(1) = r_vector(1)/d_vector(1)
		else
			r_prime_vector(1) = 0._wp
		endif
		do jl=2,nlays
			if (d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1) /= 0) then
				r_prime_vector(jl) = (r_vector(jl) - r_prime_vector(jl-1)*c_vector(jl-1))/(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
			else
				r_prime_vector(jl) = 0._wp
			endif
		enddo
		
		! upward sweep (final solution)
		solution_vector(nlays) = r_prime_vector(nlays)
		do jl=nlays-1,1,-1
			solution_vector(jl) = r_prime_vector(jl) - e_prime_vector(jl)*solution_vector(jl+1)
		enddo
	
	end subroutine thomas_algorithm

end module vertical_slice_solvers











