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
		real(wp)                 :: c_vector(nlays-2) ! needed for the vertical solver
		real(wp)                 :: d_vector(nlays-1) ! needed for the vertical solver
		real(wp)                 :: e_vector(nlays-2) ! needed for the vertical solver
		real(wp)                 :: r_vector(nlays-1) ! needed for the vertical solver
		real(wp)                 :: solution(nlays-1) ! covariant mass flux density at the interfaces (solution)
		integer                  :: ji,jk             ! loop variables

		do ji=1,nlins
			do jk=1,ncols
			
				
		
				call thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution,nlays-1)
				
			enddo
		enddo

	end subroutine three_band_solver_ver
	
	subroutine thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,solution_length)

		! This subroutine solves a system of linear equations with a three-band matrix.
		
		real(wp), intent(in)    :: c_vector(:)
		real(wp), intent(in)    :: d_vector(:)
		real(wp), intent(in)    :: e_vector(:)
		real(wp), intent(in)    :: r_vector(:)
		real(wp), intent(inout) :: solution_vector(:) ! vector containing the solution
		integer,  intent(in)    :: solution_length    ! length of the solution vector
		
		! local variables
		real(wp) :: e_prime_vector(solution_length-1) ! help vector for solving the matrix equation
		real(wp) :: r_prime_vector(solution_length)   ! help vector for solving the matrix equation
		integer  :: jl                                ! loop index
		
		! downward sweep (matrix)
		if (d_vector(1) /= 0._wp) then
			e_prime_vector(1) = e_vector(1)/d_vector(1)
		else
			e_prime_vector(1) = 0._wp
		endif
		do jl=2,solution_length-1
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
		do jl=2,solution_length
			if (d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1) /= 0) then
				r_prime_vector(jl) = (r_vector(jl) - r_prime_vector(jl-1)*c_vector(jl-1))/(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
			else
				r_prime_vector(jl) = 0._wp
			endif
		enddo
		
		! upward sweep (final solution)
		solution_vector(solution_length) = r_prime_vector(solution_length)
		do jl=solution_length-1,1,-1
			solution_vector(jl) = r_prime_vector(jl) - e_prime_vector(jl)*solution_vector(jl+1)
		enddo
	
	end subroutine thomas_algorithm

end module vertical_slice_solvers











