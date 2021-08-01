! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module io

	use definitions, only: t_state, wp

	implicit none
	
	private
	
	public :: read_init
	
	contains
	
	subroutine read_init(state)
		! reads the initial state of the model calculation
		
		type(t_state), intent(inout) :: state
		
	end subroutine read_init
	
	subroutine write_output()
		! reads out the state of the model atmosphere
		! at a single timestep to a netcdf file
		
	end subroutine write_output
	
	subroutine bc()
		! sets the boundary conditions
		
	end subroutine bc
	
	subroutine oi()
	
		! optimum interpolation
		
	
	end subroutine oi
	
	subroutine var_3d()
		
		! three-dimensional variational data assimilation
	
	end subroutine var_3d
	
	subroutine var_4d()
	
		! four-dimensional variational data assimilation
	
	end subroutine var_4d

end module io








