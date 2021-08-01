! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module io

	implicit none
	
	private
	
	public :: wp

	! setting the floating point precision
	! single precision
	integer, parameter :: ps =  6
	integer, parameter :: rs = 37
	
	! double precision
	integer, parameter :: pd = 12
	integer, parameter :: rd = 37
	
	integer, parameter :: sp = SELECTED_REAL_KIND(ps,rs) ! single precission
	integer, parameter :: dp = SELECTED_REAL_KIND(pd,rd) ! double precission
	
	integer, parameter :: wp = dp                        ! working precission
	
	contains
	
	subroutine read_init()
		! reads the initial state of the model calculation
		
	end subroutine read_init
	
	subroutine write_output()
		! reads out the state of the model atmosphere
		! at a single timestep to a netcdf file
		
	end subroutine write_output
	
	subroutine bc()
		! sets the boundary conditions
		
	end subroutine bc

end module io
