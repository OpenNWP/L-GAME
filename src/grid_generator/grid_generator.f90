! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file contains the calculation of the grid properties as well as some type definitions.

module grid_generator

	use io, only: wp

	implicit none
	
	private
	
	public :: t_scalar
	public :: t_vector
	public :: t_grid
	public :: t_state
	public :: diag
	
	type t_scalar
		
		real(wp), allocatable :: values(:,:,:)
	
	end type t_scalar
	
	type t_vector_h
	
		real(wp), allocatable :: x(:,:,:)
		real(wp), allocatable :: y(:,:,:)
	
	end type t_vector_h
	
	type t_grid
	
		t_scalar,   allocatable :: z_geo_scal
		t_scalar,   allocatable :: z_agl_scal
		t_scalar,   allocatable :: volume
		t_vector_h, allocatable :: area_h
	
	end type t_grid
	
	type t_state
	
		! type containing the state variables
		type(t_scalar),   allocatable: rho
		type(t_scalar),   allocatable: rhotheta
		type(t_scalar),   allocatable: exner
		type(t_vector_h), allocatable: wind_h
		type(t_vector_v), allocatable: wind_w
	
	end type t_state
	
	type t_diag
	
		type(t_scalar),   allocatable: theta
	
	end type t_diag
	
end module grid_generator






