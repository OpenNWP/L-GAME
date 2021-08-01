! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file contains the calculation of the grid properties as well as some type definitions.

module grid_generator

	use io, only: wp

	implicit none
	
	private
	
	public :: t_vector_h
	public :: t_grid
	public :: t_state
	public :: grid_setup
	
	type t_vector_h
	
		real(wp), allocatable :: x(:,:,:)
		real(wp), allocatable :: y(:,:,:)
	
	end type t_vector_h
	
	type t_grid
	
		real(wp),         allocatable :: z_geo_scal(:,:,:)
		real(wp),         allocatable :: z_agl_scal(:,:,:)
		real(wp),         allocatable :: volume
		type(t_vector_h), allocatable :: area_h
		real(wp),         allocatable :: dy(:)
		real(wp),         allocatable :: dx(:,:)
	
	end type t_grid
	
	type t_state
	
		! type containing the state variables
		real(wp),         allocatable :: rho(:,:,:)
		real(wp),         allocatable :: rhotheta(:,:,:)
		real(wp),         allocatable :: exner_bg(:,:,:)
		real(wp),         allocatable :: theta_bg(:,:,:)
		real(wp),         allocatable :: theta_pert(:,:,:)
		real(wp),         allocatable :: exner_pert(:,:,:)
		type(t_vector_h), allocatable :: wind_h
		real(wp),         allocatable :: wind_v(:,:,:)
	
	end type t_state
	
	contains
	
	subroutine grid_setup
	
		real(wp) :: re ! Earth radius
	
	end subroutine grid_setup
	
end module grid_generator






