! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file containss some definitions.

module definitions
	                          
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
	
	public :: t_grid
	public :: t_state
	public :: t_bg
	public :: t_diag
	
	type t_grid
	
		real(wp), allocatable :: lat_scalar(:)
		real(wp), allocatable :: lon_scalar(:)
		real(wp), allocatable :: z_geo_scal(:,:,:)
		real(wp), allocatable :: z_agl_scal(:,:,:)
		real(wp), allocatable :: z_geo_w(:,:,:)
		real(wp), allocatable :: volume(:,:,:)
		real(wp), allocatable :: area_x(:,:,:)
		real(wp), allocatable :: area_y(:,:,:)
		real(wp), allocatable :: area_z(:,:,:)
		real(wp), allocatable :: dy(:,:,:)
		real(wp), allocatable :: dx(:,:,:)
		real(wp), allocatable :: dz(:,:,:)
		real(wp), allocatable :: slope_x(:,:,:)
		real(wp), allocatable :: slope_y(:,:,:)
	
	end type t_grid
	
	type t_state
	
		! type containing the state variables
		real(wp), allocatable :: rho(:,:,:)
		real(wp), allocatable :: rhotheta(:,:,:)
		real(wp), allocatable :: theta_pert(:,:,:)
		real(wp), allocatable :: exner_pert(:,:,:)
		real(wp), allocatable :: wind_u(:,:,:)
		real(wp), allocatable :: wind_v(:,:,:)
		real(wp), allocatable :: wind_w(:,:,:)
	
	end type t_state
	
		
	type t_bg
	
		! background state
		real(wp), allocatable :: exner_bg(:,:,:)
		real(wp), allocatable :: theta_bg(:,:,:)
	
	end type t_bg
	
	type t_diag
	
		! type containing diagnostic quantities
		real(wp), allocatable :: e_kin(:,:,:)
	
	end type t_diag
	
end module definitions






