! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

! This file containss some definitions.

module definitions
	                          
	implicit none
	
	private
	
	public :: wp
	public :: t_grid
	public :: t_state
	public :: t_bg
	public :: t_diag
	public :: t_tend
	
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
	
	type t_grid
	
		real(wp), allocatable :: lat_scalar(:)
		real(wp), allocatable :: lon_scalar(:)
		real(wp), allocatable :: z_geo_scal(:,:,:)
		real(wp), allocatable :: dx(:,:,:)
		real(wp), allocatable :: dy(:,:,:)
		real(wp), allocatable :: dz(:,:,:)
		real(wp), allocatable :: z_geo_u(:,:,:)
		real(wp), allocatable :: z_geo_v(:,:,:)
		real(wp), allocatable :: z_geo_w(:,:,:)
		real(wp), allocatable :: volume(:,:,:)
		real(wp), allocatable :: area_x(:,:,:)
		real(wp), allocatable :: area_y(:,:,:)
		real(wp), allocatable :: area_z(:,:,:)
		real(wp), allocatable :: slope_x(:,:,:)
		real(wp), allocatable :: slope_y(:,:,:)
		real(wp), allocatable :: inner_product_weights(:,:,:,:)
		real(wp), allocatable :: area_dual_x(:,:,:)
		real(wp), allocatable :: area_dual_y(:,:,:)
		real(wp), allocatable :: area_dual_z(:,:,:)
		real(wp), allocatable :: z_geo_area_dual_z(:,:,:)
		real(wp), allocatable :: fvec_x(:,:)                    ! x-component of Coriolis vector
		real(wp), allocatable :: fvec_y(:,:)                    ! y-component of Coriolis vector
		real(wp), allocatable :: fvec_z(:,:)                    ! z-component of Coriolis vector
		real(wp), allocatable :: trsk_weights_u(:,:,:)
		real(wp), allocatable :: trsk_weights_v(:,:,:)
	
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
	
	type t_tend
	
		! type containing tendencies
		real(wp), allocatable :: rho(:,:,:)
		real(wp), allocatable :: rhotheta(:,:,:)
		real(wp), allocatable :: wind_u(:,:,:)
		real(wp), allocatable :: wind_v(:,:,:)
		real(wp), allocatable :: wind_w(:,:,:)
	
	end type t_tend
	
		
	type t_bg
	
		! background state
		real(wp), allocatable :: theta(:,:,:)              ! potential temperature
		real(wp), allocatable :: exner(:,:,:)              ! Exner pressure
	
	end type t_bg
	
	type t_diag
	
		! type containing diagnostic quantities
		real(wp), allocatable :: e_kin(:,:,:)              ! specific kinetic energy
		real(wp), allocatable :: p_grad_acc_l_u(:,:,:)     ! x-component of linear pressure gradient acceleration
		real(wp), allocatable :: p_grad_acc_l_v(:,:,:)     ! y-component of linear pressure gradient acceleration
		real(wp), allocatable :: p_grad_acc_l_w(:,:,:)     ! z-component of linear pressure gradient accelerationpgrad_acc_old
		real(wp), allocatable :: p_grad_acc_nl_u(:,:,:)    ! x-component of nonlinear pressure gradient acceleration
		real(wp), allocatable :: p_grad_acc_nl_v(:,:,:)    ! y-component of nonlinear pressure gradient acceleration
		real(wp), allocatable :: p_grad_acc_nl_w(:,:,:)    ! z-component of nonlinear pressure gradient acceleration
		real(wp), allocatable :: p_grad_acc_old_u(:,:,:)   ! x-component of pressure gradient at old time step
		real(wp), allocatable :: p_grad_acc_old_v(:,:,:)   ! y-component of pressure gradient at old time step
		real(wp), allocatable :: p_grad_acc_old_w(:,:,:)   ! z-component of pressure gradient at old time step
		real(wp), allocatable :: e_kin_grad_x(:,:,:)       ! x-gradient of specific kinetic energy
		real(wp), allocatable :: e_kin_grad_y(:,:,:)       ! y-gradient of specific kinetic energy
		real(wp), allocatable :: e_kin_grad_z(:,:,:)       ! z-gradient of specific kinetic energy
		real(wp), allocatable :: pot_vort_tend_x(:,:,:)    ! tendency due to the vorticity flux term in x-direction
		real(wp), allocatable :: pot_vort_tend_y(:,:,:)    ! tendency due to the vorticity flux term in y-direction
		real(wp), allocatable :: pot_vort_tend_z(:,:,:)    ! tendency due to the vorticity flux term in z-direction
		real(wp), allocatable :: mom_diff_tend_x(:,:,:)    ! tendency due to momentum diffusion in x-direction
		real(wp), allocatable :: mom_diff_tend_y(:,:,:)    ! tendency due to momentum diffusion in y-direction
		real(wp), allocatable :: mom_diff_tend_z(:,:,:)    ! tendency due to momentum diffusion in z-direction
		real(wp), allocatable :: scalar_placeholder(:,:,:) ! placeholder for scalar fields
		real(wp), allocatable :: u_placeholder(:,:,:)      ! placeholder for vector fields in x-direction
		real(wp), allocatable :: v_placeholder(:,:,:)      ! placeholder for vector fields in y-direction
		real(wp), allocatable :: u_10(:,:)                 ! 10 m wind in x direction
		real(wp), allocatable :: v_10(:,:)                 ! 10 m wind in y direction
		real(wp), allocatable :: mslp(:,:)                 ! mean sea level pressure
		real(wp), allocatable :: t_2(:,:)                  ! 2 m temperature
		real(wp), allocatable :: z_eta_x(:,:,:)             ! relative vorticity in x-direction
		real(wp), allocatable :: z_eta_y(:,:,:)             ! relative vorticity in y-direction
		real(wp), allocatable :: z_eta_z(:,:,:)             ! relative vorticity in z-direction
	
	end type t_diag
	
end module definitions






