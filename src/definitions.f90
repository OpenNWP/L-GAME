! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module definitions

  ! This file contains some definitions.
                            
  implicit none
  
  private
  
  public :: wp
  public :: t_grid
  public :: t_state
  public :: t_diag
  public :: t_tend
  public :: t_bc
  public :: t_irrev
  
  ! setting the floating point precision
  ! single precision
  integer, parameter :: ps =  6
  integer, parameter :: rs = 37
  
  ! double precision
  integer, parameter :: pd = 12
  integer, parameter :: rd = 37
  
  integer, parameter :: sp = SELECTED_REAL_KIND(ps,rs)        ! single precission
  integer, parameter :: dp = SELECTED_REAL_KIND(pd,rd)        ! double precission
  
  integer, parameter :: wp = dp                               ! working precission
  
  type t_grid
  
    ! type containing information on the model grid
    real(wp), allocatable :: lat_scalar(:)                    ! latitudes of the scalar gridpoints
    real(wp), allocatable :: lon_scalar(:)                    ! longitudes of the scalar gridpoints
    real(wp), allocatable :: lat_geo_scalar(:,:)              ! geographic latitudes of the scalar gridpoints
    real(wp), allocatable :: lon_geo_scalar(:,:)              ! geographic longitudes of the scalar gridpoints
    real(wp), allocatable :: lat_geo_u(:,:)                   ! geographic latitudes of the u-vector gridpoints
    real(wp), allocatable :: lon_geo_u(:,:)                   ! geographic longitudes of the u-vector gridpoints
    real(wp), allocatable :: dir_geo_u(:,:)                   ! geographic directions of the u-vectors
    real(wp), allocatable :: lat_geo_v(:,:)                   ! geographic latitudes of the v-vector gridpoints
    real(wp), allocatable :: lon_geo_v(:,:)                   ! geographic longitudes of the v-vector gridpoints
    real(wp), allocatable :: dir_geo_v(:,:)                   ! geographic directions of the v-vectors
    real(wp), allocatable :: dir_geo_u_scalar(:,:)            ! geographic directions of the u-vectors at the scalar points
    real(wp), allocatable :: z_geo_scal(:,:,:)                ! geometric heights of the scalar gridpoints
    real(wp), allocatable :: dx(:,:,:)                        ! grid point distance in x-direction
    real(wp), allocatable :: dy(:,:,:)                        ! grid point distance in y-direction
    real(wp), allocatable :: dz(:,:,:)                        ! grid point distance in z-direction
    real(wp), allocatable :: z_geo_u(:,:,:)                   ! geomtric height of the u-vectors
    real(wp), allocatable :: z_geo_v(:,:,:)                   ! geomtric height of the v-vectors
    real(wp), allocatable :: z_geo_w(:,:,:)                   ! geomtric height of the w-vectors
    real(wp), allocatable :: volume(:,:,:)                    ! volumes of the gridboxes
    real(wp), allocatable :: area_x(:,:,:)                    ! areas of the grid in x-direction
    real(wp), allocatable :: area_y(:,:,:)                    ! areas of the grid in y-direction
    real(wp), allocatable :: area_z(:,:,:)                    ! areas of the grid in z-direction
    real(wp), allocatable :: slope_x(:,:,:)                   ! coordinate slopes in x-direction
    real(wp), allocatable :: slope_y(:,:,:)                   ! coordinate slopes in y-direction
    real(wp), allocatable :: inner_product_weights(:,:,:,:)   ! weights for calculating the inner product
    real(wp), allocatable :: area_dual_x(:,:,:)               ! areas of the dual grid in x-direction
    real(wp), allocatable :: area_dual_y(:,:,:)               ! areas of the dual grid in y-direction
    real(wp), allocatable :: area_dual_z(:,:,:)               ! areas of the dual grid in z-direction
    real(wp), allocatable :: z_geo_area_dual_z(:,:,:)         ! geometric heights of the areas of the dual grid in z-direction
    real(wp), allocatable :: fvec_x(:,:)                      ! x-component of Coriolis vector
    real(wp), allocatable :: fvec_y(:,:)                      ! y-component of Coriolis vector
    real(wp), allocatable :: fvec_z(:,:)                      ! z-component of Coriolis vector
    real(wp), allocatable :: trsk_weights_u(:,:,:)            ! weights for computing the Coriolis acceleration in x-direction
    real(wp), allocatable :: trsk_weights_v(:,:,:)            ! weights for computing the Coriolis acceleration in y-direction
    real(wp), allocatable :: exner_bg_grad_u(:,:,:)           ! gradient of background exner pressure in x-direction
    real(wp), allocatable :: exner_bg_grad_v(:,:,:)           ! gradient of background exner pressure in y-direction
    real(wp), allocatable :: exner_bg_grad_w(:,:,:)           ! gradient of background exner pressure in z-direction
    real(wp), allocatable :: theta_bg(:,:,:)                  ! background potential temperature
    real(wp), allocatable :: exner_bg(:,:,:)                  ! background Exner pressure
    real(wp), allocatable :: sfc_albedo(:,:)                  ! albedo of the surface
    real(wp), allocatable :: sfc_rho_c(:,:)                   ! volumetric heat capacity of the surface
    real(wp), allocatable :: t_conduc_soil(:,:)               ! temperature conductivity of the soil
    real(wp), allocatable :: roughness_length(:,:)            ! roughness length of the surface
    integer,  allocatable :: is_land(:,:)                     ! land-sea-mask
    real(wp), allocatable :: z_soil_interface(:)              ! heights of the interfaces of the soil layers
    real(wp), allocatable :: z_soil_center(:)                 ! heights of the centers of the soil layers
    real(wp), allocatable :: t_const_soil(:,:)                ! temperature of the soil below the depth where it is constant
    real(wp)              :: z_t_const                        ! depth where the soil temperature is constant
    real(wp)              :: lat_center                       ! latitude of the center of the model domain
    real(wp)              :: lon_center                       ! longitude of the center of the model domain
    real(wp)              :: x_dir_center                     ! direction of the x-axis of the model
    real(wp)              :: mean_velocity_area               ! area needed for the turbulence parameterizations
    
  end type t_grid
  
  type t_state
  
    ! type containing the state variables
    real(wp), allocatable :: rho(:,:,:,:)                     ! mass densities
    real(wp), allocatable :: rhotheta(:,:,:)                  ! potential temperature density
    real(wp), allocatable :: theta_pert(:,:,:)                ! potential temperature peturbation
    real(wp), allocatable :: exner_pert(:,:,:)                ! Exner pressure peturbation
    real(wp), allocatable :: condensed_rho_t(:,:,:,:)         ! temperature densities of the condensates
    real(wp), allocatable :: wind_u(:,:,:)                    ! x-component of the wind
    real(wp), allocatable :: wind_v(:,:,:)                    ! y-component of the wind
    real(wp), allocatable :: wind_w(:,:,:)                    ! vertical wind
    real(wp), allocatable :: temperature_soil(:,:,:)          ! temperature of the soil
  
  end type t_state
  
  type t_tend
  
    ! type containing tendencies
    real(wp), allocatable :: rho(:,:,:,:)                     ! mass densities
    real(wp), allocatable :: rhotheta(:,:,:)                  ! potential temperature densities
    real(wp), allocatable :: wind_u(:,:,:)                    ! x-component of the wind
    real(wp), allocatable :: wind_v(:,:,:)                    ! y-component of the wind
    real(wp), allocatable :: wind_w(:,:,:)                    ! vertical wind
    real(wp), allocatable :: condensed_rho_t(:,:,:,:)         ! temperature densities of the constituents
  
  end type t_tend
  
  type t_bc
  
    ! type containing information on boundary conditions
    real(wp), allocatable :: rho_old(:,:,:,:)                 ! mass densities at the old timestep
    real(wp), allocatable :: rhotheta_old(:,:,:)              ! potential temperature densities at the old timestep
    real(wp), allocatable :: wind_u_old(:,:,:)                ! x-component of the wind at the old timestep
    real(wp), allocatable :: wind_v_old(:,:,:)                ! y-component of the wind at the old timestep
    real(wp), allocatable :: wind_w_old(:,:,:)                ! vertical wind at the old timestep
    real(wp), allocatable :: condensed_rho_t_old(:,:,:,:)     ! temperature densities of the constituents at the old timestep
    real(wp), allocatable :: rho_new(:,:,:,:)                 ! mass densities at the new timestep
    real(wp), allocatable :: rhotheta_new(:,:,:)              ! potential temperature densities at the new timestep
    real(wp), allocatable :: wind_u_new(:,:,:)                ! x-component of the wind at the new timestep
    real(wp), allocatable :: wind_v_new(:,:,:)                ! y-component of the wind at the new timestep
    real(wp), allocatable :: wind_w_new(:,:,:)                ! vertical wind at the new timestep
    real(wp), allocatable :: condensed_rho_t_new(:,:,:,:)     ! temperature densities of the constituents at the new timestep
    real(wp), allocatable :: scalar_bc_factor(:,:)            ! boundary conditions factor for scalar fields
    real(wp), allocatable :: u_bc_factor(:,:)                 ! boundary conditions factor for u-vector fields
    real(wp), allocatable :: v_bc_factor(:,:)                 ! boundary conditions factor for v-vector fields
    
  end type t_bc
  
  type t_diag
  
    ! type containing diagnostic quantities
    real(wp), allocatable :: v_squared(:,:,:)                 ! specific kinetic energy
    real(wp), allocatable :: p_grad_acc_neg_l_u(:,:,:)        ! x-component of linear pressure gradient acceleration
    real(wp), allocatable :: p_grad_acc_neg_l_v(:,:,:)        ! y-component of linear pressure gradient acceleration
    real(wp), allocatable :: p_grad_acc_neg_l_w(:,:,:)        ! z-component of linear pressure gradient accelerationpgrad_acc_old
    real(wp), allocatable :: p_grad_acc_neg_nl_u(:,:,:)       ! x-component of nonlinear pressure gradient acceleration
    real(wp), allocatable :: p_grad_acc_neg_nl_v(:,:,:)       ! y-component of nonlinear pressure gradient acceleration
    real(wp), allocatable :: p_grad_acc_neg_nl_w(:,:,:)       ! z-component of nonlinear pressure gradient acceleration
    real(wp), allocatable :: p_grad_acc_old_u(:,:,:)          ! x-component of pressure gradient at old time step
    real(wp), allocatable :: p_grad_acc_old_v(:,:,:)          ! y-component of pressure gradient at old time step
    real(wp), allocatable :: p_grad_acc_old_w(:,:,:)          ! z-component of pressure gradient at old time step
    real(wp), allocatable :: v_squared_grad_x(:,:,:)          ! x-gradient of specific kinetic energy
    real(wp), allocatable :: v_squared_grad_y(:,:,:)          ! y-gradient of specific kinetic energy
    real(wp), allocatable :: v_squared_grad_z(:,:,:)          ! z-gradient of specific kinetic energy
    real(wp), allocatable :: pot_vort_tend_x(:,:,:)           ! tendency due to the vorticity flux term in x-direction
    real(wp), allocatable :: pot_vort_tend_y(:,:,:)           ! tendency due to the vorticity flux term in y-direction
    real(wp), allocatable :: pot_vort_tend_z(:,:,:)           ! tendency due to the vorticity flux term in z-direction
    real(wp), allocatable :: scalar_placeholder(:,:,:)        ! placeholder for scalar fields
    real(wp), allocatable :: temperature_gas(:,:,:)           ! temperature of the gas phase
    real(wp), allocatable :: u_placeholder(:,:,:)             ! placeholder for vector fields in x-direction
    real(wp), allocatable :: v_placeholder(:,:,:)             ! placeholder for vector fields in y-direction
    real(wp), allocatable :: w_placeholder(:,:,:)             ! placeholder for vector fields in z-direction
    real(wp), allocatable :: u_10(:,:)                        ! 10 m wind in x direction
    real(wp), allocatable :: v_10(:,:)                        ! 10 m wind in y direction
    real(wp), allocatable :: gust(:,:)                        ! gusts speed 10 m AGL
    real(wp), allocatable :: mslp(:,:)                        ! mean sea level pressure
    real(wp), allocatable :: t_2(:,:)                         ! 2 m temperature
    real(wp), allocatable :: z_eta_x(:,:,:)                   ! relative vorticity in x-direction
    real(wp), allocatable :: z_eta_y(:,:,:)                   ! relative vorticity in y-direction
    real(wp), allocatable :: z_eta_z(:,:,:)                   ! relative vorticity in z-direction
    real(wp), allocatable :: radiation_tendency(:,:,:)        ! power density due to radiation
    real(wp), allocatable :: scalar_flux_resistance(:,:)      ! surface flux resistance acting on scalar quantities
    real(wp), allocatable :: monin_obukhov_length(:,:)        ! Monin-Obukhov length
    real(wp), allocatable :: power_flux_density_sensible(:,:) ! power flux density acting on the surface due to sensible heat
    real(wp), allocatable :: power_flux_density_latent(:,:)   ! power flux density acting on the surface due to phase transitions
    real(wp), allocatable :: sfc_sw_in(:,:)                   ! shortwave radiation in the surface
    real(wp), allocatable :: sfc_lw_out(:,:)                  ! longwave radiation out of the surface
    real(wp), allocatable :: roughness_velocity(:,:)          ! roughness velocity
  
  
  end type t_diag
  
  type t_irrev
    
    ! type cotaining irreversible quantities
    real(wp), allocatable :: tke(:,:,:)                       ! specific turbulent kinetic energy
    real(wp), allocatable :: viscosity_molecular(:,:,:)       ! molecular diffusion coefficient
    real(wp), allocatable :: viscosity_coeff(:,:,:)           ! efficient viscosity (Eddies + molecular)
    real(wp), allocatable :: mom_diff_tend_x(:,:,:)           ! tendency due to momentum diffusion in x-direction
    real(wp), allocatable :: mom_diff_tend_y(:,:,:)           ! tendency due to momentum diffusion in y-direction
    real(wp), allocatable :: mom_diff_tend_z(:,:,:)           ! tendency due to momentum diffusion in z-direction
    real(wp), allocatable :: heating_diss(:,:,:)              ! dissipative heating power density
    real(wp), allocatable :: mass_source_rates(:,:,:,:)       ! mass source rates due to phase transitions and cloud physics
    real(wp), allocatable :: heat_source_rates(:,:,:,:)       ! heat source rates due to phase transitions and cloud physics
    real(wp)              :: max_diff_h_coeff_turb            ! maximum horizontal diffusion coefficient
  
  end type t_irrev
  
end module definitions






