! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module io

	! This module handles everything dealing with IO.

	use definitions,    only: t_state,wp,t_diag,t_grid,t_bg
	use netcdf
	use run_nml,        only: nlins,ncols,nlays,scenario,p_0
	use thermodynamics, only: spec_heat_cap_diagnostics_v,gas_constant_diagnostics,spec_heat_cap_diagnostics_p
	use grid_generator, only: bg_temp,bg_pres,geopot

	implicit none
	
	private
	
	public :: ideal
	public :: restart
	public :: var_3d
	public :: var_4d
	public :: write_output
	
	contains
	
	subroutine ideal(state,diag,bg,grid)
	
		! sets the initial state of the model calculation i terms if analytic functions
		
		type(t_state), intent(inout) :: state ! state to write the initial state to
		type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
		type(t_bg),    intent(in)    :: bg    ! background state
		type(t_grid),  intent(in)    :: grid  ! model grid
		
		! local variables
		integer                      :: ji,jk,jl                           ! loop indices
		real(wp)                     :: pres_lowest_layer(nlins+2,ncols+2) ! pressure in the lowest layer
			
		select case (trim(scenario))
		
			case("standard")
			
			! This test case is the standard atmosphere.
			
			state%wind_u(:,:,:) = 0._wp
			state%wind_v(:,:,:) = 0._wp
			
			do ji=1,nlins+2
				do jk=1,ncols+2
					do jl=1,nlays
						diag%scalar_placeholder(ji,jk,jl) = bg_temp(grid%z_geo_scal(ji,jk,jl))
					enddo
					pres_lowest_layer(ji,jk) = bg_pres(grid%z_geo_scal(ji,jk,nlays))
				enddo
			enddo
		
		endselect
		
		call unessential_init(state,diag,bg,grid,pres_lowest_layer)
		
	end subroutine ideal
	
	subroutine restart(state,diag,bg,grid)
	
		! reads the initial state of the model calculation from a netcdf file
		
		type(t_state), intent(inout) :: state ! state to write the initial state to
		type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
		type(t_bg),    intent(in)    :: bg    ! background state
		type(t_grid),  intent(in)    :: grid  ! model grid
		
		! local variables
		real(wp)                     :: pres_lowest_layer(nlins+2,ncols+2) ! pressure in the lowest layer
		
		call unessential_init(state,diag,bg,grid,pres_lowest_layer)
		
	end subroutine restart
	
	subroutine var_3d(state,diag,bg,grid)
	
		! three-dimensional variational data assimilation
		
		type(t_state), intent(inout) :: state ! state to write the initial state to
		type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
		type(t_bg),    intent(in)    :: bg    ! background state
		type(t_grid),  intent(in)    :: grid  ! model grid
		
		! local variables
		real(wp)                     :: pres_lowest_layer(nlins+2,ncols+2) ! pressure in the lowest layer
		
		call unessential_init(state,diag,bg,grid,pres_lowest_layer)
	
	end subroutine var_3d
	
	subroutine var_4d(state,diag,bg,grid)
	
		! four-dimensional variational data assimilation
		
		type(t_state), intent(inout) :: state ! state to write the initial state to
		type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
		type(t_bg),    intent(in)    :: bg    ! background state
		type(t_grid),  intent(in)    :: grid  ! model grid
		
		! local variables
		real(wp)                     :: pres_lowest_layer(nlins+2,ncols+2) ! pressure in the lowest layer
		
		call unessential_init(state,diag,bg,grid,pres_lowest_layer)
	
	end subroutine var_4d
	
	subroutine bc()
	
		! sets the boundary conditions
		
	end subroutine bc
	
	subroutine write_output(state,diag,time_since_init_min,grid,bg)
		! reads out the state of the model atmosphere
		! at a single timestep to a netcdf file
		
		type(t_state), intent(in)    :: state               ! state to write out
		type(t_diag),  intent(inout) :: diag                ! diagnostic quantities
		integer,       intent(in)    :: time_since_init_min ! time in minutes since init
		type(t_grid),  intent(in)    :: grid                ! model grid
		type(t_bg),    intent(in)    :: bg                  ! background state
		
		! local variables
		integer                   :: ncid                      ! ID of the netcdf file
		integer                   :: x_dimid                   ! ID of the x dimension
		integer                   :: y_dimid                   ! ID of the y dimension
		integer                   :: z_dimid                   ! ID of the z dimension
		integer                   :: dimids_2d(2)              ! dimensions of horizontal fields
		integer                   :: dimids_3d(3)              ! dimensions of 3D fields
		integer                   :: varid_p                   ! variable ID of the 3D pressure field
		integer                   :: varid_t                   ! variable ID of the 3D temperature field
		integer                   :: varid_u                   ! variable ID of the 3D u wind field
		integer                   :: varid_v                   ! variable ID of the 3D v wind field
		integer                   :: varid_w                   ! variable ID of the 3D w wind field
	    character(len=32)         :: filename                  ! output filename
	    integer                   :: ji,jk,jl                  ! line indices
	    real(wp)                  :: upper_weight(nlins,ncols) ! interpolation weights
		
		! creating the netcdf file
		write(filename,"(I10,A3)") time_since_init_min,".nc"
		call check(nf90_create(trim(filename),NF90_CLOBBER,ncid))
		
		! defining the dimensions
		call check(nf90_def_dim(ncid,"x",ncols,x_dimid))
		call check(nf90_def_dim(ncid,"y",nlins,y_dimid))
		call check(nf90_def_dim(ncid,"z",nlays,  z_dimid))

		! setting the dimension ID arrays
		! 2D
		dimids_2d(1) = y_dimid
		dimids_2d(2) = x_dimid
		! 3D
		dimids_3d(1) = y_dimid
		dimids_3d(2) = x_dimid
		dimids_3d(3) = z_dimid

		! Define the variable. The type of the variable in this case is
		! NF90_INT (4-byte integer).
		call check(nf90_def_var(ncid,"p",NF90_REAL,dimids_3d,varid_p))
		call check(nf90_def_var(ncid,"T",NF90_REAL,dimids_3d,varid_t))
		call check(nf90_def_var(ncid,"u",NF90_REAL,dimids_3d,varid_u))
		call check(nf90_def_var(ncid,"v",NF90_REAL,dimids_3d,varid_v))
		call check(nf90_def_var(ncid,"w",NF90_REAL,dimids_3d,varid_w))
  
		! ending the definition section
		call check(nf90_enddef(ncid))
	
		! 3D temperature
		call check(nf90_put_var(ncid,varid_t,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
		diag%scalar_placeholder(2:nlins+1,2:ncols+1,:) =  (bg%theta(2:nlins+1,2:ncols+1,:) &
		+ state%theta_pert(2:nlins+1,2:ncols+1,:)) &
		*(bg%exner(2:nlins+1,2:ncols+1,:) + state%exner_pert(2:nlins+1,2:ncols+1,:))
		
		! writing the data to the output file
		! 3D pressure
		call check(nf90_put_var(ncid,varid_p,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
		diag%scalar_placeholder(2:nlins+1,2:ncols+1,:) = state%rho(2:nlins+1,2:ncols+1,:) &
		*gas_constant_diagnostics(1)*diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)
		
		! 3D u wind
		do jk=1,ncols
			diag%scalar_placeholder(2:nlins+1,jk+1,:) = 0.5_wp*(state%wind_u(2:nlins+1,jk,:)+state%wind_u(2:nlins+1,jk+1,:))
		enddo
		call check(nf90_put_var(ncid,varid_u,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
		
		! 3D v wind
		do ji=1,nlins
			diag%scalar_placeholder(ji+1,2:ncols+1,:) = 0.5_wp*(state%wind_v(ji,2:ncols+1,:)+state%wind_v(ji+1,2:ncols+1,:))
		enddo
		call check(nf90_put_var(ncid,varid_v,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
		
		! 3D w wind
		do jl=1,nlays
			upper_weight(:,:) = (grid%z_geo_scal(2:nlins+1,2:ncols+1,jl) -&
			grid%z_geo_w(2:nlins+1,2:ncols+1,jl+1))-(grid%z_geo_w(2:nlins+1,2:ncols+1,jl)  &
			- grid%z_geo_w(2:nlins+1,2:ncols+1,jl+1))
			diag%scalar_placeholder(2:nlins+1,2:ncols+1,jl) = &
			upper_weight*state%wind_w(2:nlins+1,2:ncols+1,jl)+(1-upper_weight)*state%wind_w(2:nlins+1,2:ncols+1,jl+1)
		enddo
		call check(nf90_put_var(ncid,varid_w,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)))
  
		! closing the netcdf file
		call check(nf90_close(ncid))
		
	end subroutine write_output
	
	subroutine unessential_init(state,diag,bg,grid,pres_lowest_layer)
	
		! setting the unessential quantities of an initial state
		
		type(t_state), intent(inout) :: state                  ! state to work with
		type(t_diag),  intent(in)    :: diag                   ! diagnostic quantities
		type(t_bg),    intent(in)    :: bg                     ! background state
		type(t_grid),  intent(in)    :: grid                   ! model grid
		real(wp),      intent(in)    :: pres_lowest_layer(:,:) ! pressure in the lowest layer
		
		! local variables
		integer                      :: ji,jk,jl          ! loop indices
		real(wp)                     :: b,c               ! abbreviations needed for the hydrostatic initialization routine
		real(wp)                     :: temperature       ! single temperature value
		real(wp)                     :: pressure          ! single pressure value
		
		! integrating the hydrostatic initial state according to the given temperature field and pressure in the lowest layer
		do ji=1,nlins+2
			do jk=1,ncols+2	
				! integrating from bottom to top
				do jl=nlays,1,-1
					temperature = diag%scalar_placeholder(ji,jk,jl)
					! lowest layer
					if (jl == nlays) then
						pressure    = pres_lowest_layer(ji,jk)
						state%theta_pert(ji,jk,jl) = temperature*(pressure/p_0)**(gas_constant_diagnostics(1)/spec_heat_cap_diagnostics_p(1))
						state%exner_pert(ji,jk,jl) = temperature/state%theta_pert(ji,jk,jl)
					! other layers
					else
						! solving a quadratic equation for the Exner pressure
						b = -0.5_wp*state%exner_pert(ji,jk,jl+1)/bg_temp(grid%z_geo_scal(ji,jk,jl+1)) &
						*(temperature - bg_temp(grid%z_geo_scal(ji,jk,jl+1)) + 2.0_wp/ &
						spec_heat_cap_diagnostics_p(1)*(geopot(grid%z_geo_scal(ji,jk,jl)) - geopot(grid%z_geo_scal(ji,jk,jl+1))))
						c = state%exner_pert(ji,jk,jl+1)**2*temperature/bg_temp(grid%z_geo_scal(ji,jk,jl+1))
						state%exner_pert(ji,jk,jl) = b+sqrt(b**2+c)
						state%theta_pert(ji,jk,jl) = temperature/state%exner_pert(ji,jk,jl)
					endif
				enddo
			enddo
		enddo
		
		! substracting the background state
		state%theta_pert(:,:,:) = state%theta_pert(:,:,:) - bg%theta(:,:,:)
		state%exner_pert(:,:,:) = state%exner_pert(:,:,:) - bg%exner(:,:,:)
		
		! potential temoerature density
		state%rhotheta(:,:,:) = p_0/gas_constant_diagnostics(1)*(bg%exner(:,:,:) + state%exner_pert(:,:,:)) &
		**(spec_heat_cap_diagnostics_v(1)/gas_constant_diagnostics(1))
		! potential temperature
		state%theta_pert(:,:,:) = diag%scalar_placeholder(:,:,:)/(bg%exner(:,:,:) + state%exner_pert(:,:,:)) - bg%theta(:,:,:)
		! density
		state%rho(:,:,:)      = state%rhotheta(:,:,:)/(bg%theta(:,:,:) + state%theta_pert(:,:,:))
		
		! vertical wind velocity
		state%wind_w(:,:,:) = 0._wp
		
	end subroutine unessential_init
	
	subroutine check(i_status)
		integer, intent(in) :: i_status

		if(i_status /= nf90_noerr) then 
			print *, trim(nf90_strerror(i_status))
			stop "Netcdf threw an error."
		end if
	end subroutine check  

end module io








