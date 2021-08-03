! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module io

	use definitions,    only: t_state,wp,t_diag,t_grid,t_bg
	use netcdf
	use run_nml,        only: nlins,ncols,nlays
	use thermodynamics, only: gas_constant_diagnostics

	implicit none
	
	private
	
	public :: read_init
	public :: write_output
	
	contains
	
	subroutine read_init(state)
		! reads the initial state of the model calculation
		
		type(t_state), intent(inout) :: state
		
	end subroutine read_init
	
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
		call check(nf90_def_dim(ncid,"z",nlays,z_dimid))

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
		call check(nf90_put_var(ncid,varid_t,diag%scalar_placeholder))
		diag%scalar_placeholder(:,:,:) =  state%theta(:,:,:)*(bg%exner(:,:,:) + state%exner_pert(:,:,:))
		
		! writing the data to the output file
		! 3D pressure
		call check(nf90_put_var(ncid,varid_p,diag%scalar_placeholder))
		diag%scalar_placeholder(:,:,:) = state%rho(:,:,:)*gas_constant_diagnostics(1)*diag%scalar_placeholder(:,:,:)
		
		! 3D u wind
		do jk=1,ncols
			diag%scalar_placeholder(:,jk,:) = 0.5_wp*(state%wind_u(:,jk,:)+state%wind_u(:,jk+1,:))
		enddo
		call check(nf90_put_var(ncid,varid_u,diag%scalar_placeholder))
		
		! 3D v wind
		do ji=1,nlins
			diag%scalar_placeholder(ji,:,:) = 0.5_wp*(state%wind_u(ji,:,:)+state%wind_u(ji+1,:,:))
		enddo
		call check(nf90_put_var(ncid,varid_v,diag%scalar_placeholder))
		
		! 3D w wind
		call check(nf90_put_var(ncid,varid_w,diag%scalar_placeholder))
		do jl=1,nlays
			upper_weight(:,:) = (grid%z_geo_scal(:,:,jl) - grid%z_geo_w(:,:,jl+1))-(grid%z_geo_w(:,:,jl) - grid%z_geo_w(:,:,jl+1))
			diag%scalar_placeholder(:,:,jl) = upper_weight*state%wind_w(:,:,jl)+(1-upper_weight)*state%wind_w(:,:,jl+1)
		enddo
  
		! closing the netcdf file
		call check(nf90_close(ncid))
		
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
	
	subroutine check(i_status)
		integer, intent(in) :: i_status

		if(i_status /= nf90_noerr) then 
			print *, trim(nf90_strerror(i_status))
			stop "Netcdf threw an error."
		end if
	end subroutine check  

end module io








