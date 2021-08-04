! This source file is pa! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module vorticities

	! This module contians the calculation of the vorticities.

	use definitions, only: t_state,t_diag,t_grid,wp
	use run_nml,     only: nlins,ncols,nlays
	
	implicit none
	
	private
	
	public :: calc_pot_vort
	
	contains

	subroutine rel_vort(state,diag,grid)
	
		! This subroutine calculates the relative vorticity.
		
		type(t_state), intent(in)    :: state    ! state to work with
		type(t_diag),  intent(inout) :: diag     ! diagnostic quantities
		type(t_grid),  intent(in)    :: grid     ! model grid
		
		! local variables
		integer                      :: ji,jl,jk ! loop indices
		
		! calculating the relative vorticity in x-direction
		do ji=1,nlins+1
			do jk=1,ncols
				do jl=2,nlays
					diag%zeta_x(ji,jk,jl) = &
					+ grid%dz(ji  ,jk+1,jl  )*state%wind_w(ji  ,jk+1,jl  ) &
					+ grid%dy(ji  ,jk+1,jl-1)*state%wind_v(ji  ,jk+1,jl-1) &
					- grid%dz(ji+1,jk+1,jl  )*state%wind_w(ji+1,jk+1,jl  ) &
					- grid%dy(ji  ,jk+1,jl  )*state%wind_v(ji  ,jk+1,jl  )
				enddo
				! At the surface, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
				diag%zeta_x(ji,jk,nlays+1) = grid%dy(ji,jk+1,jl-1)*state%wind_v(ji,jk+1,nlays)
			enddo
		enddo
		! At the TOA, the horizontal vorticity is assumed to have no vertical shear.
		diag%zeta_x(:,:,1) = diag%zeta_x(:,:,2)
		! dividing by the area
		diag%zeta_x(:,:,:) = grid%area_dual_x(:,:,:)
		
		! calculating the relative vorticity in y-direction
		do ji=1,nlins
			do jk=1,ncols+1
				do jl=2,nlays
					diag%zeta_y(ji,jk,jl) = &
					+ grid%dz(ji+1,jk+1,jl  )*state%wind_w(ji+1,jk+1,jl  ) &
					- grid%dy(ji+1,jk  ,jl-1)*state%wind_v(ji+1,jk  ,jl-1) &
					- grid%dz(ji+1,jk  ,jl  )*state%wind_w(ji+1,jk  ,jl  ) &
					+ grid%dy(ji+1,jk  ,jl  )*state%wind_v(ji+1,jk  ,jl  )
				enddo
				! At the surface, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
				diag%zeta_y(ji,jk,nlays+1) = -grid%dy(ji+1,jk,nlays)*state%wind_v(ji+1,jk,jl-1)
			enddo
		enddo
		! At the TOA, the horizontal vorticity is assumed to have no vertical shear.
		diag%zeta_y(:,:,1) = diag%zeta_y(:,:,2)
		! dividing by the area
		diag%zeta_y(:,:,:) = grid%area_dual_y(:,:,:)
		
		! calculating the relative vorticity in z-direction
		do ji=1,nlins+1
			do jk=1,ncols+1
				do jl=1,nlays
					diag%zeta_z(ji,jk,jl) = &
					+ grid%dy(ji,  jk+1,jl)*state%wind_v(ji,  jk+1,jl) &
					- grid%dx(ji+1,jk,  jl)*state%wind_u(ji+1,jk,  jl) &
					- grid%dy(ji,  jk,  jl)*state%wind_v(ji,  jk,  jl) &
					+ grid%dx(ji,  jk,  jl)*state%wind_u(ji,  jk,  jl)
				enddo
			enddo
		enddo
		! dividing by the area
		diag%zeta_z(:,:,:) = grid%area_dual_z(:,:,:)
			
	end subroutine rel_vort
	
	subroutine calc_pot_vort(state,diag,grid)
	
		! This subroutine calculates the relative vorticity.
		
		type(t_state), intent(in)    :: state    ! state to work with
		type(t_diag),  intent(inout) :: diag     ! diagnostic quantities
		type(t_grid),  intent(in)    :: grid     ! model grid
	
		! calculating the relative vorticity
		call rel_vort(state,diag,grid)
		! adding the Coriolis vector to the relative vorticity to obtain the absolute vorticity
		! dividing by the averaged density to obtain the "potential vorticity"
	
	end subroutine calc_pot_vort

end module vorticities






