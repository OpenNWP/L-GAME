! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module boundaries

  ! This module handles everything dealing with boundary conditions.

  use definitions,      only: t_state,t_bc,t_grid,wp
  use run_nml,          only: nlins,ncols,nlays
  use bc_nml,           only: n_swamp
  use constants,        only: M_PI,p_0
  use constituents_nml, only: no_of_condensed_constituents
  use dictionary,       only: spec_heat_capacities_v_gas,specific_gas_constants

  implicit none
  
  private
  
  public :: update_boundaries
  public :: rescale_tend
  public :: setup_bc_factor
  
  contains
  
  subroutine update_boundaries(state,bc,grid)
  
    ! updates the boundary conditions
    
    ! input arguments and output
    type(t_state), intent(inout) :: state
    type(t_bc),    intent(inout) :: bc
    type(t_grid),  intent(in)    :: grid
    
    ! local variables
    real(wp) :: old_weight,new_weight ! time interpolation weights
    real(wp) :: c_v                   ! specific heat capacity at constant volume
    real(wp) :: r_d                   ! individual gas constant of dry air
    integer  :: ji,jk,jl              ! loop indicies
    
    c_v = spec_heat_capacities_v_gas(0)
    r_d = specific_gas_constants(0)
    
    ! linear combination of the model state and the boundary conditions
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          state%rho(ji,jk,jl,:) = bc%scalar_bc_factor(ji,jk)*(old_weight*bc%rho(ji,jk,jl,:,1) &
          + new_weight*bc%rho(ji,jk,jl,:,2)) + (1._wp - bc%scalar_bc_factor(ji,jk))*state%rho(ji,jk,jl,:)
          state%condensed_rho_t(ji,jk,jl,:) = bc%scalar_bc_factor(ji,jk)*(old_weight*bc%condensed_rho_t(ji,jk,jl,:,1) &
          + new_weight*bc%condensed_rho_t(ji,jk,jl,:,2)) + (1._wp - bc%scalar_bc_factor(ji,jk))*state%condensed_rho_t(ji,jk,jl,:)
          state%rhotheta(ji,jk,jl) = bc%scalar_bc_factor(ji,jk)*(old_weight*bc%rhotheta(ji,jk,jl,1) &
          + new_weight*bc%rhotheta(ji,jk,jl,2)) + (1._wp - bc%scalar_bc_factor(ji,jk))*state%rhotheta(ji,jk,jl)
          state%wind_w(ji,jk,jl) = bc%scalar_bc_factor(ji,jk)*(old_weight*bc%wind_w(ji,jk,jl,1) &
          + new_weight*bc%wind_w(ji,jk,jl,2)) + (1._wp - bc%scalar_bc_factor(ji,jk))*state%wind_w(ji,jk,jl)
        enddo
        state%wind_w(ji,jk,nlays+1) = bc%scalar_bc_factor(ji,jk)*(old_weight*bc%wind_w(ji,jk,nlays+1,1) &
        + new_weight*bc%wind_w(ji,jk,nlays+1,2)) + (1._wp - bc%scalar_bc_factor(ji,jk))*state%wind_w(ji,jk,nlays+1)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols+1
        do jl=1,nlays
          state%wind_u(ji,jk,jl) = bc%u_bc_factor(ji,jk)*(old_weight*bc%wind_u(ji,jk,jl,1) &
          + new_weight*bc%wind_u(ji,jk,jl,2)) + (1._wp - bc%u_bc_factor(ji,jk))*state%wind_u(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins+1
      do jk=1,ncols
        do jl=1,nlays
          state%wind_v(ji,jk,jl) = bc%v_bc_factor(ji,jk)*(old_weight*bc%wind_v(ji,jk,jl,1) &
          + new_weight*bc%wind_v(ji,jk,jl,2)) + (1._wp - bc%v_bc_factor(ji,jk))*state%wind_v(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! setting the potential temperature peturbation
    !$OMP PARALLEL
    !$OMP WORKSHARE
    state%theta_pert = state%rhotheta/state%rho(:,:,:,no_of_condensed_constituents+1) - grid%theta_bg
    state%exner_pert = (r_d*state%rhotheta/p_0)**(r_d/c_v) - grid%exner_bg
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    
  end subroutine update_boundaries
  
  subroutine rescale_tend
  
    ! rescales the tendencies to account for boundary conditions
  
  end subroutine rescale_tend
  
  subroutine setup_bc_factor(bc)
  
    ! This subroutine calculates the boundary conditions rescale factors.
  
    ! argument and output
    type(t_bc), intent(inout) :: bc
    
    ! local variables
    real(wp) :: dist_from_boundary ! index distance from the boundaary of the domain
    integer  :: ji,jk              ! loop indices
    
    ! rescale factor for scalar fields
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins
      do jk=1,ncols
        dist_from_boundary = min(ji-1,jk-1,nlins-ji,ncols-jk)
        bc%scalar_bc_factor(ji,jk) = max(1._wp - dist_from_boundary/n_swamp,0._wp)
        bc%scalar_bc_factor(ji,jk) = sin(0.5_wp*M_PI*bc%scalar_bc_factor(ji,jk))**2
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! u rescale factor
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins
      do jk=1,ncols+1
        dist_from_boundary = min(ji-1,jk-1,nlins-ji,ncols+1-jk)
        bc%u_bc_factor(ji,jk) = max(1._wp - dist_from_boundary/n_swamp,0._wp)
        bc%u_bc_factor(ji,jk) = sin(0.5_wp*M_PI*bc%u_bc_factor(ji,jk))**2
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! v rescale factor
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins+1
      do jk=1,ncols
        dist_from_boundary = min(ji-1,jk-1,nlins+1-ji,ncols-jk)
        bc%v_bc_factor(ji,jk) = max(1._wp - dist_from_boundary/n_swamp,0._wp)
        bc%v_bc_factor(ji,jk) = sin(0.5_wp*M_PI*bc%v_bc_factor(ji,jk))**2
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine setup_bc_factor

end module boundaries








