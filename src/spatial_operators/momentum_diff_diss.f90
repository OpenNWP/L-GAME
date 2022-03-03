! This ! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module momentum_diff_diss

  ! This module handles momentum diffusion and dissipation.
  
  use definitions,           only: t_grid,t_diag,t_irrev,t_state
  use divergence_operators,  only: div_h,add_vertical_div
  use gradient_operators,    only: grad_hor,grad_vert_cov
  use run_nml,               only: nlins,ncols,nlays,wp
  use inner_product,         only: inner
  use derived_quantities,    only: density_gas
  use effective_diff_coeffs, only: hori_div_viscosity,vert_vert_mom_viscosity
  use multiplications,       only: scalar_times_scalar
  
  implicit none
  
  private
  
  public :: mom_diff_h
  public :: mom_diff_v
  public :: simple_dissipation_rate
  
  contains
  
  subroutine mom_diff_h(state,diag,irrev,grid)
  
    ! This subroutine handles horizontal momentum diffusion.
    
    type(t_state), intent(in)    :: state ! state with which to calculate the horizontal diffusion
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! calculating the divergence of the horizontal wind field
    call div_h(state%wind_u,state%wind_v,diag%scalar_placeholder,grid)
    
    ! computing the relevant diffusion coefficient
    call hori_div_viscosity(state,diag,diag%scalar_placeholder,irrev,grid)
    
    ! multiplying the divergence by the diffusion coefficient acting on divergent movements
    call scalar_times_scalar(irrev%viscosity_coeff_div,diag%scalar_placeholder,diag%scalar_placeholder)
    
    ! calculating the horizontal gradient of the divergence
    call grad_hor(diag%scalar_placeholder,irrev%mom_diff_tend_x,irrev%mom_diff_tend_y,irrev%mom_diff_tend_z,grid)
  
  end subroutine mom_diff_h

  subroutine mom_diff_v(state,diag,irrev,grid)
  
    ! This subroutine handles vertical momentum diffusion.
    
    type(t_state), intent(in)    :: state ! state with which to calculate the horizontal diffusion
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    	
    ! 2.) vertical diffusion of vertical velocity
    ! -------------------------------------------
    ! resetting the placeholder field
    !$OMP PARALLEL
    !$OMP WORKSHARE
    diag%scalar_placeholder = 0._wp
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! computing something like dw/dz
    call add_vertical_div(state%wind_w,diag%scalar_placeholder,grid)
    ! computing and multiplying by the respective diffusion coefficient
    call vert_vert_mom_viscosity(state,diag,irrev,grid)
    ! taking the second derivative to compute the diffusive tendency
    call grad_vert_cov(diag%scalar_placeholder,irrev%mom_diff_tend_z,grid)

  end subroutine mom_diff_v
  
  subroutine simple_dissipation_rate(state,irrev,grid)
  
    ! This subroutine calculates a simplified dissipation rate.
    
    ! input arguments and output
    type(t_state), intent(in)    :: state ! the state with which to calculate the dissipation rates
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! calculating the inner product of the momentum diffusion acceleration and the wind
    call inner(state%wind_u,state%wind_v,state%wind_w, &
    irrev%mom_diff_tend_x,irrev%mom_diff_tend_y,irrev%mom_diff_tend_z,irrev%heating_diss,grid)
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          irrev%heating_diss(ji,jk,jl) = -density_gas(state,ji,jk,jl)*irrev%heating_diss(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine simple_dissipation_rate

end module momentum_diff_diss
















