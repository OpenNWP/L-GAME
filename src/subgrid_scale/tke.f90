! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module tke
  
  ! This module computes everything related to the turbulent kinetic energy (TKE).
  
  use definitions,        only: wp,t_state,t_diag,t_irrev,t_grid
  use run_nml,            only: ny,nx,nlays,dtime
  use derived_quantities, only: density_gas
  use constants,          only: M_PI
  use gradient_operators, only: grad
  use inner_product,      only: inner
  
  implicit none
  
  private
  
  public :: tke_update
  
  contains
  
  subroutine tke_update(state,diag,irrev,grid)
  
    ! This subroutine updates the specific turbulent kinetic energy (TKE), unit: J/kg.
  
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    real(wp) :: decay_constant ! defines how quickly the TKE decays
    integer  :: ji,jk,jl       ! loop variables
    
    ! computing the advection
    call grad(irrev%tke,diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
    call inner(state%wind_u,state%wind_v,state%wind_w,diag%u_placeholder,diag%v_placeholder, &
    diag%w_placeholder,diag%scalar_placeholder,grid)
    
    !$omp parallel do private(ji,jk,jl,decay_constant)
    do ji=1,ny
      do jk=1,nx
        do jl=1,nlays
          ! calculating the decay constant
          decay_constant = 8._wp*M_PI**2/grid%mean_velocity_area*(irrev%viscosity_coeff_div(ji,jk,jl) &
          + irrev%viscosity_coeff_curl(ji,jk,jl))/density_gas(state,ji,jk,jl)
          
          irrev%tke(ji,jk,jl) = irrev%tke(ji,jk,jl) + dtime*( &
          ! advection
          -diag%scalar_placeholder(ji,jk,jl) &
          ! production of TKE through generation of resolved energy
          + irrev%heating_diss(ji,jk,jl)/density_gas(state,ji,jk,jl) &
          ! decay through molecular dissipation
          - decay_constant*irrev%tke(ji,jk,jl) &
          )
          
          ! clipping negative values
          if (irrev%tke(ji,jk,jl)<0._wp) then
            irrev%tke(ji,jk,jl) = 0._wp
          endif
          
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine tke_update
  
end module tke






