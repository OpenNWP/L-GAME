! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_tke
  
  ! This module computes everything related to the turbulent kinetic energy (TKE).
  
  use mo_definitions,        only: wp,t_state,t_diag,t_grid
  use mo_run_nml,            only: ny,nx,n_layers,dtime
  use mo_constituents_nml,   only: n_condensed_constituents,n_constituents
  use mo_constants,          only: M_PI
  use mo_gradient_operators, only: grad_hor,grad_vert
  use mo_inner_product,      only: inner_product
  
  implicit none
  
  contains
  
  subroutine tke_update(state,diag,grid)
    
    ! This subroutine updates the specific turbulent kinetic energy (TKE), unit: J/kg.
    
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    real(wp) :: decay_constant ! defines how quickly the TKE decays
    integer  :: ji             ! horizontal index
    integer  :: jk             ! horizontal index
    integer  :: jl             ! layer index
    
    ! computing the advection
    call grad_vert(diag%tke,diag%w_placeholder,grid)
    call grad_hor(diag%tke,diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
    call inner_product(state%wind_u,state%wind_v,state%wind_w,diag%u_placeholder,diag%v_placeholder, &
    diag%w_placeholder,diag%scalar_placeholder,grid)
    
    !$omp parallel do private(ji,jk,jl,decay_constant)
    do jl=1,n_layers
      do jk=1,nx
        do ji=1,ny
          ! calculating the decay constant
          decay_constant = 8._wp*M_PI**2/grid%mean_velocity_area*(diag%viscosity_coeff_div(ji,jk,jl) &
          + diag%viscosity_coeff_curl(ji,jk,jl))/state%rho(ji,jk,jl,n_condensed_constituents+1)
          
          diag%tke(ji,jk,jl) = diag%tke(ji,jk,jl) + dtime*( &
          ! advection
          -diag%scalar_placeholder(ji,jk,jl) &
          ! production of TKE through generation of resolved energy
          + diag%heating_diss(ji,jk,jl)/state%rho(ji,jk,jl,n_condensed_constituents+1) &
          ! decay through molecular dissipation
          - decay_constant*diag%tke(ji,jk,jl) &
          )
          
          ! clipping negative values
          if (diag%tke(ji,jk,jl)<0._wp) then
            diag%tke(ji,jk,jl) = 0._wp
          endif
          
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine tke_update
  
end module mo_tke






