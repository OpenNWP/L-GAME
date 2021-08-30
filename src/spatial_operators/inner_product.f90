! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module inner_product

  use run_nml,     only: nlins,ncols,nlays
  use definitions, only: t_state,t_diag,t_grid,wp
  
  implicit none
  
  private
  
  public :: kinetic_energy
  
  contains

  subroutine kinetic_energy(state,diag,grid)
  
    ! This subroutine calculates the specific kinetic energy.
    
    type(t_state), intent(in)    :: state ! state to use for calculating e_kin
    type(t_diag),  intent(inout) :: diag  ! diagnostic properties (e_kin is a diagnostic property)
    type(t_grid),  intent(in)    :: grid  ! grid properties
    ! local variables
    integer :: ji,jk,jl                   ! loop indices
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          diag%e_kin(ji,jk,jl) = &
          grid%inner_product_weights(ji,jk,jl,1)*state%wind_u(ji+1,jk+1,jl  )**2 &
          + grid%inner_product_weights(ji,jk,jl,2)*state%wind_v(ji+1,jk+1,jl  )**2 &
          + grid%inner_product_weights(ji,jk,jl,3)*state%wind_u(ji+1,jk  ,jl  )**2 &
          + grid%inner_product_weights(ji,jk,jl,4)*state%wind_v(ji  ,jk+1,jl  )**2 &
          + grid%inner_product_weights(ji,jk,jl,5)*state%wind_w(ji+1,jk+1,jl  )**2 &
          + grid%inner_product_weights(ji,jk,jl,6)*state%wind_w(ji+1,jk+1,jl+1)**2
          diag%e_kin(ji,jk,jl) = 0.5_wp*diag%e_kin(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine kinetic_energy

end module inner_product
