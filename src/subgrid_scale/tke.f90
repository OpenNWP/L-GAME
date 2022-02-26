! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module tke
  
  ! This module computes everything related to the turbulent kinetic energy (TKE).
  
  use definitions,        only: wp,t_state,t_irrev
  use run_nml,            only: nlins,ncols,nlays,dtime
  use derived_quantities, only: density_gas
  
  implicit none
  
  private
  
  public :: tke_update
  
  contains
  
  subroutine tke_update(state,irrev)
  
    ! This subroutine updates the specific turbulent kinetic energy (TKE), unit: J/kg.
  
    ! input arguments and output
    type(t_state), intent(in)    :: state    ! state
    type(t_irrev), intent(inout) :: irrev    ! irreversible quantities
    
    ! local variables
    integer                      :: ji,jk,jl ! loop variables
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          irrev%tke(ji,jk,jl) = irrev%tke(ji,jk,jl) + dtime*( &
          ! production of TKE through generation of resolved energy
          irrev%heating_diss(ji,jk,jl)/density_gas(state,ji,jk,jl) &
          )
          
          ! clipping negative values
          
          if (irrev%tke(ji,jk,jl) < 0._wp) then
            irrev%tke(ji,jk,jl) = 0._wp
          endif
          
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine tke_update
  
end module tke






