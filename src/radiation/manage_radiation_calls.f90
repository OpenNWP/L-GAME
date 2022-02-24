! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module manage_radiation_calls

  ! This manages the calls to RTE+RRTMGP.

  use definitions,                only: t_grid,t_state,t_diag,t_irrev
  
  implicit none
  
  private
  
  public :: call_radiation
  
  contains
  
  subroutine call_radiation(state_old,grid,diag,irrev)
  
    ! calls RTE+RRTMGP
    
    type(t_state),  intent(inout) :: state_old          ! the state at the old timestep
    type(t_grid),   intent(inout) :: grid               ! the grid of the model
    type(t_diag),   intent(inout) :: diag               ! diagnostic quantities
    type(t_irrev),  intent(inout) :: irrev              ! irreversible quantities
    
    write(*,*) "Starting update of radiative fluxes ..."
    
    write(*,*) "Update of radiative fluxes completed."
    
  end subroutine call_radiation

end module manage_radiation_calls








