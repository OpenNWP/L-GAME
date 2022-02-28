! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module manage_radiation_calls

  ! This manages the calls to RTE+RRTMGP.

  use definitions, only: t_grid,t_state,t_diag,t_irrev
  use run_nml,     only: nlins,ncols,nlays,wp
  use radiation,   only: calc_radiative_flux_convergence
  
  implicit none
  
  private
  
  public :: call_radiation
  
  contains
  
  subroutine call_radiation(state,grid,diag)
  
    ! This subroutine calls RTE+RRTMGP in a parallelized way.
    
    type(t_state),  intent(inout) :: state ! the state with which to calculate the radiative fluxes
    type(t_grid),   intent(inout) :: grid  ! the grid of the model
    type(t_diag),   intent(inout) :: diag  ! diagnostic quantities
    
    write(*,*) "Starting update of radiative fluxes ..."
    
    call calc_radiative_flux_convergence(grid%lat_geo_scalar,grid%lon_geo_scalar, &
    grid%z_scalar,grid%z_w,state%rho,diag%temperature_gas,diag%radiation_tendency, &
    state%temperature_soil(:,:,1),diag%sfc_sw_in,diag%sfc_lw_out,grid%sfc_albedo, &
    nlins*ncols,0._wp)
    
    write(*,*) "Update of radiative fluxes completed."
    
  end subroutine call_radiation

end module manage_radiation_calls








