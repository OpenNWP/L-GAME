! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_manage_radiation_calls
  
  ! This manages the calls to RTE+RRTMGP.
  
  use mo_definitions,    only: t_grid,t_state,t_diag
  use mo_run_nml,        only: ny,wp
  use mo_rad_nml,        only: dtime_rad
  use mo_rrtmgp_coupler, only: calc_radiative_flux_convergence
  
  implicit none
  
  contains
  
  subroutine update_rad_fluxes(state,grid,diag,time_coordinate)
    
    ! This subroutine calls RTE+RRTMGP in a parallelized way.
    
    type(t_state), intent(in)    :: state           ! the state with which to calculate the radiative fluxes
    type(t_grid),  intent(inout) :: grid            ! the grid of the model
    type(t_diag),  intent(inout) :: diag            ! diagnostic quantities
    real(wp),      intent(in)    :: time_coordinate ! Unix time
    
    ! local variables
    integer :: ji ! horizontal index
    
    write(*,*) "Starting update of radiative fluxes ..."
    
    !$omp parallel do private(ji)
    do ji=1,ny
      call calc_radiative_flux_convergence(grid%lat_geo_scalar(ji,:),grid%lon_geo_scalar(ji,:),grid%z_scalar(ji,:,:), &
                                           grid%z_w(ji,:,:),state%rho(ji,:,:,:),diag%temperature(ji,:,:), &
                                           diag%radiation_tendency(ji,:,:),state%temperature_soil(ji,:,1), &
                                           diag%sfc_sw_in(ji,:),diag%sfc_lw_out(ji,:),grid%sfc_albedo(ji,:), &
                                           time_coordinate+0.5_wp*dtime_rad)
    enddo
    !$omp end parallel do
    
    write(*,*) "Update of radiative fluxes completed."
    
  end subroutine update_rad_fluxes

end module mo_manage_radiation_calls








