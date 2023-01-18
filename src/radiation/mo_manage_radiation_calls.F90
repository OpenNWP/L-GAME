! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_manage_radiation_calls
  
  ! This manages the calls to RTE+RRTMGP.
  
  use mo_definitions,    only: t_grid,t_state,t_diag
  use mo_run_nml,        only: nx,ny,wp
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
    real(wp)              :: sea_fraction       ! fraction of a grid cell that is covered by sea
    real(wp), allocatable :: temp_sfc_full(:,:) ! used surface temperature field
    integer               :: ji                 ! horizontal index
    integer               :: jk                 ! horizontal index
    
    write(*,*) "Starting update of radiative fluxes ..."
    
    allocate(temp_sfc_full(ny,nx))
    
    !$omp parallel do private(ji,jk,sea_fraction)
    do jk=1,nx
      do ji=1,ny
        sea_fraction = 1._wp-grid%land_fraction(ji,jk)-grid%lake_fraction(ji,jk)
        temp_sfc_full(ji,jk) = sea_fraction*diag%sst(ji,jk) + (1._wp-sea_fraction)*state%temperature_soil(ji,jk,1)
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji)
    do ji=1,ny
      call calc_radiative_flux_convergence(grid%lat_geo_scalar(ji,:),grid%lon_geo_scalar(ji,:),grid%z_scalar(ji,:,:), &
                                           grid%z_w(ji,:,:),state%rho(ji,:,:,:),diag%temperature(ji,:,:), &
                                           diag%radiation_tendency(ji,:,:),temp_sfc_full(ji,:), &
                                           diag%sfc_sw_in(ji,:),diag%sfc_lw_out(ji,:),grid%sfc_albedo(ji,:), &
                                           time_coordinate+0.5_wp*dtime_rad)
    enddo
    !$omp end parallel do
    
    deallocate(temp_sfc_full)
    
    write(*,*) "Update of radiative fluxes completed."
    
  end subroutine update_rad_fluxes

end module mo_manage_radiation_calls








