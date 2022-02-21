! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module explicit_scalar_tendencies

  ! This module manages the calculation of the explicit component of the scalar tendencies.

  use definitions,          only: wp,t_grid,t_state,t_diag,t_irrev,t_tend
  use multiplications,      only: scalar_times_vector_h
  use divergence_operators, only: divv_h
  use run_nml,              only: nlins,ncols
  use constituents_nml,     only: no_of_constituents
  use phase_trans,          only: calc_h2otracers_source_rates
  use constituents_nml,     only: no_of_condensed_constituents,no_of_constituents

  implicit none
  
  private
  
  public :: expl_scalar_tend
  public :: moisturizer
  
  contains
  
  subroutine expl_scalar_tend(grid,state,tend,diag)
  
    type(t_grid),  intent(in)    :: grid  ! model grid
    type(t_state), intent(in)    :: state ! state with which to calculate the divergence
    type(t_tend),  intent(inout) :: tend  ! state which will contain the tendencies
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    
    ! local variables
    integer                      :: j_constituent ! loop variable
    
    ! explicit mass densities integration
    do j_constituent=1,no_of_constituents
      ! calculating the mass density flux
      call scalar_times_vector_h(state%rho(:,:,:,j_constituent),state%wind_u,state%wind_v,diag%u_placeholder,diag%v_placeholder)
      ! calculating the divergence of the mass density flux
      call divv_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:),grid)
      tend%rho(:,:,:,j_constituent) = -diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)
    enddo
    
    ! explicit potential temperature density integration
    ! calculating the potential temperature density flux
    diag%scalar_placeholder(:,:,:) = grid%theta_bg(:,:,:) + state%theta_pert(:,:,:)
    call scalar_times_vector_h(diag%scalar_placeholder,diag%u_placeholder,diag%v_placeholder, &
    diag%u_placeholder,diag%v_placeholder)
    ! calculating the divergence of the potential temperature density flux
    call divv_h(diag%u_placeholder,diag%v_placeholder,diag%scalar_placeholder(2:nlins+1,2:ncols+1,:),grid)
    tend%rhotheta(:,:,:) = -diag%scalar_placeholder(2:nlins+1,2:ncols+1,:)
        
  end subroutine
  
  subroutine moisturizer(state,irrev,dtime)
  
    ! This subroutine manages the calculation of the phase transition rates.
    
    type(t_state), intent(inout) :: state
    type(t_irrev), intent(inout) :: irrev
    real(wp),      intent(in)    :: dtime
      
    if (no_of_constituents > 1) then
    
      ! calculating the source rates
      call calc_h2otracers_source_rates()
      
      ! condensates
      state%rho(:,:,:,1:no_of_condensed_constituents) = state%rho(:,:,:,1:no_of_condensed_constituents) &
      + dtime*irrev%mass_source_rates(:,:,:,1:no_of_condensed_constituents)
      ! water vapour
      state%rho(:,:,:,no_of_constituents) = state%rho(:,:,:,no_of_constituents) &
      + dtime*irrev%mass_source_rates(:,:,:,no_of_constituents-1)
      
    endif
  
  end subroutine moisturizer

end module explicit_scalar_tendencies

















