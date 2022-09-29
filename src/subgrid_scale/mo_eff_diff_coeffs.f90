! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_eff_diff_coeffs
  
  ! This module computes the effective diffusion coefficients.
  
  use mo_run_nml,              only: ny,nx,n_layers,n_levels,dtime,eff_hor_res
  use mo_definitions,          only: wp,t_state,t_diag,t_grid
  use mo_diff_nml,             only: lmom_diff_h,ltemp_diff_h
  use mo_derived,              only: calc_diffusion_coeff
  use mo_constituents_nml,     only: n_constituents,n_condensed_constituents
  use mo_tke,                  only: tke_update
  use mo_divergence_operators, only: div_h
  use mo_gradient_operators,   only: grad_vert
  use mo_multiplications,      only: scalar_times_vector_v
  use mo_derived,              only: spec_heat_cap_diagnostics_v
  use mo_bc_nml,               only: lperiodic
  
  implicit none
  
  contains
  
  subroutine hor_viscosity(state,diag)
  
    ! This subroutine computes the effective diffusion coefficient (molecular + turbulent) acting on horizontal divergent movements.
    
    ! input arguments and output
    type(t_state), intent(in)    :: state ! the state variables of the model atmosphere
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! computing the eddy viscosity
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx
        do jl=1,n_layers
          diag%viscosity_coeff_div(ji,jk,jl) = tke2hor_diff_coeff(diag%tke(ji,jk,jl),eff_hor_res)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! calculation of the molecular diffusion coefficient
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx
        do jl=1,n_layers
          diag%viscosity_molecular(ji,jk,jl) = calc_diffusion_coeff(diag%temperature(ji,jk,jl), &
          state%rho(ji,jk,jl,n_condensed_constituents+1))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! adding the molecular diffusion coefficient
    !$omp parallel workshare
    diag%viscosity_coeff_div = diag%viscosity_molecular + diag%viscosity_coeff_div
    !$omp end parallel workshare
    
    ! multiplying by the density
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx
        do jl=1,n_layers
          diag%viscosity_coeff_div(ji,jk,jl) = state%rho(ji,jk,jl,n_condensed_constituents+1) &
                                               *diag%viscosity_coeff_div(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! initialization with zeros
    !$omp parallel workshare
     diag%viscosity_coeff_curl_dual = 0._wp
    !$omp end parallel workshare
    
    ! molecular component
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do ji=2,ny
        do jk=2,nx
          diag%viscosity_coeff_curl_dual(ji,jk,jl) = 0.25_wp*sum(diag%viscosity_molecular(ji-1:ji,jk-1:jk,jl))
        enddo
      enddo
      
      ! periodic boundary conditions
      if (lperiodic) then
      
        do jk=2,nx
          diag%viscosity_coeff_curl_dual(1,jk,jl) = 0.25_wp*( &
          diag%viscosity_molecular(1,jk-1,jl) + diag%viscosity_molecular(ny,jk-1,jl) &
          + diag%viscosity_molecular(1,jk,jl) + diag%viscosity_molecular(ny,jk,jl))
          diag%viscosity_coeff_curl_dual(ny+1,jk,jl) = diag%viscosity_coeff_curl_dual(1,jk,jl)
        enddo
        do ji=2,ny
          diag%viscosity_coeff_curl_dual(ji,1,jl) = 0.25_wp*( &
          diag%viscosity_molecular(ji-1,nx,jl) + diag%viscosity_molecular(ji-1,1,jl) &
          + diag%viscosity_molecular(ji,nx,jl) + diag%viscosity_molecular(ji,1,jl))
          diag%viscosity_coeff_curl_dual(ji,nx+1,jl) = diag%viscosity_coeff_curl_dual(ji,1,jl)
        enddo
      
        ! corners
        diag%viscosity_coeff_curl_dual(1,1,jl) = 0.25*(diag%viscosity_molecular(1,1,jl) + diag%viscosity_molecular(ny,1,jl) &
        + diag%viscosity_molecular(1,nx,jl) + diag%viscosity_molecular(ny,nx,jl))
        diag%viscosity_coeff_curl_dual(1,nx+1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl)
        diag%viscosity_coeff_curl_dual(ny+1,1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl)
        diag%viscosity_coeff_curl_dual(ny+1,nx+1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl)
        
      endif
    
    enddo
    !$omp end parallel do
        
    ! turbulent component
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do ji=2,ny
        do jk=2,nx
          diag%viscosity_coeff_curl_dual(ji,jk,jl) = 0.25_wp*( &
          tke2hor_diff_coeff(diag%tke(ji-1,jk-1,jl),eff_hor_res) &
          + tke2hor_diff_coeff(diag%tke(ji-1,jk,jl),eff_hor_res) &
          + tke2hor_diff_coeff(diag%tke(ji,jk-1,jl),eff_hor_res) &
          + tke2hor_diff_coeff(diag%tke(ji,jk,jl),eff_hor_res))
        enddo
      enddo
      
      ! periodic boundary conditions
      if (lperiodic) then
      
        do jk=2,nx
          diag%viscosity_coeff_curl_dual(1,jk,jl) = 0.25_wp*( &
          tke2hor_diff_coeff(diag%tke(1,jk-1,jl),eff_hor_res) + tke2hor_diff_coeff(diag%tke(ny,jk-1,jl),eff_hor_res) &
          + tke2hor_diff_coeff(diag%tke(1,jk,jl),eff_hor_res) + tke2hor_diff_coeff(diag%tke(ny,jk,jl),eff_hor_res))
          diag%viscosity_coeff_curl_dual(ny+1,jk,jl) = diag%viscosity_coeff_curl_dual(1,jk,jl)
        enddo
        do ji=2,ny
          diag%viscosity_coeff_curl_dual(ji,1,jl) = 0.25_wp*( &
          tke2hor_diff_coeff(diag%tke(ji-1,nx,jl),eff_hor_res) + tke2hor_diff_coeff(diag%tke(ji-1,1,jl),eff_hor_res) &
          + tke2hor_diff_coeff(diag%tke(ji,nx,jl),eff_hor_res) + tke2hor_diff_coeff(diag%tke(ji,1,jl),eff_hor_res))
          diag%viscosity_coeff_curl_dual(ji,nx+1,jl) = diag%viscosity_coeff_curl_dual(ji,1,jl)
        enddo
      
        ! corners
        diag%viscosity_coeff_curl_dual(1,1,jl) = 0.25*( &
        tke2hor_diff_coeff(diag%tke(1,1,jl),eff_hor_res) + tke2hor_diff_coeff(diag%tke(ny,1,jl),eff_hor_res) &
        + tke2hor_diff_coeff(diag%tke(1,nx,jl),eff_hor_res) + tke2hor_diff_coeff(diag%tke(ny,nx,jl),eff_hor_res))
        diag%viscosity_coeff_curl_dual(1,nx+1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl)
        diag%viscosity_coeff_curl_dual(ny+1,1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl)
        diag%viscosity_coeff_curl_dual(ny+1,nx+1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl)
        
      endif
    
    enddo
    !$omp end parallel do
    
    ! multiplication by the density
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do ji=2,ny
        do jk=2,nx
          diag%viscosity_coeff_curl_dual(ji,jk,jl) = 0.25_wp*( &
          state%rho(ji-1,jk-1,jl,n_condensed_constituents+1) &
          + state%rho(ji-1,jk,jl,n_condensed_constituents+1) &
          + state%rho(ji,jk-1,jl,n_condensed_constituents+1) &
          + state%rho(ji,jk,jl,n_condensed_constituents+1)) &
          *diag%viscosity_coeff_curl_dual(ji,jk,jl)
        enddo
      enddo
      
      ! periodic boundary conditions
      if (lperiodic) then
      
        do jk=2,nx
          diag%viscosity_coeff_curl_dual(1,jk,jl) = diag%viscosity_coeff_curl_dual(1,jk,jl) &
          *(0.25_wp*( &
          state%rho(1,jk-1,jl,n_condensed_constituents+1) &
          + state%rho(ny,jk-1,jl,n_condensed_constituents+1) &
          + state%rho(1,jk,jl,n_condensed_constituents+1) &
          + state%rho(ny,jk,jl,n_condensed_constituents+1)))
          diag%viscosity_coeff_curl_dual(ny+1,jk,jl) = diag%viscosity_coeff_curl_dual(1,jk,jl)
        enddo
        do ji=2,ny
          diag%viscosity_coeff_curl_dual(ji,1,jl) = diag%viscosity_coeff_curl_dual(ji,1,jl) &
          *(0.25_wp*( &
          state%rho(ji-1,nx,jl,n_condensed_constituents+1) &
          + state%rho(ji-1,1,jl,n_condensed_constituents+1) &
          + state%rho(ji,nx,jl,n_condensed_constituents+1) &
          + state%rho(ji,1,jl,n_condensed_constituents+1)))
          diag%viscosity_coeff_curl_dual(ji,nx+1,jl) = diag%viscosity_coeff_curl_dual(ji,1,jl)
        enddo
      
        ! corners
        diag%viscosity_coeff_curl_dual(1,1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl) &
        *(0.25_wp*( &
        state%rho(ny,nx,jl,n_condensed_constituents+1) &
        + state%rho(ny,1,jl,n_condensed_constituents+1) &
        + state%rho(1,nx,jl,n_condensed_constituents+1) &
        + state%rho(1,1,jl,n_condensed_constituents+1)))
        diag%viscosity_coeff_curl_dual(1,nx+1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl)
        diag%viscosity_coeff_curl_dual(ny+1,1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl)
        diag%viscosity_coeff_curl_dual(ny+1,nx+1,jl) = diag%viscosity_coeff_curl_dual(1,1,jl)
        
      endif
    
    enddo
    !$omp end parallel do
    
    ! averaging the curl diffusion coefficient to the cell centers
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx
        do jl=1,n_layers
          diag%viscosity_coeff_curl(ji,jk,jl) = 0.25_wp*sum(diag%viscosity_coeff_curl_dual(ji:ji+1,jk:jk+1,jl))
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine hor_viscosity
  
  subroutine vert_hor_mom_viscosity(state,diag,grid)
  
    ! This subroutine computes the effective viscosity (eddy + molecular viscosity) for the vertical diffusion of horizontal velocity.
	! This quantity is located at the half level edges.
	! To obey the symmetry of the stress tensor, the same coefficient must be used for the horizontal diffusion of vertical velocity.
  
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl               ! loop indices
	
	!  updating the TKE
    call tke_update(state,diag,grid)
    
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jl=2,n_layers
        do jk=2,nx
          diag%vert_hor_viscosity_u(ji,jk,jl) = &
          ! molecular component
          0.25_wp*( &
          diag%viscosity_molecular(ji,jk-1,jl-1) + diag%viscosity_molecular(ji,jk,jl-1) &
          + diag%viscosity_molecular(ji,jk-1,jl) + diag%viscosity_molecular(ji,jk,jl)) &
          ! turbulent component
          + 0.25_wp*( &
          tke2vert_diff_coeff(diag%tke(ji,jk-1,jl-1),diag%n_squared(ji,jk-1,jl-1),grid%layer_thickness(ji,jk-1,jl-1)) + &
          tke2vert_diff_coeff(diag%tke(ji,jk,jl-1),diag%n_squared(ji,jk,jl-1),grid%layer_thickness(ji,jk,jl-1)) + &
          tke2vert_diff_coeff(diag%tke(ji,jk-1,jl),diag%n_squared(ji,jk-1,jl),grid%layer_thickness(ji,jk-1,jl)) + &
          tke2vert_diff_coeff(diag%tke(ji,jk,jl),diag%n_squared(ji,jk,jl),grid%layer_thickness(ji,jk,jl)))
          
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          diag%vert_hor_viscosity_u(ji,1,jl) = &
          ! molecular component
          0.25_wp*( &
          diag%viscosity_molecular(ji,nx,jl-1) + diag%viscosity_molecular(ji,1,jl-1) &
          + diag%viscosity_molecular(ji,nx,jl) + diag%viscosity_molecular(ji,1,jl)) &
          ! turbulent component
          + 0.25_wp*( &
          tke2vert_diff_coeff(diag%tke(ji,nx,jl-1),diag%n_squared(ji,nx,jl-1),grid%layer_thickness(ji,nx,jl-1)) + &
          tke2vert_diff_coeff(diag%tke(ji,1,jl-1),diag%n_squared(ji,1,jl-1),grid%layer_thickness(ji,1,jl-1)) + &
          tke2vert_diff_coeff(diag%tke(ji,nx,jl),diag%n_squared(ji,nx,jl),grid%layer_thickness(ji,nx,jl)) + &
          tke2vert_diff_coeff(diag%tke(ji,1,jl),diag%n_squared(ji,1,jl),grid%layer_thickness(ji,1,jl)))
          
          diag%vert_hor_viscosity_u(ji,nx+1,jl) = diag%vert_hor_viscosity_u(ji,1,jl)
          
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl)
    do jk=1,nx
      do jl=2,n_layers
        do ji=2,ny
          diag%vert_hor_viscosity_v(ji,jk,jl) = &
          ! molecular component
          0.25_wp*( &
          diag%viscosity_molecular(ji-1,jk,jl-1) + diag%viscosity_molecular(ji,jk,jl-1) &
          + diag%viscosity_molecular(ji-1,jk,jl) + diag%viscosity_molecular(ji,jk,jl)) &
          ! turbulent component
          + 0.25_wp*( &
          tke2vert_diff_coeff(diag%tke(ji-1,jk,jl-1),diag%n_squared(ji-1,jk,jl-1),grid%layer_thickness(ji-1,jk,jl-1)) + &
          tke2vert_diff_coeff(diag%tke(ji,jk,jl-1),diag%n_squared(ji,jk,jl-1),grid%layer_thickness(ji,jk,jl-1)) + &
          tke2vert_diff_coeff(diag%tke(ji-1,jk,jl),diag%n_squared(ji-1,jk,jl),grid%layer_thickness(ji-1,jk,jl)) + &
          tke2vert_diff_coeff(diag%tke(ji,jk,jl),diag%n_squared(ji,jk,jl),grid%layer_thickness(ji,jk,jl)))
          
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          diag%vert_hor_viscosity_v(1,jk,jl) = &
          ! molecular component
          0.25_wp*( &
          diag%viscosity_molecular(ny,jk,jl-1) + diag%viscosity_molecular(1,jk,jl-1) &
          + diag%viscosity_molecular(ny,jk,jl) + diag%viscosity_molecular(1,jk,jl)) &
          ! turbulent component
          + 0.25_wp*( &
          tke2vert_diff_coeff(diag%tke(ny,jk,jl-1),diag%n_squared(ny,jk,jl-1),grid%layer_thickness(ny,jk,jl-1)) + &
          tke2vert_diff_coeff(diag%tke(1,jk,jl-1),diag%n_squared(1,jk,jl-1),grid%layer_thickness(1,jk,jl-1)) + &
          tke2vert_diff_coeff(diag%tke(ny,jk,jl),diag%n_squared(ny,jk,jl),grid%layer_thickness(ny,jk,jl)) + &
          tke2vert_diff_coeff(diag%tke(1,jk,jl),diag%n_squared(1,jk,jl),grid%layer_thickness(1,jk,jl)))
          
          diag%vert_hor_viscosity_v(ny+1,jk,jl) = diag%vert_hor_viscosity_v(1,jk,jl)
          
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    ! multiplication by the density
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jl=2,n_layers
        do jk=2,nx
          diag%vert_hor_viscosity_u(ji,jk,jl) = diag%vert_hor_viscosity_u(ji,jk,jl) &
          ! molecular component
          *0.25_wp*( &
          state%rho(ji,jk-1,jl-1,n_condensed_constituents+1) &
          + state%rho(ji,jk,jl-1,n_condensed_constituents+1) &
          + state%rho(ji,jk-1,jl,n_condensed_constituents+1) &
          + state%rho(ji,jk,jl,n_condensed_constituents+1))
          
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          diag%vert_hor_viscosity_u(ji,1,jl) = diag%vert_hor_viscosity_u(ji,1,jl) &
          ! molecular component
          *0.25_wp*( &
          state%rho(ji,nx,jl-1,n_condensed_constituents+1) &
          + state%rho(ji,1,jl-1,n_condensed_constituents+1) &
          + state%rho(ji,nx,jl,n_condensed_constituents+1) &
          + state%rho(ji,1,jl,n_condensed_constituents+1))
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl)
    do jk=1,nx
      do jl=2,n_layers
        do ji=2,ny
          diag%vert_hor_viscosity_v(ji,jk,jl) = diag%vert_hor_viscosity_v(ji,jk,jl) &
          ! molecular component
          *0.25_wp*( &
          state%rho(ji-1,jk,jl-1,n_condensed_constituents+1) &
          + state%rho(ji,jk,jl-1,n_condensed_constituents+1) &
          + state%rho(ji-1,jk,jl,n_condensed_constituents+1) &
          + state%rho(ji,jk,jl,n_condensed_constituents+1))
          
        enddo
        
        ! periodic boundary conditions
        if (lperiodic) then
          diag%vert_hor_viscosity_v(1,jk,jl) = diag%vert_hor_viscosity_v(1,jk,jl) &
          ! molecular component
          *0.25_wp*( &
          state%rho(ny,jk,jl-1,n_condensed_constituents+1) &
          + state%rho(1,jk,jl-1,n_condensed_constituents+1) &
          + state%rho(ny,jk,jl,n_condensed_constituents+1) &
          + state%rho(1,jk,jl,n_condensed_constituents+1))
          
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel workshare
    ! for now, we set the vertical diffusion coefficient at the TOA equal to the vertical diffusion coefficient in the layer below
    diag%vert_hor_viscosity_u(:,:,1) = diag%vert_hor_viscosity_u(:,:,2)
    ! for now, we set the vertical diffusion coefficient at the surface equal to the vertical diffusion coefficient in the layer above
    diag%vert_hor_viscosity_u(:,:,n_levels) = diag%vert_hor_viscosity_u(:,:,n_layers)
    !$omp end parallel workshare
    !$omp parallel workshare
    ! for now, we set the vertical diffusion coefficient at the TOA equal to the vertical diffusion coefficient in the layer below
    diag%vert_hor_viscosity_v(:,:,1) = diag%vert_hor_viscosity_v(:,:,2)
    ! for now, we set the vertical diffusion coefficient at the surface equal to the vertical diffusion coefficient in the layer above
    diag%vert_hor_viscosity_v(:,:,n_levels) = diag%vert_hor_viscosity_v(:,:,n_layers)
    !$omp end parallel workshare
  
  end subroutine vert_hor_mom_viscosity
  
  subroutine vert_vert_mom_viscosity(state,diag,grid)
  
    ! This subroutine multiplies scalar_placeholder (containing dw/dz) by the diffusion coefficient acting on w because of w.
   
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    real(wp) :: mom_diff_coeff        ! the diffusion coefficient
    integer  :: ji,jk,jl              ! loop indices
    
    !$omp parallel do private(ji,jk,jl,mom_diff_coeff)
    do ji=1,ny
      do jk=1,nx
        do jl=1,n_layers
    
          mom_diff_coeff &
          ! molecular viscosity
          = diag%viscosity_molecular(ji,jk,jl) &
          ! turbulent component
          + tke2vert_diff_coeff(diag%tke(ji,jk,jl),diag%n_squared(ji,jk,jl),grid%layer_thickness(ji,jk,jl))

          diag%scalar_placeholder(ji,jk,jl) = state%rho(ji,jk,jl,n_condensed_constituents+1) &
                                              *mom_diff_coeff*diag%scalar_placeholder(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine vert_vert_mom_viscosity
  
  subroutine scalar_diffusion_coeffs(state,diag,grid)
  
    ! This subroutine computes the scalar diffusion coefficients (including eddies).
  
    ! input arguments and output
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji,jk,jl ! loop variables
    real(wp) :: c_g_v    ! specific heat capacity
    
    ! The eddy viscosity coefficient and the TKE only has to be calculated if it has not yet been done.
    if (.not. lmom_diff_h .and. .not. ltemp_diff_h) then
    
      call hor_viscosity(state,diag)
      call tke_update(state,diag,grid)
      
      ! molecular viscosity
      !$omp parallel do private(ji,jk,jl)
      do ji=1,ny
        do jk=1,nx
          do jl=1,n_layers
            diag%viscosity_molecular(ji,jk,jl) = calc_diffusion_coeff(diag%temperature(ji,jk,jl), &
            state%rho(ji,jk,jl,n_condensed_constituents+1))
          enddo
        enddo
      enddo
      !$omp end parallel do
    
    endif
    
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx
        do jl=1,n_layers
          ! horizontal diffusion coefficient
          diag%mass_diffusion_coeff_numerical_h(ji,jk,jl) &
          = 0.5_wp*(diag%viscosity_coeff_div(ji,jk,jl) + diag%viscosity_coeff_curl(ji,jk,jl)) &
          /state%rho(ji,jk,jl,n_condensed_constituents+1)
          ! vertical diffusion coefficient
          diag%mass_diffusion_coeff_numerical_v(ji,jk,jl) &
          ! molecular component
          = diag%viscosity_molecular(ji,jk,jl) &
          ! turbulent component
          + tke2vert_diff_coeff(diag%tke(ji,jk,jl),diag%n_squared(ji,jk,jl),grid%layer_thickness(ji,jk,jl))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! The eddy viscosity coefficient and the TKE only has to be calculated if it has not yet been done.
    if (.not. lmom_diff_h) then
    
      call hor_viscosity(state,diag)
      call tke_update(state,diag,grid)
      
      ! molecular viscosity
      !$omp parallel do private(ji,jk,jl)
      do ji=1,ny
        do jk=1,nx
          do jl=1,n_layers
            diag%viscosity_molecular(ji,jk,jl) = calc_diffusion_coeff(diag%temperature(ji,jk,jl), &
            state%rho(ji,jk,jl,n_condensed_constituents+1))
          enddo
        enddo
      enddo
      !$omp end parallel do
    
    endif
    
    !$omp parallel do private(ji,jk,jl,c_g_v)
    do ji=1,ny
      do jk=1,nx
        do jl=1,n_layers
          c_g_v = spec_heat_cap_diagnostics_v(state,ji,jk,jl)
          ! horizontal diffusion coefficient
          diag%temp_diffusion_coeff_numerical_h(ji,jk,jl) = c_g_v &
          *0.5_wp*(diag%viscosity_coeff_div(ji,jk,jl) + diag%viscosity_coeff_curl(ji,jk,jl))
          ! vertical diffusion coefficient
          diag%temp_diffusion_coeff_numerical_v(ji,jk,jl) &
          ! molecular component
          = state%rho(ji,jk,jl,n_condensed_constituents+1) &
          *c_g_v*(diag%viscosity_molecular(ji,jk,jl) &
          ! turbulent component
          + tke2vert_diff_coeff(diag%tke(ji,jk,jl),diag%n_squared(ji,jk,jl),grid%layer_thickness(ji,jk,jl)))
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine scalar_diffusion_coeffs
  
  subroutine update_n_squared(state,diag,grid)
    
    ! This subroutine calculates the Brunt-Väisälä frequency.
    
    type(t_state), intent(in)    :: state ! state which to use for the calculation
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer :: ji,jk,jl
    
    ! calculating the full virtual potential temperature
    !$omp parallel workshare
    diag%scalar_placeholder = grid%theta_v_bg+state%theta_v_pert
    !$omp end parallel workshare
    ! vertical gradient of the full virtual potential temperature
    call grad_vert(diag%scalar_placeholder,diag%w_placeholder,grid)
    ! calculating the inverse full virtual potential temperature
    !$omp parallel workshare
    diag%scalar_placeholder = 1.0/diag%scalar_placeholder
    !$omp end parallel workshare
    call scalar_times_vector_v(diag%scalar_placeholder,diag%w_placeholder,diag%w_placeholder)
    
    ! multiplying by the gravity acceleration
    call scalar_times_vector_v(diag%w_placeholder,grid%gravity_m_v,diag%w_placeholder)
    
    ! averaging vertically to the scalar points
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx
        do jl=1,n_layers
          if (jl==1) then
            diag%n_squared(ji,jk,jl) = diag%w_placeholder(ji,jk,jl)
          elseif (jl==n_layers) then
            diag%n_squared(ji,jk,jl) = diag%w_placeholder(ji,jk,jl+1)
          else
            diag%n_squared(ji,jk,jl) &
            = grid%inner_product_weights(ji,jk,jl,5)*diag%w_placeholder(ji,jk,jl) &
            + grid%inner_product_weights(ji,jk,jl,6)*diag%w_placeholder(ji,jk,jl+1)
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine update_n_squared
  
  function tke2hor_diff_coeff(tke,effective_resolution)
  
    ! This function returns the horizontal kinematic eddy viscosity as a function of the specific TKE.
    
    real(wp), intent(in)  :: tke,effective_resolution
    real(wp)              :: tke2hor_diff_coeff
    
    ! local variables
    real(wp) :: mean_velocity,mean_free_path
    
    mean_velocity = (2._wp*tke)**0.5_wp
    mean_free_path = effective_resolution/6._wp
    tke2hor_diff_coeff = 1._wp/6._wp*mean_free_path*mean_velocity
  
  end function tke2hor_diff_coeff

  function tke2vert_diff_coeff(tke,n_squared,layer_thickness)

    ! This function returns the vertical kinematic eddy viscosity as a function of the specific TKE and the Brunt-Väisälä frequency.
    
    real(wp), intent(in)  :: tke,n_squared,layer_thickness
    real(wp)              :: tke2vert_diff_coeff
    
    ! local variables
    real(wp) :: tke_vert,mean_velocity,n_used,mean_free_path
  
    ! vertical component of the turbulent kinetic energy
    tke_vert = 3._wp*1e-3_wp*tke
  
    mean_velocity = (2._wp*tke_vert)**0.5_wp
    ! used Brunt-Väisälä frequency
    n_used = (max(n_squared,1e-4_wp))**0.5_wp
    mean_free_path = (2._wp*tke_vert)**0.5_wp/n_used
    mean_free_path = min(mean_free_path,layer_thickness)
    tke2vert_diff_coeff = 1._wp/6._wp*mean_free_path*mean_velocity
    
  end function tke2vert_diff_coeff
  
end module mo_eff_diff_coeffs







