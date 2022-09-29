! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_vorticities

  ! This module contains the calculation of the vorticities.

  use mo_definitions,      only: t_state,t_diag,t_grid,wp
  use mo_run_nml,          only: ny,nx,n_layers,n_levels,n_oro_layers,toa,lcorio,llinear,n_flat_layers
  use mo_constants,        only: r_e
  use mo_constituents_nml, only: n_condensed_constituents
  use mo_averaging,        only: horizontal_covariant_x,horizontal_covariant_y
  use mo_bc_nml,           only: lperiodic
  
  implicit none
  
  contains

  subroutine rel_vort(state,diag,grid)
  
    ! This subroutine calculates the relative vorticity.
    
    type(t_state), intent(in)    :: state ! state to work with
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! model grid
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! calculating the relative vorticity in x-direction
    !$omp parallel do private(ji,jk,jl)
    do jk=1,nx
      do jl=2,n_layers
        do ji=2,ny
          diag%zeta_x(ji,jk,jl) = &
          grid%dz(ji-1,jk,jl)*state%wind_w(ji-1,jk,jl) &
          - grid%dy(ji,jk,jl-1)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,ji,jk,jl-1) &
          - grid%dz(ji,jk,jl)*state%wind_w(ji,jk,jl) &
          + grid%dy(ji,jk,jl)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,ji,jk,jl)
        enddo
        
        ! boundary conditions
        if (lperiodic) then
          diag%zeta_x(1,jk,jl) = &
          grid%dz(ny,jk,jl)*state%wind_w(ny,jk,jl) &
          - grid%dy(1,jk,jl-1)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,1,jk,jl-1) &
          - grid%dz(1,jk,jl)*state%wind_w(1,jk,jl) &
          + grid%dy(1,jk,jl)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,1,jk,jl)
          diag%zeta_x(ny+1,jk,jl) = diag%zeta_x(1,jk,jl)
        else
          diag%zeta_x(1,jk,jl) = 0._wp
          diag%zeta_x(ny+1,jk,jl) = 0._wp
        endif
        
        do ji=1,ny+1
          ! At the surface, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
          diag%zeta_x(ji,jk,n_levels) = -grid%dy(ji,jk,n_layers) &
                                        *horizontal_covariant_y(state%wind_v,state%wind_w,grid,ji,jk,n_layers)
        enddo
      enddo
    enddo
    !$omp end parallel do
    ! At the TOA, the horizontal vorticity is assumed to have no vertical shear.
    !$omp parallel workshare
    diag%zeta_x(:,:,1) = diag%zeta_x(:,:,2)
    !$omp end parallel workshare
    ! dividing by the area
    !$omp parallel workshare
    diag%zeta_x = diag%zeta_x/grid%area_dual_x
    !$omp end parallel workshare
    
    ! calculating the relative vorticity in y-direction
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jl=2,n_layers
        do jk=2,nx
          diag%zeta_y(ji,jk,jl) = &
          -grid%dz(ji,jk,jl)*state%wind_w(ji,jk,jl) &
          - grid%dx(ji,jk,jl)*horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,jk,jl) &
          + grid%dz(ji,jk-1,jl)*state%wind_w(ji,jk-1,jl) &
          + grid%dx(ji,jk,jl-1)*horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,jk,jl-1)
        enddo
        
        ! boundaries
        if (lperiodic) then
          diag%zeta_y(ji,1,jl) = &
          -grid%dz(ji,1,jl)*state%wind_w(ji,1,jl) &
          - grid%dx(ji,1,jl)*horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,1,jl) &
          + grid%dz(ji,nx,jl)*state%wind_w(ji,nx,jl) &
          + grid%dx(ji,1,jl-1)*horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,1,jl-1)
          diag%zeta_y(ji,nx+1,jl) = diag%zeta_y(ji,1,jl)
        else
          diag%zeta_y(ji,1,jl) = 0._wp
          diag%zeta_y(ji,nx+1,jl) = 0._wp
        endif
        
        do jk=1,nx+1
          ! At the surface, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
          diag%zeta_y(ji,jk,n_levels) = grid%dx(ji,jk,n_layers) &
                                        *horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,jk,n_layers)
        enddo
        
      enddo
    enddo
    !$omp end parallel do
    ! At the TOA, the horizontal vorticity is assumed to have no vertical shear.
    !$omp parallel workshare
    diag%zeta_y(:,:,1) = diag%zeta_y(:,:,2)
    !$omp end parallel workshare
    ! dividing by the area
    !$omp parallel workshare
    diag%zeta_y = diag%zeta_y/grid%area_dual_y
    !$omp end parallel workshare
      
    ! calculating the relative vorticity in z-direction
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny+1
      do jk=1,nx+1
        do jl=1,n_layers
          diag%zeta_z(ji,jk,jl) = rel_vort_z_local(state,grid,ji,jk,jl)/grid%area_dual_z(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
      
  end subroutine rel_vort
  
  subroutine calc_pot_vort(state,diag,grid)
  
    ! This subroutine calculates the potential vorticity.
    
    type(t_state), intent(in)    :: state ! state to work with
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! model grid
  
    ! local variables
    integer :: ji,jk,jl ! loop indices
  
    ! calculating the relative vorticity
    if (.not. llinear) then
      call rel_vort(state,diag,grid)
    else
      !$omp parallel workshare
      diag%zeta_x = 0._wp
      diag%zeta_y = 0._wp
      diag%zeta_z = 0._wp
      !$omp end parallel workshare
    endif
    
    ! adding the Coriolis vector to the relative vorticity to obtain the absolute vorticity
    if (lcorio) then
      !$omp parallel do private(jl)
      do jl=1,n_levels
        diag%eta_x(:,:,jl) = diag%zeta_x(:,:,jl) + grid%fvec_x(:,:)
        diag%eta_y(:,:,jl) = diag%zeta_y(:,:,jl) + grid%fvec_y(:,:)
      enddo
      !$omp end parallel do
      !$omp parallel do private(jl)
      do jl=1,n_layers
        diag%eta_z(:,:,jl) = diag%zeta_z(:,:,jl) + grid%fvec_z(:,:)
      enddo
      !$omp end parallel do
    else
      !$omp parallel do private(jl)
      do jl=1,n_levels
        diag%eta_x(:,:,jl) = diag%zeta_x(:,:,jl)
        diag%eta_y(:,:,jl) = diag%zeta_y(:,:,jl)
      enddo
      !$omp end parallel do
      !$omp parallel do private(jl)
      do jl=1,n_layers
        diag%eta_z(:,:,jl) = diag%zeta_z(:,:,jl)
      enddo
      !$omp end parallel do
    endif
    
    ! dividing by the averaged density to obtain the "potential vorticity"
    ! horizontal vorticity in x-direction
    !$omp parallel do private(ji,jk,jl)
    do jk=1,nx
      do jl=1,n_levels
        do ji=2,ny
          if (jl==1) then
            diag%eta_x(ji,jk,jl) = diag%eta_x(ji,jk,jl)/(0.5_wp*(state%rho(ji-1,jk,jl,n_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji-1,jk,jl))*(state%rho(ji-1,jk,jl,n_condensed_constituents+1) &
            -state%rho(ji-1,jk,jl+1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji-1,jk,jl)-grid%z_scalar(ji-1,jk,jl+1)) &
            +state%rho(ji,jk,jl,n_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,jk,jl))*(state%rho(ji,jk,jl,n_condensed_constituents+1) &
            -state%rho(ji,jk,jl+1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))))
          elseif (jl==n_levels) then
            diag%eta_x(ji,jk,jl) = diag%eta_x(ji,jk,jl)/(0.5_wp*(state%rho(ji-1,jk,jl-1,n_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji-1,jk,jl)-grid%z_scalar(ji-1,jk,jl-1))*(state%rho(ji-1,jk,jl-2,n_condensed_constituents+1) &
            -state%rho(ji-1,jk,jl-1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji-1,jk,jl-2)-grid%z_scalar(ji-1,jk,jl-1)) &
            +state%rho(ji,jk,jl-1,n_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,jk,jl)-grid%z_scalar(ji,jk,jl-1)) &
            *(state%rho(ji,jk,jl-2,n_condensed_constituents+1) &
            -state%rho(ji,jk,jl-1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk,jl-2)-grid%z_scalar(ji,jk,jl-1))))
          else
            diag%eta_x(ji,jk,jl) = diag%eta_x(ji,jk,jl)/(0.25_wp*(state%rho(ji-1,jk,jl-1,n_condensed_constituents+1) &
            +state%rho(ji,jk,jl-1,n_condensed_constituents+1)+ &
            state%rho(ji-1,jk,jl,n_condensed_constituents+1)+state%rho(ji,jk,jl,n_condensed_constituents+1)))
          endif
        enddo
        
        if (lperiodic) then
          if (jl==1) then
            diag%eta_x(1,jk,jl) = diag%eta_x(1,jk,jl)/(0.5_wp*(state%rho(ny,jk,jl,n_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ny,jk,jl))*(state%rho(ny,jk,jl,n_condensed_constituents+1) &
            -state%rho(ny,jk,jl+1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ny,jk,jl)-grid%z_scalar(ny,jk,jl+1)) &
            +state%rho(1,jk,jl,n_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(1,jk,jl))*(state%rho(1,jk,jl,n_condensed_constituents+1) &
            -state%rho(1,jk,jl+1,n_condensed_constituents+1))/ &
            (grid%z_scalar(1,jk,jl)-grid%z_scalar(1,jk,jl+1))))
          elseif (jl==n_levels) then
            diag%eta_x(1,jk,jl) = diag%eta_x(1,jk,jl)/(0.5_wp*(state%rho(ny,jk,jl-1,n_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ny,jk,jl)-grid%z_scalar(ny,jk,jl-1))*(state%rho(ny,jk,jl-2,n_condensed_constituents+1) &
            -state%rho(ny,jk,jl-1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ny,jk,jl-2)-grid%z_scalar(ny,jk,jl-1)) &
            +state%rho(1,jk,jl-1,n_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(1,jk,jl)-grid%z_scalar(1,jk,jl-1)) &
            *(state%rho(1,jk,jl-2,n_condensed_constituents+1) &
            -state%rho(1,jk,jl-1,n_condensed_constituents+1))/ &
            (grid%z_scalar(1,jk,jl-2)-grid%z_scalar(1,jk,jl-1))))
          else
            diag%eta_x(1,jk,jl) = diag%eta_x(1,jk,jl)/(0.25_wp*(state%rho(ny,jk,jl-1,n_condensed_constituents+1) &
            +state%rho(1,jk,jl-1,n_condensed_constituents+1)+ &
            state%rho(ny,jk,jl,n_condensed_constituents+1)+state%rho(1,jk,jl,n_condensed_constituents+1)))
          endif
          diag%eta_x(ny+1,jk,jl) = diag%eta_x(1,jk,jl)
        endif
        
      enddo
    enddo
    !$omp end parallel do
  
    ! horizontal vorticity in y-direction
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jl=1,n_levels
        do jk=2,nx
          if (jl==1) then
            diag%eta_y(ji,jk,jl) = diag%eta_y(ji,jk,jl)/(0.5_wp*(state%rho(ji,jk-1,jl,n_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,jk-1,jl))*(state%rho(ji,jk-1,jl,n_condensed_constituents+1) &
            -state%rho(ji,jk-1,jl+1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk-1,jl)-grid%z_scalar(ji,jk-1,jl+1)) &
            +state%rho(ji,jk,jl,n_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,jk,jl))*(state%rho(ji,jk,jl,n_condensed_constituents+1) &
            -state%rho(ji,jk,jl+1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))))
          elseif (jl==n_levels) then
            diag%eta_y(ji,jk,jl) = diag%eta_y(ji,jk,jl)/(0.5_wp*(state%rho(ji,jk-1,jl-1,n_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,jk-1,jl)-grid%z_scalar(ji,jk-1,jl-1))*(state%rho(ji,jk-1,jl-2,n_condensed_constituents+1) &
            -state%rho(ji,jk-1,jl-1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk-1,jl-2)-grid%z_scalar(ji,jk-1,jl-1)) &
            +state%rho(ji,jk,jl-1,n_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,jk,jl)-grid%z_scalar(ji,jk,jl-1)) &
            *(state%rho(ji,jk,jl-2,n_condensed_constituents+1) &
            -state%rho(ji,jk,jl-1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk,jl-2)-grid%z_scalar(ji,jk,jl-1))))
          else
            diag%eta_y(ji,jk,jl) = diag%eta_y(ji,jk,jl)/(0.25_wp*(state%rho(ji,jk-1,jl-1,n_condensed_constituents+1) &
            +state%rho(ji,jk,jl-1,n_condensed_constituents+1)+ &
            state%rho(ji,jk-1,jl,n_condensed_constituents+1)+state%rho(ji,jk,jl,n_condensed_constituents+1)))
          endif
        enddo
        
        if (lperiodic) then
          if (jl==1) then
            diag%eta_y(ji,1,jl) = diag%eta_y(ji,1,jl)/(0.5_wp*(state%rho(ji,nx,jl,n_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,nx,jl))*(state%rho(ji,nx,jl,n_condensed_constituents+1) &
            -state%rho(ji,nx,jl+1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,nx,jl)-grid%z_scalar(ji,nx,jl+1)) &
            +state%rho(ji,1,jl,n_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,1,jl))*(state%rho(ji,1,jl,n_condensed_constituents+1) &
            -state%rho(ji,1,jl+1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,1,jl)-grid%z_scalar(ji,1,jl+1))))
          elseif (jl==n_levels) then
            diag%eta_y(ji,1,jl) = diag%eta_y(ji,1,jl)/(0.5_wp*(state%rho(ji,nx,jl-1,n_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,nx,jl)-grid%z_scalar(ji,nx,jl-1))*(state%rho(ji,nx,jl-2,n_condensed_constituents+1) &
            -state%rho(ji,nx,jl-1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,nx,jl-2)-grid%z_scalar(ji,nx,jl-1)) &
            +state%rho(ji,1,jl-1,n_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,1,jl)-grid%z_scalar(ji,1,jl-1)) &
            *(state%rho(ji,1,jl-2,n_condensed_constituents+1) &
            -state%rho(ji,1,jl-1,n_condensed_constituents+1))/ &
            (grid%z_scalar(ji,1,jl-2)-grid%z_scalar(ji,1,jl-1))))
          else
            diag%eta_y(ji,1,jl) = diag%eta_y(ji,1,jl)/(0.25_wp*(state%rho(ji,nx,jl-1,n_condensed_constituents+1) &
            +state%rho(ji,1,jl-1,n_condensed_constituents+1)+ &
            state%rho(ji,nx,jl,n_condensed_constituents+1)+state%rho(ji,1,jl,n_condensed_constituents+1)))
          endif
          diag%eta_y(ji,nx+1,jl) = diag%eta_y(ji,1,jl)
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    ! vertical vorticity
    !$omp parallel do private(ji,jk)
    do ji=2,ny
      do jk=2,nx
        diag%eta_z(ji,jk,:) = diag%eta_z(ji,jk,:)/(0.25_wp*(state%rho(ji,jk-1,:,n_condensed_constituents+1) &
        +state%rho(ji,jk,:,n_condensed_constituents+1) &
        +state%rho(ji-1,jk,:,n_condensed_constituents+1)+state%rho(ji-1,jk-1,:,n_condensed_constituents+1)))
      enddo
      
      ! periodic boundary conditions
      if (lperiodic) then
        diag%eta_z(ji,1,:) = diag%eta_z(ji,1,:)/(0.25_wp*(state%rho(ji,nx,:,n_condensed_constituents+1) &
        +state%rho(ji,1,:,n_condensed_constituents+1) &
        +state%rho(ji-1,1,:,n_condensed_constituents+1)+state%rho(ji-1,nx,:,n_condensed_constituents+1)))
        diag%eta_z(ji,nx+1,:) = diag%eta_z(ji,1,:)
      endif
      
    enddo
    !$omp end parallel do
    
    ! periodic boundary conditions
    if (lperiodic) then
      !$omp parallel do private(jk)
      do jk=2,nx
        diag%eta_z(1,jk,:) = diag%eta_z(1,jk,:)/(0.25_wp*(state%rho(1,jk-1,:,n_condensed_constituents+1) &
        +state%rho(1,jk,:,n_condensed_constituents+1) &
        +state%rho(ny,jk,:,n_condensed_constituents+1)+state%rho(ny,jk-1,:,n_condensed_constituents+1)))
        diag%eta_z(ny+1,jk,:) = diag%eta_z(1,jk,:)
      enddo
      !$omp end parallel do
    endif
    
    ! corners under periodic boundary conditions
    if (lperiodic) then
      !$omp parallel workshare
      diag%eta_z(1,1,:) = diag%eta_z(1,1,:)/(0.25_wp*(state%rho(1,nx,:,n_condensed_constituents+1) &
      +state%rho(1,1,:,n_condensed_constituents+1) &
      +state%rho(ny,1,:,n_condensed_constituents+1)+state%rho(ny,nx,:,n_condensed_constituents+1)))
      diag%eta_z(1,nx+1,:) = diag%eta_z(1,1,:)
      diag%eta_z(ny+1,1,:) = diag%eta_z(1,1,:)
      diag%eta_z(ny+1,nx+1,:) = diag%eta_z(1,1,:)
      !$omp end parallel workshare
    endif
    
  end subroutine calc_pot_vort
  
  function rel_vort_z_local(state,grid,ji,jk,jl)
  
    ! This function returns the vertical relative vorticity at a grindpoint.
    
    ! input arguments
    type(t_state) :: state    ! state with which to calculate the relative vorticity
    type(t_grid)  :: grid     ! grid properties
    integer       :: ji,jk,jl ! indices of the gridpoint
    ! result
    real(wp)      :: rel_vort_z_local
    
    ! local variables
    real(wp) :: delta_z,l_rescale,vertical_gradient       ! needed for terrain handling
    integer  :: ind_shift,j_i(4),j_k(4),jm,sign_vector(4) ! helper variables containing indices
    
    ! setting the indices
    j_i(1) = ji
    j_k(1) = jk
    j_i(2) = ji-1
    j_k(2) = jk
    j_i(3) = ji
    j_k(3) = jk-1
    j_i(4) = ji
    j_k(4) = jk
    sign_vector(1) = 1
    sign_vector(2) = -1
    sign_vector(3) = -1
    sign_vector(4) = 1
    
    ! initializing the result with zero
    rel_vort_z_local = 0._wp
    
    ! boundary handling
    if (ji==1 .or. jk==1 .or. ji==ny+1 .or. jk==nx+1) then
      if (lperiodic) then
        j_i(1) = ji
        if (jk==nx+1) then
          j_k(1) = 1
        else
          j_k(1) = jk
        endif
        
        if (ji==1) then
          j_i(2) = ny
        else
          j_i(2) = ji-1
        endif
        j_k(2) = jk
        
        j_i(3) = ji
        if (jk==1) then
          j_k(3) = nx
        else
          j_k(3) = jk-1
        endif
        
        if (ji==ny+1) then
          j_i(4) = 1
        else
          j_i(4) = ji
        endif
        j_k(4) = jk
      else
        return
      endif
    endif
    
    ! flat layers
    if (jl<=n_flat_layers) then
      rel_vort_z_local = rel_vort_z_local &
      + sign_vector(1)*grid%dy(j_i(1),j_k(1),jl)*state%wind_v(j_i(1),j_k(1),jl) &
      + sign_vector(2)*grid%dx(j_i(2),j_k(2),jl)*state%wind_u(j_i(2),j_k(2),jl) &
      + sign_vector(3)*grid%dy(j_i(3),j_k(3),jl)*state%wind_v(j_i(3),j_k(3),jl) &
      + sign_vector(4)*grid%dx(j_i(4),j_k(4),jl)*state%wind_u(j_i(4),j_k(4),jl)
    ! layers which follow the orography
    else
      do jm=1,4
        if (jm==1 .or. jm==3) then
          l_rescale = (r_e + grid%z_area_dual_z(ji,jk,jl))/(r_e + grid%z_v(j_i(jm),j_k(jm),jl))
          delta_z = grid%z_area_dual_z(ji,jk,jl) - grid%z_v(j_i(jm),j_k(jm),jl)
          ind_shift = 1
          if (delta_z>0._wp .or. jl==n_layers) then
            ind_shift = -1
          endif
          if (jl==1) then
            ind_shift = 1
          endif
          vertical_gradient = (state%wind_v(j_i(jm),j_k(jm),jl) - state%wind_v(j_i(jm),j_k(jm),jl+ind_shift))/ &
          (grid%z_v(j_i(jm),j_k(jm),jl) - grid%z_v(j_i(jm),j_k(jm),jl+ind_shift))
          rel_vort_z_local = rel_vort_z_local + sign_vector(jm)*l_rescale*grid%dy(j_i(jm),j_k(jm),jl)* &
          (state%wind_v(j_i(jm),j_k(jm),jl) + delta_z*vertical_gradient)
        else
          l_rescale = (r_e + grid%z_area_dual_z(ji,jk,jl))/(r_e + grid%z_u(j_i(jm),j_k(jm),jl))
          delta_z = grid%z_area_dual_z(ji,jk,jl) - grid%z_u(j_i(jm),j_k(jm),jl)
          ind_shift = 1
          if (delta_z>0._wp .or. jl==n_layers) then
            ind_shift = -1
          endif
          if (jl==1) then
            ind_shift = 1
          endif
          vertical_gradient = (state%wind_u(j_i(jm),j_k(jm),jl) - state%wind_u(j_i(jm),j_k(jm),jl+ind_shift))/ &
          (grid%z_u(j_i(jm),j_k(jm),jl) - grid%z_u(j_i(jm),j_k(jm),jl+ind_shift))
          rel_vort_z_local = rel_vort_z_local + sign_vector(jm)*l_rescale*grid%dx(j_i(jm),j_k(jm),jl)* &
          (state%wind_u(j_i(jm),j_k(jm),jl) + delta_z*vertical_gradient)
        endif
      enddo
    endif
  
  end function rel_vort_z_local

end module mo_vorticities








