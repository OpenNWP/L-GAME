! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module vorticities

  ! This module contains the calculation of the vorticities.

  use definitions,      only: t_state,t_diag,t_grid,wp
  use run_nml,          only: nlins,ncols,nlays,nlays_oro,toa,lcorio,llinear
  use constants,        only: re
  use constituents_nml, only: no_of_condensed_constituents
  use averaging,        only: horizontal_covariant_x,horizontal_covariant_y
  use bc_nml,           only: lperiodic
  
  implicit none
  
  private
  
  public :: rel_vort
  public :: calc_pot_vort
  
  contains

  subroutine rel_vort(state,diag,grid)
  
    ! This subroutine calculates the relative vorticity.
    
    type(t_state), intent(in)    :: state ! state to work with
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! model grid
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! calculating the relative vorticity in x-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do jk=1,ncols
      do jl=2,nlays
        do ji=2,nlins
          diag%zeta_x(ji,jk,jl) = &
          grid%dz(ji-1,jk,jl)*state%wind_w(ji-1,jk,jl) &
          - grid%dy(ji,jk,jl-1)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,ji,jk,jl-1) &
          - grid%dz(ji,jk,jl)*state%wind_w(ji,jk,jl) &
          + grid%dy(ji,jk,jl)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,ji,jk,jl)
        enddo
        
        ! boundary conditions
        if (lperiodic) then
          diag%zeta_x(1,jk,jl) = &
          grid%dz(nlins,jk,jl)*state%wind_w(nlins,jk,jl) &
          - grid%dy(1,jk,jl-1)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,1,jk,jl-1) &
          - grid%dz(1,jk,jl)*state%wind_w(1,jk,jl) &
          + grid%dy(1,jk,jl)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,1,jk,jl)
          diag%zeta_x(nlins+1,jk,jl) = diag%zeta_x(1,jk,jl)
        else
          diag%zeta_x(1,jk,jl) = 0._wp
          diag%zeta_x(nlins+1,jk,jl) = 0._wp
        endif
        
        do ji=1,nlins+1
          ! At the surface, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
          diag%zeta_x(ji,jk,nlays+1) = -grid%dy(ji,jk,nlays)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,ji,jk,nlays)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! At the TOA, the horizontal vorticity is assumed to have no vertical shear.
    !$OMP PARALLEL
    !$OMP WORKSHARE
    diag%zeta_x(:,:,1) = diag%zeta_x(:,:,2)
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! dividing by the area
    !$OMP PARALLEL
    !$OMP WORKSHARE
    diag%zeta_x = diag%zeta_x/grid%area_dual_x
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    
    ! calculating the relative vorticity in y-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jl=2,nlays
        do jk=2,ncols
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
          + grid%dz(ji,ncols,jl)*state%wind_w(ji,ncols,jl) &
          + grid%dx(ji,1,jl-1)*horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,1,jl-1)
          diag%zeta_y(ji,ncols+1,jl) = diag%zeta_y(ji,1,jl)
        else
          diag%zeta_y(ji,1,jl) = 0._wp
          diag%zeta_y(ji,ncols+1,jl) = 0._wp
        endif
        
        do jk=1,ncols+1
          ! At the surface, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
          diag%zeta_y(ji,jk,nlays+1) = grid%dx(ji,jk,nlays)*horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,jk,nlays)
        enddo
        
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! At the TOA, the horizontal vorticity is assumed to have no vertical shear.
    !$OMP PARALLEL
    !$OMP WORKSHARE
    diag%zeta_y(:,:,1) = diag%zeta_y(:,:,2)
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
    ! dividing by the area
    !$OMP PARALLEL
    !$OMP WORKSHARE
    diag%zeta_y = diag%zeta_y/grid%area_dual_y
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
      
    ! calculating the relative vorticity in z-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins+1
      do jk=1,ncols+1
        do jl=1,nlays
          diag%zeta_z(ji,jk,jl) = rel_vort_z_local(state,grid,ji,jk,jl)/grid%area_dual_z(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
      
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
      !$OMP PARALLEL
      !$OMP WORKSHARE
      diag%zeta_x = 0._wp
      diag%zeta_y = 0._wp
      diag%zeta_z = 0._wp
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
    
    ! adding the Coriolis vector to the relative vorticity to obtain the absolute vorticity
    if (lcorio) then
      !$OMP PARALLEL
      !$OMP DO PRIVATE(jl)
      do jl=1,nlays+1
        diag%eta_x(:,:,jl) = diag%zeta_x(:,:,jl) + grid%fvec_x(:,:)
        diag%eta_y(:,:,jl) = diag%zeta_y(:,:,jl) + grid%fvec_y(:,:)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP DO PRIVATE(jl)
      do jl=1,nlays
        diag%eta_z(:,:,jl) = diag%zeta_z(:,:,jl) + grid%fvec_z(:,:)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    else
      !$OMP PARALLEL
      !$OMP DO PRIVATE(jl)
      do jl=1,nlays+1
        diag%eta_x(:,:,jl) = diag%zeta_x(:,:,jl)
        diag%eta_y(:,:,jl) = diag%zeta_y(:,:,jl)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP DO PRIVATE(jl)
      do jl=1,nlays
        diag%eta_z(:,:,jl) = diag%zeta_z(:,:,jl)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    endif
    
    ! dividing by the averaged density to obtain the "potential vorticity"
    ! horizontal vorticity in x-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do jk=1,ncols
      do jl=1,nlays+1
        do ji=2,nlins
          if (jl==1) then
            diag%eta_x(ji,jk,jl) = diag%eta_x(ji,jk,jl)/(0.5_wp*(state%rho(ji-1,jk,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji-1,jk,jl))*(state%rho(ji-1,jk,jl,no_of_condensed_constituents+1) &
            -state%rho(ji-1,jk,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji-1,jk,jl)-grid%z_scalar(ji-1,jk,jl+1)) &
            +state%rho(ji,jk,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,jk,jl))*(state%rho(ji,jk,jl,no_of_condensed_constituents+1) &
            -state%rho(ji,jk,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))))
          elseif (jl==nlays+1) then
            diag%eta_x(ji,jk,jl) = diag%eta_x(ji,jk,jl)/(0.5_wp*(state%rho(ji-1,jk,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji-1,jk,jl)-grid%z_scalar(ji-1,jk,jl-1))*(state%rho(ji-1,jk,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji-1,jk,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji-1,jk,jl-2)-grid%z_scalar(ji-1,jk,jl-1)) &
            +state%rho(ji,jk,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,jk,jl)-grid%z_scalar(ji,jk,jl-1)) &
            *(state%rho(ji,jk,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji,jk,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk,jl-2)-grid%z_scalar(ji,jk,jl-1))))
          else
            diag%eta_x(ji,jk,jl) = diag%eta_x(ji,jk,jl)/(0.25_wp*(state%rho(ji-1,jk,jl-1,no_of_condensed_constituents+1) &
            +state%rho(ji,jk,jl-1,no_of_condensed_constituents+1)+ &
            state%rho(ji-1,jk,jl,no_of_condensed_constituents+1)+state%rho(ji,jk,jl,no_of_condensed_constituents+1)))
          endif
        enddo
        
        if (lperiodic) then
          if (jl==1) then
            diag%eta_x(1,jk,jl) = diag%eta_x(1,jk,jl)/(0.5_wp*(state%rho(nlins,jk,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(nlins,jk,jl))*(state%rho(nlins,jk,jl,no_of_condensed_constituents+1) &
            -state%rho(nlins,jk,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(nlins,jk,jl)-grid%z_scalar(nlins,jk,jl+1)) &
            +state%rho(1,jk,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(1,jk,jl))*(state%rho(1,jk,jl,no_of_condensed_constituents+1) &
            -state%rho(1,jk,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(1,jk,jl)-grid%z_scalar(1,jk,jl+1))))
          elseif (jl==nlays+1) then
            diag%eta_x(1,jk,jl) = diag%eta_x(1,jk,jl)/(0.5_wp*(state%rho(nlins,jk,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(nlins,jk,jl)-grid%z_scalar(nlins,jk,jl-1))*(state%rho(nlins,jk,jl-2,no_of_condensed_constituents+1) &
            -state%rho(nlins,jk,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(nlins,jk,jl-2)-grid%z_scalar(nlins,jk,jl-1)) &
            +state%rho(1,jk,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(1,jk,jl)-grid%z_scalar(1,jk,jl-1)) &
            *(state%rho(1,jk,jl-2,no_of_condensed_constituents+1) &
            -state%rho(1,jk,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(1,jk,jl-2)-grid%z_scalar(1,jk,jl-1))))
          else
            diag%eta_x(1,jk,jl) = diag%eta_x(1,jk,jl)/(0.25_wp*(state%rho(nlins,jk,jl-1,no_of_condensed_constituents+1) &
            +state%rho(1,jk,jl-1,no_of_condensed_constituents+1)+ &
            state%rho(nlins,jk,jl,no_of_condensed_constituents+1)+state%rho(1,jk,jl,no_of_condensed_constituents+1)))
          endif
          diag%eta_x(nlins+1,jk,jl) = diag%eta_x(1,jk,jl)
        endif
        
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    ! horizontal vorticity in y-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jl=1,nlays+1
        do jk=2,ncols
          if (jl==1) then
            diag%eta_y(ji,jk,jl) = diag%eta_y(ji,jk,jl)/(0.5_wp*(state%rho(ji,jk-1,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,jk-1,jl))*(state%rho(ji,jk-1,jl,no_of_condensed_constituents+1) &
            -state%rho(ji,jk-1,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk-1,jl)-grid%z_scalar(ji,jk-1,jl+1)) &
            +state%rho(ji,jk,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,jk,jl))*(state%rho(ji,jk,jl,no_of_condensed_constituents+1) &
            -state%rho(ji,jk,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1))))
          elseif (jl==nlays+1) then
            diag%eta_y(ji,jk,jl) = diag%eta_y(ji,jk,jl)/(0.5_wp*(state%rho(ji,jk-1,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,jk-1,jl)-grid%z_scalar(ji,jk-1,jl-1))*(state%rho(ji,jk-1,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji,jk-1,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk-1,jl-2)-grid%z_scalar(ji,jk-1,jl-1)) &
            +state%rho(ji,jk,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,jk,jl)-grid%z_scalar(ji,jk,jl-1)) &
            *(state%rho(ji,jk,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji,jk,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,jk,jl-2)-grid%z_scalar(ji,jk,jl-1))))
          else
            diag%eta_y(ji,jk,jl) = diag%eta_y(ji,jk,jl)/(0.25_wp*(state%rho(ji,jk-1,jl-1,no_of_condensed_constituents+1) &
            +state%rho(ji,jk,jl-1,no_of_condensed_constituents+1)+ &
            state%rho(ji,jk-1,jl,no_of_condensed_constituents+1)+state%rho(ji,jk,jl,no_of_condensed_constituents+1)))
          endif
        enddo
        
        if (lperiodic) then
          if (jl==1) then
            diag%eta_y(ji,1,jl) = diag%eta_y(ji,1,jl)/(0.5_wp*(state%rho(ji,ncols,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,ncols,jl))*(state%rho(ji,ncols,jl,no_of_condensed_constituents+1) &
            -state%rho(ji,ncols,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,ncols,jl)-grid%z_scalar(ji,ncols,jl+1)) &
            +state%rho(ji,1,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_scalar(ji,1,jl))*(state%rho(ji,1,jl,no_of_condensed_constituents+1) &
            -state%rho(ji,1,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,1,jl)-grid%z_scalar(ji,1,jl+1))))
          elseif (jl==nlays+1) then
            diag%eta_y(ji,1,jl) = diag%eta_y(ji,1,jl)/(0.5_wp*(state%rho(ji,ncols,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,ncols,jl)-grid%z_scalar(ji,ncols,jl-1))*(state%rho(ji,ncols,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji,ncols,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,ncols,jl-2)-grid%z_scalar(ji,ncols,jl-1)) &
            +state%rho(ji,1,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_w(ji,1,jl)-grid%z_scalar(ji,1,jl-1)) &
            *(state%rho(ji,1,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji,1,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_scalar(ji,1,jl-2)-grid%z_scalar(ji,1,jl-1))))
          else
            diag%eta_y(ji,1,jl) = diag%eta_y(ji,1,jl)/(0.25_wp*(state%rho(ji,ncols,jl-1,no_of_condensed_constituents+1) &
            +state%rho(ji,1,jl-1,no_of_condensed_constituents+1)+ &
            state%rho(ji,ncols,jl,no_of_condensed_constituents+1)+state%rho(ji,1,jl,no_of_condensed_constituents+1)))
          endif
          diag%eta_y(ji,ncols+1,jl) = diag%eta_y(ji,1,jl)
        endif
        
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! vertical vorticity
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=2,nlins
      do jk=2,ncols
        diag%eta_z(ji,jk,:) = diag%eta_z(ji,jk,:)/(0.25_wp*(state%rho(ji,jk-1,:,no_of_condensed_constituents+1) &
        +state%rho(ji,jk,:,no_of_condensed_constituents+1) &
        +state%rho(ji-1,jk,:,no_of_condensed_constituents+1)+state%rho(ji-1,jk-1,:,no_of_condensed_constituents+1)))
      enddo
      
      ! periodic boundary conditions
      if (lperiodic) then
        diag%eta_z(ji,1,:) = diag%eta_z(ji,1,:)/(0.25_wp*(state%rho(ji,ncols,:,no_of_condensed_constituents+1) &
        +state%rho(ji,1,:,no_of_condensed_constituents+1) &
        +state%rho(ji-1,1,:,no_of_condensed_constituents+1)+state%rho(ji-1,ncols,:,no_of_condensed_constituents+1)))
        diag%eta_z(ji,ncols+1,:) = diag%eta_z(ji,1,:)
      endif
      
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! periodic boundary conditions
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP DO PRIVATE(jk)
      do jk=2,ncols
        diag%eta_z(1,jk,:) = diag%eta_z(1,jk,:)/(0.25_wp*(state%rho(1,jk-1,:,no_of_condensed_constituents+1) &
        +state%rho(1,jk,:,no_of_condensed_constituents+1) &
        +state%rho(nlins,jk,:,no_of_condensed_constituents+1)+state%rho(nlins,jk-1,:,no_of_condensed_constituents+1)))
        diag%eta_z(nlins+1,jk,:) = diag%eta_z(1,jk,:)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    endif
    
    ! corners under periodic boundary conditions
    if (lperiodic) then
      diag%eta_z(1,1,:) = diag%eta_z(1,1,:)/(0.25_wp*(state%rho(1,ncols,:,no_of_condensed_constituents+1) &
      +state%rho(1,1,:,no_of_condensed_constituents+1) &
      +state%rho(nlins,1,:,no_of_condensed_constituents+1)+state%rho(nlins,ncols,:,no_of_condensed_constituents+1)))
      diag%eta_z(1,ncols+1,:) = diag%eta_z(1,1,:)
      diag%eta_z(nlins+1,1,:) = diag%eta_z(1,1,:)
      diag%eta_z(nlins+1,ncols+1,:) = diag%eta_z(1,1,:)
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
    if (ji==1 .or. jk==1 .or. ji==nlins+1 .or. jk==ncols+1) then
      if (lperiodic) then
        j_i(1) = ji
        if (jk==ncols+1) then
          j_k(1) = 1
        else
          j_k(1) = jk
        endif
        
        if (ji==1) then
          j_i(2) = nlins
        else
          j_i(2) = ji-1
        endif
        j_k(2) = jk
        
        j_i(3) = ji
        if (jk==1) then
          j_k(3) = ncols
        else
          j_k(3) = jk-1
        endif
        
        if (ji==nlins+1) then
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
    if (jl<=nlays-nlays_oro) then
      rel_vort_z_local = rel_vort_z_local &
      + sign_vector(1)*grid%dy(j_i(1),j_k(1),jl)*state%wind_v(j_i(1),j_k(1),jl) &
      + sign_vector(2)*grid%dx(j_i(2),j_k(2),jl)*state%wind_u(j_i(2),j_k(2),jl) &
      + sign_vector(3)*grid%dy(j_i(3),j_k(3),jl)*state%wind_v(j_i(3),j_k(3),jl) &
      + sign_vector(4)*grid%dx(j_i(4),j_k(4),jl)*state%wind_u(j_i(4),j_k(4),jl)
    ! layers which follow the orography
    else
      do jm=1,4
        if (jm==1 .or. jm==3) then
          l_rescale = (re + grid%z_area_dual_z(ji,jk,jl))/(re + grid%z_v(j_i(jm),j_k(jm),jl))
          delta_z = grid%z_area_dual_z(ji,jk,jl) - grid%z_v(j_i(jm),j_k(jm),jl)
          ind_shift = 1
          if (delta_z>0._wp .or. jl==nlays) then
            ind_shift = -1
          endif
          vertical_gradient = (state%wind_v(j_i(jm),j_k(jm),jl) - state%wind_v(j_i(jm),j_k(jm),jl+ind_shift))/ &
          (grid%z_v(j_i(jm),j_k(jm),jl) - grid%z_v(j_i(jm),j_k(jm),jl+ind_shift))
          rel_vort_z_local = rel_vort_z_local + sign_vector(jm)*l_rescale*grid%dy(j_i(jm),j_k(jm),jl)* &
          (state%wind_v(j_i(jm),j_k(jm),jl) + delta_z*vertical_gradient)
        else
          l_rescale = (re + grid%z_area_dual_z(ji,jk,jl))/(re + grid%z_u(j_i(jm),j_k(jm),jl))
          delta_z = grid%z_area_dual_z(ji,jk,jl) - grid%z_u(j_i(jm),j_k(jm),jl)
          ind_shift = 1
          if (delta_z>0._wp .or. jl==nlays) then
            ind_shift = -1
          endif
          vertical_gradient = (state%wind_u(j_i(jm),j_k(jm),jl) - state%wind_u(j_i(jm),j_k(jm),jl+ind_shift))/ &
          (grid%z_u(j_i(jm),j_k(jm),jl) - grid%z_u(j_i(jm),j_k(jm),jl+ind_shift))
          rel_vort_z_local = rel_vort_z_local + sign_vector(jm)*l_rescale*grid%dx(j_i(jm),j_k(jm),jl)* &
          (state%wind_u(j_i(jm),j_k(jm),jl) + delta_z*vertical_gradient)
        endif
      enddo
    endif
  
  end function rel_vort_z_local

end module vorticities








