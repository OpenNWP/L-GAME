! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module vorticities

  ! This module contians the calculation of the vorticities.

  use definitions,      only: t_state,t_diag,t_grid,wp
  use run_nml,          only: nlins,ncols,nlays,nlays_oro,toa,lcorio,llinear
  use constants,        only: re
  use constituents_nml, only: no_of_condensed_constituents
  use averaging,        only: horizontal_covariant_x,horizontal_covariant_y
  
  implicit none
  
  private
  
  public :: calc_pot_vort
  
  contains

  subroutine rel_vort(state,diag,grid)
  
    ! This subroutine calculates the relative vorticity.
    
    type(t_state), intent(in)    :: state ! state to work with
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! model grid
    
    ! local variables
    real(wp) :: l_rescale         ! length rescale factor in orography
    real(wp) :: delta_z           ! needed for handling terrain
    real(wp) :: vertical_gradient ! needed for handling terrain
    integer  :: ind_shift         ! needed for handling terrain
    integer  :: ji,jk,jl          ! loop indices
    
    ! calculating the relative vorticity in x-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=2,nlins
      do jk=1,ncols
        do jl=2,nlays
          diag%z_eta_x(ji,jk,jl) = &
          + grid%dz(ji,jk,jl)*state%wind_w(ji,jk,jl) &
          - grid%dy(ji,jk,jl-1)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,ji,jk,jl-1) &
          - grid%dz(ji-1,jk,jl)*state%wind_w(ji-1,jk,jl) &
          + grid%dy(ji,jk,jl)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,ji,jk,jl)
        enddo
        ! At the surface, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
        diag%z_eta_x(ji,jk,nlays+1) = -grid%dy(ji,jk,nlays)*horizontal_covariant_y(state%wind_v,state%wind_w,grid,ji,jk,nlays)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! At the TOA, the horizontal vorticity is assumed to have no vertical shear.
    diag%z_eta_x(:,:,1) = diag%z_eta_x(:,:,2)
    ! dividing by the area
    diag%z_eta_x = diag%z_eta_x/grid%area_dual_x
    
    ! calculating the relative vorticity in y-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=2,ncols
        do jl=2,nlays
          diag%z_eta_y(ji,jk,jl) = &
          - grid%dz(ji,jk,jl)*state%wind_w(ji,jk,jl) &
          - grid%dx(ji,jk,jl)*horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,jk,jl) &
          + grid%dz(ji,jk-1,jl)*state%wind_w(ji,jk-1,jl) &
          + grid%dx(ji,jk,jl-1)*horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,jk,jl-1)
        enddo
        ! At the surface, w vanishes. Furthermore, the covariant velocity below the surface is also zero.
        diag%z_eta_y(ji,jk,nlays+1) = grid%dx(ji,jk,nlays)*horizontal_covariant_x(state%wind_u,state%wind_w,grid,ji,jk,nlays)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! At the TOA, the horizontal vorticity is assumed to have no vertical shear.
    diag%z_eta_y(:,:,1) = diag%z_eta_y(:,:,2)
    ! dividing by the area
    diag%z_eta_y = diag%z_eta_y/grid%area_dual_y
      
    ! calculating the relative vorticity in z-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,l_rescale,delta_z,ind_shift,vertical_gradient)
    do ji=2,nlins
      do jk=2,ncols
        ! layers which do not follow the orography
        do jl=1,nlays-nlays_oro
          diag%z_eta_z(ji,jk,jl) = &
          + grid%dy(ji,jk,jl)*state%wind_v(ji,jk,jl) &
          - grid%dx(ji,jk,jl)*state%wind_u(ji,jk,jl) &
          - grid%dy(ji,jk-1,jl)*state%wind_v(ji,jk-1,jl) &
          + grid%dx(ji-1,jk,jl)*state%wind_u(ji-1,jk,jl)
        enddo
        ! layers which follow the orography
        do jl=nlays-nlays_oro+1,nlays
          ! first
          l_rescale = (re + grid%z_geo_area_dual_z(ji,jk,jl))/(re + grid%z_geo_v(ji,jk,jl))
          delta_z = grid%z_geo_area_dual_z(ji,jk,jl) - grid%z_geo_v(ji,jk,jl)
          ind_shift = 1
          if (delta_z>0._wp .or. jl==nlays) then
            ind_shift = -1
          endif
          vertical_gradient = (state%wind_v(ji,jk,jl) - state%wind_v(ji,jk,jl+ind_shift))/ &
          (grid%z_geo_v(ji,jk,jl) - grid%z_geo_v(ji,jk,jl+ind_shift))
          diag%z_eta_z(ji,jk,jl) = l_rescale*grid%dy(ji,jk,jl)* &
          (state%wind_v(ji,jk,jl) + delta_z*vertical_gradient)
          ! second
          l_rescale = (re + grid%z_geo_area_dual_z(ji,jk,jl))/(re + grid%z_geo_u(ji,jk,jl))
          delta_z = grid%z_geo_area_dual_z(ji,jk,jl) - grid%z_geo_u(ji,jk,jl)
          ind_shift = 1
          if (delta_z>0._wp .or. jl==nlays) then
            ind_shift = -1
          endif
          vertical_gradient = (state%wind_u(ji,jk,jl) - state%wind_u(ji,jk,jl+ind_shift))/ &
          (grid%z_geo_u(ji,jk,jl) - grid%z_geo_u(ji,jk,jl+ind_shift))
          diag%z_eta_z(ji,jk,jl) = diag%z_eta_z(ji,jk,jl) - l_rescale*grid%dx(ji,jk,jl)* &
          (state%wind_u(ji,jk,jl) + delta_z*vertical_gradient)
          ! third
          l_rescale = (re + grid%z_geo_area_dual_z(ji,jk-1,jl))/(re + grid%z_geo_v(ji,jk-1,jl))
          delta_z = grid%z_geo_area_dual_z(ji,jk-1,jl) - grid%z_geo_v(ji,jk-1,jl)
          ind_shift = 1
          if (delta_z>0._wp .or. jl==nlays) then
            ind_shift = -1
          endif
          vertical_gradient = (state%wind_v(ji,jk-1,jl) - state%wind_v(ji,jk-1,jl+ind_shift))/ &
          (grid%z_geo_v(ji,jk-1,jl) - grid%z_geo_v(ji,jk-1,jl+ind_shift))
          diag%z_eta_z(ji,jk-1,jl) = diag%z_eta_z(ji,jk-1,jl) - l_rescale*grid%dy(ji,jk-1,jl)* &
          (state%wind_v(ji,jk-1,jl) + delta_z*vertical_gradient)
          ! fourth
          l_rescale = (re + grid%z_geo_area_dual_z(ji-1,jk,jl))/(re + grid%z_geo_u(ji-1,jk,jl))
          delta_z = grid%z_geo_area_dual_z(ji-1,jk,jl) - grid%z_geo_u(ji-1,jk,jl)
          ind_shift = 1
          if (delta_z>0._wp .or. jl==nlays) then
            ind_shift = -1
          endif
          vertical_gradient = (state%wind_u(ji-1,jk,jl) - state%wind_u(ji-1,jk,jl+ind_shift))/ &
          (grid%z_geo_u(ji-1,jk,jl) - grid%z_geo_u(ji-1,jk,jl+ind_shift))
          diag%z_eta_z(ji-1,jk,jl) = diag%z_eta_z(ji-1,jk,jl) + l_rescale*grid%dx(ji-1,jk,jl)* &
          (state%wind_u(ji-1,jk,jl) + delta_z*vertical_gradient)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! dividing by the area
    diag%z_eta_z = diag%z_eta_z/grid%area_dual_z
      
  end subroutine rel_vort
  
  subroutine calc_pot_vort(state,diag,grid)
  
    ! This subroutine calculates the potential vorticity.
    
    type(t_state), intent(in)    :: state ! state to work with
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! model grid
  
    ! local variables
    integer :: ji,jk,jl ! loop index
  
    ! calculating the relative vorticity
    if (.not. llinear) then
      call rel_vort(state,diag,grid)
    else
      diag%z_eta_x = 0._wp
      diag%z_eta_y = 0._wp
      diag%z_eta_z = 0._wp
    endif
    
    ! adding the Coriolis vector to the relative vorticity to obtain the absolute vorticity
    if (lcorio) then
      !$OMP PARALLEL
      !$OMP DO PRIVATE(jl)
      do jl=1,nlays+1
        diag%z_eta_x(:,:,jl) = diag%z_eta_x(:,:,jl) + grid%fvec_x(:,:)
        diag%z_eta_y(:,:,jl) = diag%z_eta_y(:,:,jl) + grid%fvec_y(:,:)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP DO PRIVATE(jl)
      do jl=1,nlays
        diag%z_eta_z(:,:,jl) = diag%z_eta_z(:,:,jl) + grid%fvec_z(:,:)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
    endif
    
    ! dividing by the averaged density to obtain the "potential vorticity"
    ! horizontal vorticity in x-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins+1
      do jk=1,ncols
        do jl=1,nlays+1
          if (jl == 1) then
            diag%z_eta_x(ji,jk,jl) = diag%z_eta_x(ji,jk,jl)/(0.5_wp*(state%rho(ji,jk+1,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_geo_scal(ji,jk+1,jl))*(state%rho(ji,jk+1,jl,no_of_condensed_constituents+1) &
            -state%rho(ji,jk+1,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_geo_scal(ji,jk+1,jl)-grid%z_geo_scal(ji,jk+1,jl+1)) &
            +state%rho(ji+1,jk+1,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_geo_scal(ji+1,jk+1,jl))*(state%rho(ji+1,jk+1,jl,no_of_condensed_constituents+1) &
            -state%rho(ji+1,jk+1,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_geo_scal(ji+1,jk+1,jl)-grid%z_geo_scal(ji+1,jk+1,jl+1))))
          elseif (jl==nlays+1) then
            diag%z_eta_x(ji,jk,jl) = diag%z_eta_x(ji,jk,jl)/(0.5_wp*(state%rho(ji,jk+1,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_geo_w(ji,jk+1,jl)-grid%z_geo_scal(ji,jk+1,jl-1))*(state%rho(ji,jk+1,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji,jk+1,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_geo_scal(ji,jk+1,jl-2)-grid%z_geo_scal(ji,jk+1,jl-1)) &
            +state%rho(ji+1,jk+1,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_geo_w(ji+1,jk+1,jl)-grid%z_geo_scal(ji+1,jk+1,jl-1)) &
            *(state%rho(ji+1,jk+1,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji+1,jk+1,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_geo_scal(ji+1,jk+1,jl-2)-grid%z_geo_scal(ji+1,jk+1,jl-1))))
          else
            diag%z_eta_x(ji,jk,jl) = diag%z_eta_x(ji,jk,jl)/(0.25_wp*(state%rho(ji,jk+1,jl-1,no_of_condensed_constituents+1) &
            +state%rho(ji+1,jk+1,jl-1,no_of_condensed_constituents+1)+ &
            state%rho(ji,jk+1,jl,no_of_condensed_constituents+1)+state%rho(ji+1,jk+1,jl,no_of_condensed_constituents+1)))
          endif
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    ! horizontal vorticity in y-direction
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols+1
        do jl=1,nlays+1
          if (jl == 1) then
            diag%z_eta_y(ji,jk,jl) = diag%z_eta_y(ji,jk,jl)/(0.5_wp*(state%rho(ji+1,jk,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_geo_scal(ji+1,jk,jl))*(state%rho(ji+1,jk,jl,no_of_condensed_constituents+1) &
            -state%rho(ji+1,jk,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_geo_scal(ji+1,jk,jl)-grid%z_geo_scal(ji+1,jk,jl+1)) &
            +state%rho(ji+1,jk+1,jl,no_of_condensed_constituents+1) &
            ! linear extrapolation to the TOA
            +(toa-grid%z_geo_scal(ji+1,jk+1,jl))*(state%rho(ji+1,jk+1,jl,no_of_condensed_constituents+1) &
            -state%rho(ji+1,jk+1,jl+1,no_of_condensed_constituents+1))/ &
            (grid%z_geo_scal(ji+1,jk+1,jl)-grid%z_geo_scal(ji+1,jk+1,jl+1))))
          elseif (jl==nlays+1) then
            diag%z_eta_y(ji,jk,jl) = diag%z_eta_y(ji,jk,jl)/(0.5_wp*(state%rho(ji+1,jk,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_geo_w(ji+1,jk,jl)-grid%z_geo_scal(ji+1,jk,jl-1))*(state%rho(ji+1,jk,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji+1,jk,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_geo_scal(ji+1,jk,jl-2)-grid%z_geo_scal(ji+1,jk,jl-1)) &
            +state%rho(ji+1,jk+1,jl-1,no_of_condensed_constituents+1) &
            ! linear extrapolation to the surface
            +(grid%z_geo_w(ji+1,jk+1,jl)-grid%z_geo_scal(ji+1,jk+1,jl-1)) &
            *(state%rho(ji+1,jk+1,jl-2,no_of_condensed_constituents+1) &
            -state%rho(ji+1,jk+1,jl-1,no_of_condensed_constituents+1))/ &
            (grid%z_geo_scal(ji+1,jk+1,jl-2)-grid%z_geo_scal(ji+1,jk+1,jl-1))))
          else
            diag%z_eta_y(ji,jk,jl) = diag%z_eta_y(ji,jk,jl)/(0.25_wp*(state%rho(ji+1,jk,jl-1,no_of_condensed_constituents+1) &
            +state%rho(ji+1,jk+1,jl-1,no_of_condensed_constituents+1)+ &
            state%rho(ji+1,jk,jl,no_of_condensed_constituents+1)+state%rho(ji+1,jk+1,jl,no_of_condensed_constituents+1)))
          endif
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! vertical vorticity
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins+1
      do jk=1,ncols+1
        diag%z_eta_z(ji,jk,:) = diag%z_eta_z(ji,jk,:)/(0.25_wp*(state%rho(ji,jk,:,no_of_condensed_constituents+1) &
        +state%rho(ji+1,jk,:,no_of_condensed_constituents+1) &
        +state%rho(ji+1,jk+1,:,no_of_condensed_constituents+1)+state%rho(ji,jk+1,:,no_of_condensed_constituents+1)))
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
  end subroutine calc_pot_vort

end module vorticities








