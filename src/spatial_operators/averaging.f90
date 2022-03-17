! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module averaging

  ! This module contains averaging operators.

  use definitions, only: t_grid,wp
  use run_nml,     only: nlays,nlays_oro,nlins,ncols
  use bc_nml,      only: lperiodic
  
  implicit none
  
  private
  
  public :: vertical_contravariant_corr
  public :: hor_cov_to_con
  public :: horizontal_covariant_x
  public :: horizontal_covariant_y
  
  contains
  
  function vertical_contravariant_corr(vector_field_x,vector_field_y,ji,jk,jl,grid)
  
    ! calculates (the vertical contravariant component - the vertical covariant component)
    ! of a vector field out of the horizontal contravariant components
    
    real(wp),     intent(in) :: vector_field_x(:,:,:) ! x-component of vector field to work with
    real(wp),     intent(in) :: vector_field_y(:,:,:) ! y-component of vector field to work with
    type(t_grid), intent(in) :: grid                  ! model grid
    integer,      intent(in) :: ji,jk,jl              ! spatial indices
    
    real(wp) :: vertical_contravariant_corr
    vertical_contravariant_corr = 0._wp
    
    if (jl>=nlays-nlays_oro+1) then
      ! highest level following orography
      if (jl==nlays-nlays_oro+1) then
        vertical_contravariant_corr = vertical_contravariant_corr &
        - 0.5_wp*vector_field_x(ji,jk+1,jl)*grid%slope_x(ji,jk+1,jl)*grid%inner_product_weights(ji,jk,jl,1) &
        - 0.5_wp*vector_field_y(ji,jk,jl)*grid%slope_y(ji,jk,jl)*grid%inner_product_weights(ji,jk,jl,2) &
        - 0.5_wp*vector_field_x(ji,jk,jl)*grid%slope_x(ji,jk,jl)*grid%inner_product_weights(ji,jk,jl,3) &
        - 0.5_wp*vector_field_y(ji+1,jk,jl)*grid%slope_y(ji+1,jk,jl)*grid%inner_product_weights(ji,jk,jl,4)
      ! levels below
      else
        vertical_contravariant_corr = vertical_contravariant_corr &
        - 0.5_wp*vector_field_x(ji,jk+1,jl-1)*grid%slope_x(ji,jk+1,jl-1)*grid%inner_product_weights(ji,jk,jl-1,1) &
        - 0.5_wp*vector_field_y(ji,jk,jl-1)*grid%slope_y(ji,jk,jl-1)*grid%inner_product_weights(ji,jk,jl-1,2) &
        - 0.5_wp*vector_field_x(ji,jk,jl-1)*grid%slope_x(ji,jk,jl-1)*grid%inner_product_weights(ji,jk,jl-1,3) &
        - 0.5_wp*vector_field_y(ji+1,jk,jl-1)*grid%slope_y(ji+1,jk,jl-1)*grid%inner_product_weights(ji,jk,jl-1,4)
        vertical_contravariant_corr = vertical_contravariant_corr &
        - 0.5_wp*vector_field_x(ji,jk+1,jl)*grid%slope_x(ji,jk+1,jl)*grid%inner_product_weights(ji,jk,jl,1) &
        - 0.5_wp*vector_field_y(ji,jk,jl)*grid%slope_y(ji,jk,jl)*grid%inner_product_weights(ji,jk,jl,2) &
        - 0.5_wp*vector_field_x(ji,jk,jl)*grid%slope_x(ji,jk,jl)*grid%inner_product_weights(ji,jk,jl,3) &
        - 0.5_wp*vector_field_y(ji+1,jk,jl)*grid%slope_y(ji+1,jk,jl)*grid%inner_product_weights(ji,jk,jl,4)
      endif
    endif
  
  end function vertical_contravariant_corr
  
  subroutine hor_cov_to_con(result_field_x,result_field_y,result_field_z,grid)
  
    ! This subroutine computes the terrain correction of the gradient.
    real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
    real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
    real(wp),     intent(in)    :: result_field_z(:,:,:) ! z-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
  
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! correction to the x-component
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols+1
        do jl=1,nlays
          result_field_x(ji,jk,jl) = result_field_x(ji,jk,jl) &
          - grid%slope_x(ji,jk,jl)*remap_ver2hor_x(result_field_z,grid,ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! correction to the y-component
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins+1
      do jk=1,ncols
        do jl=1,nlays
          result_field_y(ji,jk,jl) = result_field_y(ji,jk,jl) &
          - grid%slope_y(ji,jk,jl)*remap_ver2hor_y(result_field_z,grid,ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine hor_cov_to_con
  
  function remap_ver2hor_x(vertical_cov,grid,ji,jk,jl)
  
    ! This function remaps a vertical covariant component of
    ! a vector field to a position of a vector component in x-direction.
  
    real(wp),     intent(in) :: vertical_cov(:,:,:) ! z-component of vector field to work with
    type(t_grid), intent(in) :: grid                ! the grid properties
    integer,      intent(in) :: ji,jk,jl            ! positional indices
  
    real(wp) :: remap_ver2hor_x ! the result
    
    ! initialization with zero
    remap_ver2hor_x = 0._wp
    
    if (jk==1 .or. jk==ncols+1) then
      if (lperiodic) then
        remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,ncols,jl,5)*vertical_cov(ji,ncols,jl)
        remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,1,jl,5)*vertical_cov(ji,1,jl)
        ! layer below
        if (jl<nlays) then
          remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,ncols,jl,6)*vertical_cov(ji,ncols,jl+1)
          remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,1,jl,6)*vertical_cov(ji,1,jl+1)
        endif
      else
        return
      endif
    else
      remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,jk-1,jl,5)*vertical_cov(ji,jk-1,jl)
      remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,jk,jl,5)*vertical_cov(ji,jk,jl)
      ! layer below
      if (jl<nlays) then
        remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,jk-1,jl,6)*vertical_cov(ji,jk-1,jl+1)
        remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,jk,jl,6)*vertical_cov(ji,jk,jl+1)
      endif
    endif
    
    ! horizontal average
    remap_ver2hor_x = 0.5_wp*remap_ver2hor_x
  
  end function remap_ver2hor_x
  
  function remap_ver2hor_y(vertical_cov,grid,ji,jk,jl)
  
    ! This function remaps a vertical covariant component of
    ! a vector field to a position of a vector component in y-direction.
  
    real(wp),     intent(in) :: vertical_cov(:,:,:) ! z-component of vector field to work with
    type(t_grid), intent(in) :: grid                ! the grid properties
    integer,      intent(in) :: ji,jk,jl            ! positional indices
  
    real(wp) :: remap_ver2hor_y ! the result
    
    ! initialization with zero
    remap_ver2hor_y = 0._wp
    
    if (ji==1 .or. ji==nlins+1) then
      if (lperiodic) then
        remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(nlins,jk,jl,5)*vertical_cov(nlins,jk,jl)
        remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(1,jk,jl,5)*vertical_cov(1,jk,jl)
        ! layer below
        if (jl<nlays) then
          remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(nlins,jk,jl,6)*vertical_cov(nlins,jk,jl+1)
          remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(1,jk,jl,6)*vertical_cov(1,jk,jl+1)
        endif
      else
        return
      endif
    else
      remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ji-1,jk,jl,5)*vertical_cov(ji-1,jk,jl)
      remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ji,jk,jl,5)*vertical_cov(ji,jk,jl)
      ! layer below
      if (jl<nlays) then
        remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ji-1,jk,jl,6)*vertical_cov(ji-1,jk,jl+1)
        remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ji,jk,jl,6)*vertical_cov(ji,jk,jl+1)
      endif
    endif
    
    ! horizontal average
    remap_ver2hor_y = 0.5_wp*remap_ver2hor_y
  
  end function remap_ver2hor_y
  
  function horizontal_covariant_x(hor_comp_x,vert_comp,grid,ji,jk,jl)
  
    ! This function calculates the horizontal covariant component of a vector field in x-direction.
  
    real(wp),     intent(in) :: hor_comp_x(:,:,:) ! horizontal component in x-direction of vector field to work with
    real(wp),     intent(in) :: vert_comp(:,:,:)  ! vertical component of vector field to work with
    type(t_grid), intent(in) :: grid              ! model grid
    integer,      intent(in) :: ji,jk,jl          ! positional indices
    
    ! output
    real(wp) :: horizontal_covariant_x
    
    horizontal_covariant_x = hor_comp_x(ji,jk,jl)
    
    if (jl>nlays-nlays_oro) then
      horizontal_covariant_x = horizontal_covariant_x + grid%slope_x(ji,jk,jl)*remap_ver2hor_x(vert_comp,grid,ji,jk,jl)
    endif
    
  end function horizontal_covariant_x
  
  function horizontal_covariant_y(hor_comp_y,vert_comp,grid,ji,jk,jl)
  
    ! This function calculates the horizontal covariant component of a vector field in y-direction.
  
    real(wp),     intent(in) :: hor_comp_y(:,:,:) ! horizontal component in x-direction of vector field to work with
    real(wp),     intent(in) :: vert_comp(:,:,:)  ! vertical component of vector field to work with
    type(t_grid), intent(in) :: grid              ! model grid
    integer,      intent(in) :: ji,jk,jl          ! positional indices
    
    ! output
    real(wp) :: horizontal_covariant_y
  
    horizontal_covariant_y = hor_comp_y(ji,jk,jl)
    
    if (jl>nlays-nlays_oro) then
      horizontal_covariant_y = horizontal_covariant_y + grid%slope_y(ji,jk,jl)*remap_ver2hor_y(vert_comp,grid,ji,jk,jl)
    endif
  
  end function horizontal_covariant_y

end module averaging







