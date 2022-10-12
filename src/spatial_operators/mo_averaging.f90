! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_averaging

  ! This module contains averaging operators.

  use mo_definitions, only: t_grid,wp
  use mo_run_nml,     only: n_layers,n_oro_layers,ny,nx,n_flat_layers
  use mo_bc_nml,      only: lperiodic
  
  implicit none
  
  contains
  
  function vertical_contravariant_corr(vector_field_x,vector_field_y,ji,jk,jl,grid)
  
    ! This function calculates (the vertical contravariant component - the vertical covariant component)
    ! of a vector field out of the horizontal contravariant components.
    
    real(wp),     intent(in) :: vector_field_x(:,:,:) ! x-component of vector field to work with
    real(wp),     intent(in) :: vector_field_y(:,:,:) ! y-component of vector field to work with
    type(t_grid), intent(in) :: grid                  ! model grid
    integer,      intent(in) :: ji,jk,jl              ! spatial indices
    real(wp) :: vertical_contravariant_corr           ! result
        
    vertical_contravariant_corr = 0._wp
    
    if (jl>=n_flat_layers+1) then
      ! highest level following orography
      if (jl==n_flat_layers+1) then
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
  
  function remap_ver2hor_x(vertical_cov,grid,ji,jk,jl)
  
    ! This function remaps a vertical covariant component of a vector field to a position of a vector component in x-direction.
  
    real(wp),     intent(in) :: vertical_cov(:,:,:) ! z-component of vector field to work with
    type(t_grid), intent(in) :: grid                ! the grid properties
    integer,      intent(in) :: ji,jk,jl            ! positional indices
    real(wp)                 :: remap_ver2hor_x     ! result
    
    ! initialization with zero
    remap_ver2hor_x = 0._wp
    
    if (jk==1 .or. jk==nx+1) then
      if (lperiodic) then
        remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,nx,jl,5)*vertical_cov(ji,nx,jl)
        remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,1,jl,5)*vertical_cov(ji,1,jl)
        ! layer below
        if (jl<n_layers) then
          remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,nx,jl,6)*vertical_cov(ji,nx,jl+1)
          remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,1,jl,6)*vertical_cov(ji,1,jl+1)
        endif
      else
        return
      endif
    else
      remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,jk-1,jl,5)*vertical_cov(ji,jk-1,jl)
      remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,jk,jl,5)*vertical_cov(ji,jk,jl)
      ! layer below
      if (jl<n_layers) then
        remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,jk-1,jl,6)*vertical_cov(ji,jk-1,jl+1)
        remap_ver2hor_x = remap_ver2hor_x + grid%inner_product_weights(ji,jk,jl,6)*vertical_cov(ji,jk,jl+1)
      endif
    endif
    
    ! horizontal average
    remap_ver2hor_x = 0.5_wp*remap_ver2hor_x
  
  end function remap_ver2hor_x
  
  function remap_ver2hor_y(vertical_cov,grid,ji,jk,jl)
  
    ! This function remaps a vertical covariant component of a vector field to a position of a vector component in y-direction.
  
    real(wp),     intent(in) :: vertical_cov(:,:,:) ! z-component of vector field to work with
    type(t_grid), intent(in) :: grid                ! the grid properties
    integer,      intent(in) :: ji,jk,jl            ! positional indices
    real(wp)                 :: remap_ver2hor_y     ! result
    
    ! initialization with zero
    remap_ver2hor_y = 0._wp
    
    if (ji==1 .or. ji==ny+1) then
      if (lperiodic) then
        remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ny,jk,jl,5)*vertical_cov(ny,jk,jl)
        remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(1,jk,jl,5)*vertical_cov(1,jk,jl)
        ! layer below
        if (jl<n_layers) then
          remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ny,jk,jl,6)*vertical_cov(ny,jk,jl+1)
          remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(1,jk,jl,6)*vertical_cov(1,jk,jl+1)
        endif
      else
        return
      endif
    else
      remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ji-1,jk,jl,5)*vertical_cov(ji-1,jk,jl)
      remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ji,jk,jl,5)*vertical_cov(ji,jk,jl)
      ! layer below
      if (jl<n_layers) then
        remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ji-1,jk,jl,6)*vertical_cov(ji-1,jk,jl+1)
        remap_ver2hor_y = remap_ver2hor_y + grid%inner_product_weights(ji,jk,jl,6)*vertical_cov(ji,jk,jl+1)
      endif
    endif
    
    ! horizontal average
    remap_ver2hor_y = 0.5_wp*remap_ver2hor_y
  
  end function remap_ver2hor_y
  
  function horizontal_covariant_x(hor_comp_x,vert_comp,grid,ji,jk,jl)
  
    ! This function calculates the horizontal covariant component of a vector field in x-direction.
  
    real(wp),     intent(in) :: hor_comp_x(:,:,:)      ! horizontal component in x-direction of vector field to work with
    real(wp),     intent(in) :: vert_comp(:,:,:)       ! vertical component of vector field to work with
    type(t_grid), intent(in) :: grid                   ! model grid
    integer,      intent(in) :: ji,jk,jl               ! positional indices
    real(wp)                 :: horizontal_covariant_x ! result
    
    horizontal_covariant_x = hor_comp_x(ji,jk,jl)
    
    if (jl>n_flat_layers) then
      horizontal_covariant_x = horizontal_covariant_x + grid%slope_x(ji,jk,jl)*remap_ver2hor_x(vert_comp,grid,ji,jk,jl)
    endif
    
  end function horizontal_covariant_x
  
  function horizontal_covariant_y(hor_comp_y,vert_comp,grid,ji,jk,jl)
  
    ! This function calculates the horizontal covariant component of a vector field in y-direction.
  
    real(wp),     intent(in) :: hor_comp_y(:,:,:)      ! horizontal component in x-direction of vector field to work with
    real(wp),     intent(in) :: vert_comp(:,:,:)       ! vertical component of vector field to work with
    type(t_grid), intent(in) :: grid                   ! model grid
    integer,      intent(in) :: ji,jk,jl               ! positional indices
    real(wp)                 :: horizontal_covariant_y ! result
  
    horizontal_covariant_y = hor_comp_y(ji,jk,jl)
    
    if (jl>n_flat_layers) then
      horizontal_covariant_y = horizontal_covariant_y + grid%slope_y(ji,jk,jl)*remap_ver2hor_y(vert_comp,grid,ji,jk,jl)
    endif
  
  end function horizontal_covariant_y

end module mo_averaging







