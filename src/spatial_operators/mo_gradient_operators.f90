! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_gradient_operators

  ! This module is a collection of gradient operators.

  use mo_definitions, only: t_grid,wp
  use mo_run_nml,     only: ny,nx,n_layers,toa,n_flat_layers
  use mo_averaging,   only: remap_ver2hor_x,remap_ver2hor_y
  use mo_bc_nml,      only: lperiodic
    
  implicit none
  
  contains
  
  subroutine grad_hor_cov(scalar_field,result_field_x,result_field_y,grid)

    ! This subroutine computes the horizontal covariant gradient of a scalar field.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
    real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    
    ! local variables
    integer :: ji,jk ! spatial indices

    ! inner domain
    ! calculating the x-component of the gradient
    !$omp parallel do private(jk)
    do jk=2,nx
      result_field_x(:,jk,:) = (scalar_field(:,jk,:) - scalar_field(:,jk-1,:))/grid%dx(:,jk,:)
    enddo
    !$omp end parallel do

    ! calculating the y-component of the gradient
    !$omp parallel do private(ji)
    do ji=2,ny
      result_field_y(ji,:,:) = (scalar_field(ji-1,:,:) - scalar_field(ji,:,:))/grid%dy(ji,:,:)
    enddo
    !$omp end parallel do

    ! periodic boundary conditions
    if (lperiodic) then
      !$omp parallel workshare
      result_field_x(:,1,:) = (scalar_field(:,1,:) - scalar_field(:,nx,:))/grid%dx(:,1,:)
      !$omp end parallel workshare
      
      !$omp parallel workshare
      result_field_x(:,nx+1,:) = result_field_x(:,1,:)
      !$omp end parallel workshare
      
      !$omp parallel workshare
      result_field_y(1,:,:) = (scalar_field(ny,:,:) - scalar_field(1,:,:))/grid%dy(1,:,:)
      !$omp end parallel workshare
      !$omp parallel workshare
      result_field_y(ny+1,:,:) = result_field_y(1,:,:)
      !$omp end parallel workshare
    endif

  end subroutine grad_hor_cov
  
  subroutine grad_vert(scalar_field,result_field,grid)

    ! This subroutine computes the vertical covariant gradient of a scalar field.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:) ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field(:,:,:) ! z-component of resulting vector field
    type(t_grid), intent(in)    :: grid                ! the grid properties
    
    ! local variables
    integer :: jl ! spatial indices

    ! calculating the vertical gradient in the inner levels
    !$omp parallel do private(jl)
    do jl=2,n_layers
      result_field(:,:,jl) = (scalar_field(:,:,jl-1) - scalar_field(:,:,jl))/grid%dz(:,:,jl)
    enddo
    !$omp end parallel do

  end subroutine grad_vert
  
  subroutine grad_hor(scalar_field,result_field_x,result_field_y,result_field_z,grid)
  
    ! This subroutine computes the covariant gradient of a scalar field.
    ! result_field_z must be computed already
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
    real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
    real(wp),     intent(in)    :: result_field_z(:,:,:) ! z-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    
    ! local variables
    integer :: ji,jk,jl ! spatial indices
    
    ! computing the horizontal covariant gradient
    call grad_hor_cov(scalar_field,result_field_x,result_field_y,grid)

    ! correction for terrain
    
    ! correction to the x-component
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx+1
        do jl=n_flat_layers+1,n_layers
          result_field_x(ji,jk,jl) = result_field_x(ji,jk,jl) &
          - grid%slope_x(ji,jk,jl)*remap_ver2hor_x(result_field_z,grid,ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! correction to the y-component
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny+1
      do jk=1,nx
        do jl=n_flat_layers+1,n_layers
          result_field_y(ji,jk,jl) = result_field_y(ji,jk,jl) &
          - grid%slope_y(ji,jk,jl)*remap_ver2hor_y(result_field_z,grid,ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine grad_hor

end module mo_gradient_operators









