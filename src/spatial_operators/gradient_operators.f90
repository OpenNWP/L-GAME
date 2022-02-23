! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

! This file is a collection of gradient operators.

module gradient_operators

  use definitions, only: t_grid,wp
  use run_nml,     only: nlins,ncols,nlays,toa
  use averaging,   only: hor_cov_to_con
    
  implicit none
  
  private
  
  public :: grad_hor_cov
  public :: grad
  public :: grad_hor
  
  contains
  
  subroutine grad_hor_cov(scalar_field,result_field_x,result_field_y,grid)

    ! This subroutine computes the horizontal covariant gradient of a scalar field.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
    real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    ! local variables
    integer                     :: ji,jk                 ! loop variables

    ! inner domain
    ! calculating the x-component of the gradient
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins
      do jk=2,ncols
        result_field_x(ji,jk,:) = (scalar_field(ji,jk,:) - scalar_field(ji,jk-1,:))/grid%dx(ji,jk,:)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! calculating the y-component of the gradient
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=2,nlins
      do jk=1,ncols
        result_field_y(ji,jk,:) = (scalar_field(ji,jk,:) - scalar_field(ji-1,jk,:))/grid%dy(ji,jk,:)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! boundaries

  end subroutine grad_hor_cov
  
  subroutine grad_vert_cov(scalar_field,result_field,grid)

    ! This subroutine computes the vertical covariant gradient of a scalar field.
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field(:,:,:)   ! z-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    ! local variables
    integer                     :: jl                    ! loop variables

    ! calculating the vertical gradient in the inner levels
    do jl=2,nlays
      result_field(:,:,jl) = (scalar_field(:,:,jl-1) - scalar_field(:,:,jl))/grid%dz(:,:,jl)
    enddo
    result_field(:,:,1) = 0._wp
    result_field(:,:,nlays+1) = 0._wp

  end subroutine grad_vert_cov
  
  subroutine grad_cov(scalar_field,result_field_x,result_field_y,result_field_z,grid)
  
    ! This subroutine computes the covariant gradient of a scalar field.
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
    real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
    real(wp),     intent(inout) :: result_field_z(:,:,:) ! z-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    
    call grad_hor_cov(scalar_field,result_field_x,result_field_y,grid)
    call grad_vert_cov(scalar_field,result_field_z,grid)
  
  end subroutine grad_cov
  
  subroutine grad(scalar_field,result_field_x,result_field_y,result_field_z,grid)
  
    ! This subroutine computes the covariant gradient of a scalar field.
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
    real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
    real(wp),     intent(inout) :: result_field_z(:,:,:) ! z-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    
    ! covariant gradient
    call grad_cov(scalar_field,result_field_x,result_field_y,result_field_z,grid)
    ! correction for terrain
    call hor_cov_to_con(result_field_x,result_field_y,result_field_z,grid)
  
  end subroutine grad
  
  subroutine grad_hor(scalar_field,result_field_x,result_field_y,result_field_z,grid)
  
    ! This subroutine computes the covariant gradient of a scalar field.
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
    real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
    real(wp),     intent(inout) :: result_field_z(:,:,:) ! z-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    
    ! calling the gradient
    call grad(scalar_field,result_field_x,result_field_y,result_field_z,grid)
    ! setting the vertical component to zero
    result_field_z(:,:,:) = 0._wp
  
  end subroutine grad_hor

end module gradient_operators









