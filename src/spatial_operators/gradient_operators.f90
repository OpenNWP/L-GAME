! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module gradient_operators

  ! This module is a collection of gradient operators.

  use definitions, only: t_grid,wp
  use run_nml,     only: nlins,ncols,nlays,toa
  use averaging,   only: hor_cov_to_con
  use bc_nml,      only: lperiodic
    
  implicit none
  
  private
  
  public :: grad_hor_cov
  public :: grad
  public :: grad_hor
  
  contains
  
  subroutine grad_hor_cov(scalar_field,result_field_x,result_field_y,grid)

    ! This subroutine computes the horizontal covariant gradient of a scalar field.
    
    ! input arguments and output
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
    real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    
    ! local variables
    integer :: ji,jk ! loop variables

    ! inner domain
    ! calculating the x-component of the gradient
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=2,ncols
      result_field_x(:,jk,:) = (scalar_field(:,jk,:) - scalar_field(:,jk-1,:))/grid%dx(:,jk,:)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! calculating the y-component of the gradient
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji)
    do ji=2,nlins
      result_field_y(ji,:,:) = (scalar_field(ji-1,:,:) - scalar_field(ji,:,:))/grid%dy(ji,:,:)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! periodic boundary conditions
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_x(:,1,:) = (scalar_field(:,1,:) - scalar_field(:,ncols,:))/grid%dx(:,1,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_x(:,ncols+1,:) = result_field_x(:,1,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_y(1,:,:) = (scalar_field(1,:,:) - scalar_field(nlins,:,:))/grid%dy(1,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_y(nlins+1,:,:) = result_field_y(1,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif

  end subroutine grad_hor_cov
  
  subroutine grad_vert_cov(scalar_field,result_field,grid)

    ! This subroutine computes the vertical covariant gradient of a scalar field.
    
    ! input arguments and output
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field(:,:,:)   ! z-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    
    ! local variables
    integer :: jl ! loop variables

    ! calculating the vertical gradient in the inner levels
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl)
    do jl=2,nlays
      result_field(:,:,jl) = (scalar_field(:,:,jl-1) - scalar_field(:,:,jl))/grid%dz(:,:,jl)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    !$OMP PARALLEL
    !$OMP WORKSHARE
    result_field(:,:,1) = 0._wp
    result_field(:,:,nlays+1) = 0._wp
    !$OMP END WORKSHARE
    !$OMP END PARALLEL

  end subroutine grad_vert_cov
  
  subroutine grad_cov(scalar_field,result_field_x,result_field_y,result_field_z,grid)
  
    ! This subroutine computes the covariant gradient of a scalar field.
    
    ! input arguments and output
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
    
    ! input arguments and output
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
    
    ! input arguments and output
    real(wp),     intent(in)    :: scalar_field(:,:,:)   ! scalar field of which to calculate the gradient
    real(wp),     intent(inout) :: result_field_x(:,:,:) ! x-component of resulting vector field
    real(wp),     intent(inout) :: result_field_y(:,:,:) ! y-component of resulting vector field
    real(wp),     intent(inout) :: result_field_z(:,:,:) ! z-component of resulting vector field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    
    ! calling the gradient
    call grad(scalar_field,result_field_x,result_field_y,result_field_z,grid)
    ! setting the vertical component to zero
    !$OMP PARALLEL
    !$OMP WORKSHARE
    result_field_z = 0._wp
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
  
  end subroutine grad_hor

end module gradient_operators









