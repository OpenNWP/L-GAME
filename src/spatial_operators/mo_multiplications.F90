! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_multiplications
  
  ! This module is a collection of various multiplications of vector and/or scalar fields.
  
  use mo_definitions, only: wp
  use mo_run_nml,     only: ny,nx,n_layers
  use mo_bc_nml,      only: lperiodic
  
  implicit none
  
  contains

  subroutine scalar_times_vector_h(scalar_field,in_vector_x,in_vector_y,result_field_x,result_field_y)
    
    ! This subroutine multiplies a scalar with a vector field at horizontal points.
    
    real(wp), intent(in)  :: scalar_field(:,:,:)   ! input scalar field
    real(wp), intent(in)  :: in_vector_x(:,:,:)    ! input vector field, x-component
    real(wp), intent(in)  :: in_vector_y(:,:,:)    ! input vector field, y-component
    real(wp), intent(out) :: result_field_x(:,:,:) ! output vector field, x-component
    real(wp), intent(out) :: result_field_y(:,:,:) ! output vector field, y-component
  
    ! local variables
    integer :: ji,jk ! horizontal indices
    
    ! inner domain
    ! x
    !$omp parallel do private(jk)
    do jk=2,nx
      result_field_x(:,jk,:) = 0.5_wp*(scalar_field(:,jk-1,:) + scalar_field(:,jk,:))*in_vector_x(:,jk,:)
    enddo
    !$omp end parallel do
    
    ! y
    !$omp parallel do private(ji)
    do ji=2,ny
      result_field_y(ji,:,:) = 0.5_wp*(scalar_field(ji-1,:,:) + scalar_field(ji,:,:))*in_vector_y(ji,:,:)
    enddo
    !$omp end parallel do

    ! periodic boundary conditions
    if (lperiodic) then
      !$omp parallel workshare
      result_field_x(:,1,:) = 0.5_wp*(scalar_field(:,1,:) + scalar_field(:,nx,:))*in_vector_x(:,1,:)
      !$omp end parallel workshare
      
      !$omp parallel workshare
      result_field_x(:,nx+1,:) = result_field_x(:,1,:)
      !$omp end parallel workshare
      
      !$omp parallel workshare
      result_field_y(1,:,:) = 0.5_wp*(scalar_field(1,:,:) + scalar_field(ny,:,:))*in_vector_y(1,:,:)
      !$omp end parallel workshare
      
      !$omp parallel workshare
      result_field_y(ny+1,:,:) = result_field_y(1,:,:)
      !$omp end parallel workshare
    endif
    
  end subroutine scalar_times_vector_h

  subroutine scalar_times_vector_h2(scalar_field,vector_x,vector_y)
    
    ! This subroutine multiplies a scalar with a vector field at horizontal points and writes the result to the vector field.
    
    real(wp), intent(in)    :: scalar_field(:,:,:) ! input scalar field
    real(wp), intent(inout) :: vector_x(:,:,:)     ! vector field, x-component
    real(wp), intent(inout) :: vector_y(:,:,:)     ! vector field, y-component
  
    ! local variables
    integer :: ji,jk ! horizontal indices
    
    ! inner domain
    ! x
    !$omp parallel do private(jk)
    do jk=2,nx
      vector_x(:,jk,:) = 0.5_wp*(scalar_field(:,jk-1,:) + scalar_field(:,jk,:))*vector_x(:,jk,:)
    enddo
    !$omp end parallel do
    
    ! y
    !$omp parallel do private(ji)
    do ji=2,ny
      vector_y(ji,:,:) = 0.5_wp*(scalar_field(ji-1,:,:) + scalar_field(ji,:,:))*vector_y(ji,:,:)
    enddo
    !$omp end parallel do

    ! periodic boundary conditions
    if (lperiodic) then
      !$omp parallel workshare
      vector_x(:,1,:) = 0.5_wp*(scalar_field(:,1,:) + scalar_field(:,nx,:))*vector_x(:,1,:)
      !$omp end parallel workshare
      
      !$omp parallel workshare
      vector_x(:,nx+1,:) = vector_x(:,1,:)
      !$omp end parallel workshare
      
      !$omp parallel workshare
      vector_y(1,:,:) = 0.5_wp*(scalar_field(1,:,:) + scalar_field(ny,:,:))*vector_y(1,:,:)
      !$omp end parallel workshare
      
      !$omp parallel workshare
      vector_y(ny+1,:,:) = vector_y(1,:,:)
      !$omp end parallel workshare
    endif
    
  end subroutine scalar_times_vector_h2
  
  subroutine scalar_times_vector_h_upstream(scalar_field,in_vector_x,in_vector_y,result_field_x,result_field_y)
    
    ! This subroutine multiplies a scalar with a vector field at horizontal points.
    
    real(wp), intent(in)  :: scalar_field(:,:,:)   ! input scalar field
    real(wp), intent(in)  :: in_vector_x(:,:,:)    ! input vector field, x-component
    real(wp), intent(in)  :: in_vector_y(:,:,:)    ! input vector field, y-component
    real(wp), intent(out) :: result_field_x(:,:,:) ! output vector field, x-component
    real(wp), intent(out) :: result_field_y(:,:,:) ! output vector field, y-component
    
    ! local variables
    integer :: ji ! horizontal index
    integer :: jk ! horizontal index
    integer :: jl ! layer index
    
    ! inner domain
    ! x
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do jk=2,nx
        do ji=1,ny
          if (in_vector_x(ji,jk,jl)>=0._wp) then
            result_field_x(ji,jk,jl) = scalar_field(ji,jk-1,jl)*in_vector_x(ji,jk,jl)
          else
            result_field_x(ji,jk,jl) = scalar_field(ji,jk,jl)*in_vector_x(ji,jk,jl)
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! y
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do jk=1,nx
        do ji=2,ny
          if (in_vector_y(ji,jk,jl)>=0._wp) then
            result_field_y(ji,jk,jl) = scalar_field(ji,jk,jl)*in_vector_y(ji,jk,jl)
          else
            result_field_y(ji,jk,jl) = scalar_field(ji-1,jk,jl)*in_vector_y(ji,jk,jl)
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do

    ! periodic boundary conditions
    if (lperiodic) then
      !$omp parallel do private(ji,jl)
      do ji=1,ny
        do jl=1,n_layers
          if (in_vector_x(ji,1,jl)>=0._wp) then
            result_field_x(ji,1,jl) = scalar_field(ji,nx,jl)*in_vector_x(ji,1,jl)
          else
            result_field_x(ji,1,jl) = scalar_field(ji,1,jl)*in_vector_x(ji,1,jl)
          endif
        enddo
      enddo
      !$omp end parallel do
      
      !$omp parallel workshare
      result_field_x(:,nx+1,:) = result_field_x(:,1,:)
      !$omp end parallel workshare
      
      !$omp parallel do private(jk,jl)
      do jk=1,nx
        do jl=1,n_layers
          if (in_vector_y(1,jk,jl)>=0._wp) then
            result_field_y(1,jk,jl) = scalar_field(1,jk,jl)*in_vector_y(1,jk,jl)
          else
            result_field_y(1,jk,jl) = scalar_field(ny,jk,jl)*in_vector_y(1,jk,jl)
          endif
        enddo
      enddo
      !$omp end parallel do
      
      !$omp parallel workshare
      result_field_y(ny+1,:,:) = result_field_y(1,:,:)
      !$omp end parallel workshare
    endif
    
  end subroutine scalar_times_vector_h_upstream
  
  subroutine scalar_times_vector_v(scalar_field,in_vector_z,result_field_z)
  
    ! This subroutine multiplies of an extended scalar with a  vector field at vertical points.
    
    real(wp), intent(in)  :: scalar_field(:,:,:)   ! input scalar field
    real(wp), intent(in)  :: in_vector_z(:,:,:)    ! input vector field, z-component
    real(wp), intent(out) :: result_field_z(:,:,:) ! output vector field, z-component
    
    ! local variables
    integer :: ji ! horizontal index
    integer :: jk ! horizontal index
    integer :: jl ! layer index
    
    !$omp parallel do private(ji,jk,jl)
    do jl=2,n_layers
      do jk=1,nx
        do ji=1,ny
          result_field_z(ji,jk,jl) = 0.5_wp*(scalar_field(ji,jk,jl-1) + scalar_field(ji,jk,jl))*in_vector_z(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine scalar_times_vector_v
  
  subroutine scalar_times_vector_v2(scalar_field,vector_z)
  
    ! This subroutine multiplies of an extended scalar with a vector field at vertical points and writes the result to the vector field.
    
    real(wp), intent(in)    :: scalar_field(:,:,:) ! input scalar field
    real(wp), intent(inout) :: vector_z(:,:,:)     ! vector field, z-component
    
    ! local variables
    integer :: ji ! horizontal index
    integer :: jk ! horizontal index
    integer :: jl ! layer index
    
    !$omp parallel do private(ji,jk,jl)
    do jl=2,n_layers
      do jk=1,nx
        do ji=1,ny
          vector_z(ji,jk,jl) = 0.5_wp*(scalar_field(ji,jk,jl-1) + scalar_field(ji,jk,jl))*vector_z(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine scalar_times_vector_v2

end module mo_multiplications









