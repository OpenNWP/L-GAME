! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module multiplications

  ! This module is a collection of various multiplications of vector and/or scalar fields.
  
  use definitions, only: t_grid,wp
  use run_nml,     only: nlins,ncols,nlays
  
  implicit none
  
  private
  
  public :: scalar_times_vector_h
  public :: scalar_times_vector
  public :: scalar_times_vector_for_gradient
  
  contains

  subroutine scalar_times_vector_h(scalar_field,in_vector_x,in_vector_y, &
                   result_vector_x,result_vector_y)
  
    ! Multiplication of a scalar with a vector field at horizontal points.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)
    real(wp),     intent(in)    :: in_vector_x(:,:,:)
    real(wp),     intent(in)    :: in_vector_y(:,:,:)
    real(wp),     intent(inout) :: result_vector_x(:,:,:)
    real(wp),     intent(inout) :: result_vector_y(:,:,:)
  
    ! local variables
    integer                     :: ji,jk ! loop indices
    
    do jk=1,ncols+1
      result_vector_x(:,jk,:) = 0.5_wp*(scalar_field(:,jk,:) + scalar_field(:,jk+1,:))*in_vector_x(:,jk,:)
    enddo
    
    do ji=1,nlins+1
      result_vector_y(ji,:,:) = 0.5_wp*(scalar_field(ji,:,:) + scalar_field(ji+1,:,:))*in_vector_y(ji,:,:)
    enddo
  
  end subroutine scalar_times_vector_h
  
  subroutine scalar_times_vector_h_for_gradient(scalar_field,in_vector_x,in_vector_y, &
                   result_vector_x,result_vector_y)
  
    ! Multiplication of a scalar with a vector field at horizontal points.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)
    real(wp),     intent(in)    :: in_vector_x(:,:,:)
    real(wp),     intent(in)    :: in_vector_y(:,:,:)
    real(wp),     intent(inout) :: result_vector_x(:,:,:)
    real(wp),     intent(inout) :: result_vector_y(:,:,:)
  
    ! local variables
    integer                     :: ji,jk ! loop indices
    
    do jk=1,ncols-1
      result_vector_x(:,jk,:) = 0.5_wp*(scalar_field(2:nlins+1,jk+1,:) + scalar_field(2:nlins+1,jk+2,:))*in_vector_x(:,jk,:)
    enddo
    
    do ji=1,nlins-1
      result_vector_y(ji,:,:) = 0.5_wp*(scalar_field(ji+1,2:ncols+1,:) + scalar_field(ji+2,2:ncols+1,:))*in_vector_y(ji,:,:)
    enddo
  
  end subroutine scalar_times_vector_h_for_gradient
  
  subroutine scalar_times_vector_v(scalar_field,in_vector_z,result_vector_z,grid)
  
    ! Multiplication of an extended scalar with an inner vector field at vertical points.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)
    real(wp),     intent(in)    :: in_vector_z(:,:,:)
    real(wp),     intent(inout) :: result_vector_z(:,:,:)
    type(t_grid), intent(in)    :: grid
  
    ! local variables
    integer                     :: jl ! loop index
    
    do jl=2,nlays
      result_vector_z(:,:,jl) = 0.5_wp*(scalar_field(:,:,jl-1) + scalar_field(:,:,jl))*in_vector_z(:,:,jl)
    enddo
  
  end subroutine scalar_times_vector_v
  
  subroutine scalar_times_vector_v_for_gradient(scalar_field,in_vector_z,result_vector_z,grid)
  
    ! Multiplication of an extended scalar with an inner vector field at vertical points.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)
    real(wp),     intent(in)    :: in_vector_z(:,:,:)
    real(wp),     intent(inout) :: result_vector_z(:,:,:)
    type(t_grid), intent(in)    :: grid
  
    ! local variables
    integer                     :: jl ! loop index
    
    do jl=2,nlays
      result_vector_z(:,:,jl) = 0.5_wp*(scalar_field(2:nlins+1,2:ncols+1,jl-1) + scalar_field(2:nlins+1,2:ncols+1,jl)) &
      *in_vector_z(:,:,jl)
    enddo
  
  end subroutine scalar_times_vector_v_for_gradient
  
  subroutine scalar_times_vector(scalar_field,in_vector_x,in_vector_y,in_vector_z, &
                   result_vector_x,result_vector_y,result_vector_z,grid)
  
    ! Multiplication of a scalar with a vector field.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)
    real(wp),     intent(in)    :: in_vector_x(:,:,:)
    real(wp),     intent(in)    :: in_vector_y(:,:,:)
    real(wp),     intent(in)    :: in_vector_z(:,:,:)
    real(wp),     intent(inout) :: result_vector_x(:,:,:)
    real(wp),     intent(inout) :: result_vector_y(:,:,:)
    real(wp),     intent(inout) :: result_vector_z(:,:,:)
    type(t_grid), intent(in)    :: grid
    
    call scalar_times_vector_h(scalar_field,in_vector_x,in_vector_y,result_vector_x,result_vector_y)
    call scalar_times_vector_v(scalar_field,in_vector_z,result_vector_z,grid)
  
  end subroutine scalar_times_vector

  subroutine scalar_times_vector_for_gradient(scalar_field,in_vector_x,in_vector_y,in_vector_z, &
                   result_vector_x,result_vector_y,result_vector_z,grid)
  
    ! Multiplication of an extended scalar with an inner vector field.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)
    real(wp),     intent(in)    :: in_vector_x(:,:,:)
    real(wp),     intent(in)    :: in_vector_y(:,:,:)
    real(wp),     intent(in)    :: in_vector_z(:,:,:)
    real(wp),     intent(inout) :: result_vector_x(:,:,:)
    real(wp),     intent(inout) :: result_vector_y(:,:,:)
    real(wp),     intent(inout) :: result_vector_z(:,:,:)
    type(t_grid), intent(in)    :: grid
    
    call scalar_times_vector_h_for_gradient(scalar_field,in_vector_x,in_vector_y,result_vector_x,result_vector_y)
    call scalar_times_vector_v_for_gradient(scalar_field,in_vector_z,result_vector_z,grid)
  
  end subroutine scalar_times_vector_for_gradient

end module multiplications









