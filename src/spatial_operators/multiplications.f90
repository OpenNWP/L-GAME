! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module multiplications

  ! This module is a collection of various multiplications of vector and/or scalar fields.
  
  use definitions, only: wp
  use run_nml,     only: nlins,ncols,nlays
  use bc_nml,      only: lperiodic
  
  implicit none
  
  private
  
  public :: scalar_times_vector_h
  public :: scalar_times_vector
  
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
    
    ! inner domain
    ! x
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=2,ncols
      result_vector_x(:,jk,:) = 0.5_wp*(scalar_field(:,jk-1,:) + scalar_field(:,jk,:))*in_vector_x(:,jk,:)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! y
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji)
    do ji=2,nlins
      result_vector_y(ji,:,:) = 0.5_wp*(scalar_field(ji-1,:,:) + scalar_field(ji,:,:))*in_vector_y(ji,:,:)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
    ! boundaries
  
  end subroutine scalar_times_vector_h
  
  subroutine scalar_times_vector_v(scalar_field,in_vector_z,result_vector_z)
  
    ! Multiplication of an extended scalar with an inner vector field at vertical points.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)
    real(wp),     intent(in)    :: in_vector_z(:,:,:)
    real(wp),     intent(inout) :: result_vector_z(:,:,:)
  
    ! local variables
    integer                     :: jl ! loop index
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl)
    do jl=2,nlays
      result_vector_z(:,:,jl) = 0.5_wp*(scalar_field(:,:,jl-1) + scalar_field(:,:,jl))*in_vector_z(:,:,jl)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    result_vector_z(:,:,1) = 0._wp
    result_vector_z(:,:,nlays+1) = 0._wp
  
  end subroutine scalar_times_vector_v
  
  subroutine scalar_times_vector(scalar_field,in_vector_x,in_vector_y,in_vector_z, &
                   result_vector_x,result_vector_y,result_vector_z)
  
    ! Multiplication of a scalar with a vector field.
    
    real(wp),     intent(in)    :: scalar_field(:,:,:)
    real(wp),     intent(in)    :: in_vector_x(:,:,:)
    real(wp),     intent(in)    :: in_vector_y(:,:,:)
    real(wp),     intent(in)    :: in_vector_z(:,:,:)
    real(wp),     intent(inout) :: result_vector_x(:,:,:)
    real(wp),     intent(inout) :: result_vector_y(:,:,:)
    real(wp),     intent(inout) :: result_vector_z(:,:,:)
    
    call scalar_times_vector_h(scalar_field,in_vector_x,in_vector_y,result_vector_x,result_vector_y)
    call scalar_times_vector_v(scalar_field,in_vector_z,result_vector_z)
  
  end subroutine scalar_times_vector

end module multiplications









