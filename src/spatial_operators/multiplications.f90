! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module multiplications

  ! This module is a collection of various multiplications of vector and/or scalar fields.
  
  use definitions, only: wp
  use run_nml,     only: nlins,ncols,nlays
  use bc_nml,      only: lperiodic
  
  implicit none
  
  private
  
  public :: scalar_times_scalar
  public :: scalar_times_vector_h
  public :: scalar_times_vector_h_upstream
  public :: scalar_times_vector
  public :: scalar_times_vector_v
  
  contains

  subroutine scalar_times_scalar(scalar_field_1,scalar_field_2,scalar_field_out)
  
    ! This subroutine multiplies of a scalar with another scalar field.
    
    ! input arguments and the result
    real(wp), intent(in)    :: scalar_field_1(:,:,:)   ! input field 1
    real(wp), intent(in)    :: scalar_field_2(:,:,:)   ! input field 2
    real(wp), intent(inout) :: scalar_field_out(:,:,:) ! the result
    
    ! actual calculation
    !$OMP PARALLEL
    !$OMP WORKSHARE
    scalar_field_out = scalar_field_1*scalar_field_2
    !$OMP END WORKSHARE
    !$OMP END PARALLEL
  
  end subroutine scalar_times_scalar

  subroutine scalar_times_vector_h(scalar_field,in_vector_x,in_vector_y, &
                   result_field_x,result_field_y)
  
    ! This subroutine multiplies of a scalar with a vector field at horizontal points.
    
    real(wp), intent(in)    :: scalar_field(:,:,:)   ! input scalar field
    real(wp), intent(in)    :: in_vector_x(:,:,:)    ! input vector field, x-component
    real(wp), intent(in)    :: in_vector_y(:,:,:)    ! input vector field, y-component
    real(wp), intent(inout) :: result_field_x(:,:,:) ! output vector field, x-component
    real(wp), intent(inout) :: result_field_y(:,:,:) ! output vector field, y-component
  
    ! local variables
    integer :: ji,jk ! loop indices
    
    ! inner domain
    ! x
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jk)
    do jk=2,ncols
      result_field_x(:,jk,:) = 0.5_wp*(scalar_field(:,jk-1,:) + scalar_field(:,jk,:))*in_vector_x(:,jk,:)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! y
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji)
    do ji=2,nlins
      result_field_y(ji,:,:) = 0.5_wp*(scalar_field(ji-1,:,:) + scalar_field(ji,:,:))*in_vector_y(ji,:,:)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! periodic boundary conditions
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_x(:,1,:) = 0.5_wp*(scalar_field(:,1,:) + scalar_field(:,ncols,:))*in_vector_x(:,1,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_x(:,ncols+1,:) = result_field_x(:,1,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_y(1,:,:) = 0.5_wp*(scalar_field(1,:,:) + scalar_field(nlins,:,:))*in_vector_y(1,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_y(nlins+1,:,:) = result_field_y(1,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
  
  end subroutine scalar_times_vector_h
  
  subroutine scalar_times_vector_h_upstream(scalar_field,in_vector_x,in_vector_y, &
                   result_field_x,result_field_y)
  
    ! This subroutine multiplies of a scalar with a vector field at horizontal points.
    
    real(wp), intent(in)    :: scalar_field(:,:,:)   ! input scalar field
    real(wp), intent(in)    :: in_vector_x(:,:,:)    ! input vector field, x-component
    real(wp), intent(in)    :: in_vector_y(:,:,:)    ! input vector field, y-component
    real(wp), intent(inout) :: result_field_x(:,:,:) ! output vector field, x-component
    real(wp), intent(inout) :: result_field_y(:,:,:) ! output vector field, y-component
  
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! inner domain
    ! x
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=2,ncols
        do jl=1,nlays
          result_field_x(ji,jk,jl) = 0.5_wp*(scalar_field(ji,jk-1,jl) + scalar_field(ji,jk,jl))*in_vector_x(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! y
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=2,nlins
      do jk=1,ncols
        do jl=1,nlays
          result_field_y(ji,jk,jl) = 0.5_wp*(scalar_field(ji-1,jk,jl) + scalar_field(ji,jk,jl))*in_vector_y(ji,jk,jl)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    ! periodic boundary conditions
    if (lperiodic) then
      !$OMP PARALLEL
      !$OMP DO PRIVATE(ji,jl)
      do ji=1,nlins
        do jl=1,nlays
          result_field_x(ji,1,jl) = 0.5_wp*(scalar_field(ji,1,jl) + scalar_field(ji,ncols,jl))*in_vector_x(ji,1,jl)
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_x(:,ncols+1,:) = result_field_x(:,1,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
      
      !$OMP PARALLEL
      !$OMP DO PRIVATE(jk,jl)
      do jk=1,ncols
        do jl=1,nlays
          result_field_y(1,jk,jl) = 0.5_wp*(scalar_field(1,jk,jl) + scalar_field(nlins,jk,jl))*in_vector_y(1,jk,jl)
        enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
      
      !$OMP PARALLEL
      !$OMP WORKSHARE
      result_field_y(nlins+1,:,:) = result_field_y(1,:,:)
      !$OMP END WORKSHARE
      !$OMP END PARALLEL
    endif
  
  end subroutine scalar_times_vector_h_upstream
  
  subroutine scalar_times_vector_v(scalar_field,in_vector_z,result_field_z)
  
    ! This subroutine multiplies of an extended scalar with an inner vector field at vertical points.
    
    real(wp), intent(in)    :: scalar_field(:,:,:)   ! input scalar field
    real(wp), intent(in)    :: in_vector_z(:,:,:)    ! input vector field, z-component
    real(wp), intent(inout) :: result_field_z(:,:,:) ! output vector field, z-component
  
    ! local variables
    integer :: jl ! loop index
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(jl)
    do jl=2,nlays
      result_field_z(:,:,jl) = 0.5_wp*(scalar_field(:,:,jl-1) + scalar_field(:,:,jl))*in_vector_z(:,:,jl)
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine scalar_times_vector_v
  
  subroutine scalar_times_vector(scalar_field,in_vector_x,in_vector_y,in_vector_z, &
                   result_field_x,result_field_y,result_field_z)
  
    ! This subroutine of a scalar with a vector field.
    
    real(wp), intent(in)    :: scalar_field(:,:,:)   ! input scalar field
    real(wp), intent(in)    :: in_vector_x(:,:,:)    ! input vector field, x-component
    real(wp), intent(in)    :: in_vector_y(:,:,:)    ! input vector field, y-component
    real(wp), intent(in)    :: in_vector_z(:,:,:)    ! input vector field, z-component
    real(wp), intent(inout) :: result_field_x(:,:,:) ! output vector field, x-component
    real(wp), intent(inout) :: result_field_y(:,:,:) ! output vector field, y-component
    real(wp), intent(inout) :: result_field_z(:,:,:) ! output vector field, z-component
    
    ! this subroutine just calls the horizontal and vertical subroutines one after the other
    call scalar_times_vector_h(scalar_field,in_vector_x,in_vector_y,result_field_x,result_field_y)
    call scalar_times_vector_v(scalar_field,in_vector_z,result_field_z)
  
  end subroutine scalar_times_vector

end module multiplications









