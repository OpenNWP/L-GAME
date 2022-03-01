! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module inner_product

  ! The calculation of the inner product is executed in this module.

  use run_nml,     only: nlins,ncols,nlays
  use definitions, only: t_grid,wp
  
  implicit none
  
  private
  
  public :: inner
  
  contains

  subroutine inner(u_vector_1,v_vector_1,w_vector_1,u_vector_2,v_vector_2,w_vector_2,output_scalar,grid)
  
    ! This subroutine calculates the inner product of two vector fields.
    
    ! input arguments and output
    real(wp), intent(in)     :: u_vector_1(:,:,:)    ! vectorfield 1 in x-direction
    real(wp), intent(in)     :: v_vector_1(:,:,:)    ! vectorfield 1 in y-direction
    real(wp), intent(in)     :: w_vector_1(:,:,:)    ! vectorfield 1 in z-direction
    real(wp), intent(in)     :: u_vector_2(:,:,:)    ! vectorfield 2 in x-direction
    real(wp), intent(in)     :: v_vector_2(:,:,:)    ! vectorfield 2 in y-direction
    real(wp), intent(in)     :: w_vector_2(:,:,:)    ! vectorfield 2 in z-direction
    real(wp), intent(inout)  :: output_scalar(:,:,:) ! result
    type(t_grid), intent(in) :: grid                 ! grid properties
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          output_scalar(ji,jk,jl) = &
          grid%inner_product_weights(ji,jk,jl,1)*u_vector_1(ji,jk+1,jl)*u_vector_2(ji,jk+1,jl) &
          + grid%inner_product_weights(ji,jk,jl,2)*v_vector_1(ji,jk,jl)*v_vector_2(ji,jk,jl) &
          + grid%inner_product_weights(ji,jk,jl,3)*u_vector_1(ji,jk,jl)*u_vector_2(ji,jk,jl) &
          + grid%inner_product_weights(ji,jk,jl,4)*v_vector_1(ji+1,jk,jl)*v_vector_2(ji+1,jk,jl) &
          + grid%inner_product_weights(ji,jk,jl,5)*w_vector_1(ji,jk,jl)*w_vector_2(ji,jk,jl) &
          + grid%inner_product_weights(ji,jk,jl,6)*w_vector_1(ji,jk,jl+1)*w_vector_2(ji,jk,jl+1)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine inner

end module inner_product






