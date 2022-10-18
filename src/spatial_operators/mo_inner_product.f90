! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_inner_product

  ! The calculation of the inner product is executed in this module.

  use mo_run_nml,     only: ny,nx,n_layers
  use mo_definitions, only: t_grid,wp
  
  implicit none
  
  contains

  subroutine inner_product(u_vector_1,v_vector_1,w_vector_1,u_vector_2,v_vector_2,w_vector_2,output_scalar,grid)
  
    ! This subroutine calculates the inner product of two vector fields.
    
    ! input arguments and output
    real(wp),     intent(in)    :: u_vector_1(:,:,:)    ! vectorfield 1 in x-direction
    real(wp),     intent(in)    :: v_vector_1(:,:,:)    ! vectorfield 1 in y-direction
    real(wp),     intent(in)    :: w_vector_1(:,:,:)    ! vectorfield 1 in z-direction
    real(wp),     intent(in)    :: u_vector_2(:,:,:)    ! vectorfield 2 in x-direction
    real(wp),     intent(in)    :: v_vector_2(:,:,:)    ! vectorfield 2 in y-direction
    real(wp),     intent(in)    :: w_vector_2(:,:,:)    ! vectorfield 2 in z-direction
    real(wp),     intent(inout) :: output_scalar(:,:,:) ! result
    type(t_grid), intent(in)    :: grid                 ! grid properties
    
    ! local variables
    integer :: ji,jk,jl ! spatial indices
    
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do jk=1,nx
        do ji=1,ny
          output_scalar(ji,jk,jl) = &
          grid%inner_product_weights(1,ji,jk,jl)*u_vector_1(ji,jk+1,jl)*u_vector_2(ji,jk+1,jl) &
          + grid%inner_product_weights(2,ji,jk,jl)*v_vector_1(ji,jk,jl)*v_vector_2(ji,jk,jl) &
          + grid%inner_product_weights(3,ji,jk,jl)*u_vector_1(ji,jk,jl)*u_vector_2(ji,jk,jl) &
          + grid%inner_product_weights(4,ji,jk,jl)*v_vector_1(ji+1,jk,jl)*v_vector_2(ji+1,jk,jl) &
          + grid%inner_product_weights(5,ji,jk,jl)*w_vector_1(ji,jk,jl)*w_vector_2(ji,jk,jl) &
          + grid%inner_product_weights(6,ji,jk,jl)*w_vector_1(ji,jk,jl+1)*w_vector_2(ji,jk,jl+1)
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine inner_product

end module mo_inner_product






