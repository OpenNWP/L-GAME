! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_inner_product

  ! The calculation of the inner product is executed in this module.

  use mo_run_nml,              only: ny,nx,n_layers,n_levels,dy,dx
  use mo_definitions,          only: t_grid,t_state,t_diag,wp
  use mo_gradient_operators,   only: grad_hor,grad_vert
  use mo_divergence_operators, only: div_h
  use mo_bc_nml,               only: lperiodic
  
  implicit none
  
  contains

  subroutine inner_product(u_vector_1,v_vector_1,w_vector_1,u_vector_2,v_vector_2,w_vector_2,output_scalar,grid)
    
    ! This subroutine calculates the inner product of two vector fields.
    
    real(wp),     intent(in)  :: u_vector_1(:,:,:)    ! vectorfield 1 in x-direction
    real(wp),     intent(in)  :: v_vector_1(:,:,:)    ! vectorfield 1 in y-direction
    real(wp),     intent(in)  :: w_vector_1(:,:,:)    ! vectorfield 1 in z-direction
    real(wp),     intent(in)  :: u_vector_2(:,:,:)    ! vectorfield 2 in x-direction
    real(wp),     intent(in)  :: v_vector_2(:,:,:)    ! vectorfield 2 in y-direction
    real(wp),     intent(in)  :: w_vector_2(:,:,:)    ! vectorfield 2 in z-direction
    real(wp),     intent(out) :: output_scalar(:,:,:) ! result
    type(t_grid), intent(in)  :: grid                 ! grid properties
    
    ! local variables
    integer :: ji ! horizontal index
    integer :: jk ! horizontal index
    integer :: jl ! layer index
    
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
  
  subroutine theta_v_adv_3rd_order(state,diag,grid)
    
    ! This subroutine computes the virtual potential temperature at the edges for third-order upwind advection.
    ! It is assumed that diag%scalar_placeholder holds the full virtual potential temperature.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid properties
    
    ! local variables
    integer :: ji ! horizontal index
    integer :: jk ! horizontal index
    integer :: jl ! layer index
    
    call grad_vert(diag%scalar_placeholder,diag%w_placeholder,grid)
    call grad_hor(diag%scalar_placeholder,diag%u_placeholder,diag%v_placeholder,diag%w_placeholder,grid)
    call div_h(diag%u_placeholder,diag%v_placeholder,diag%flux_density_div,grid)
    
    ! inner domain
    ! x
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do jk=2,nx
        do ji=1,ny
          diag%theta_v_u(ji,jk,jl) = 0.5_wp*(diag%scalar_placeholder(ji,jk-1,jl) + diag%scalar_placeholder(ji,jk,jl))
          if (state%wind_u(ji,jk,jl)>=0._wp) then
            diag%theta_v_u(ji,jk,jl) = diag%theta_v_u(ji,jk,jl) - 1._wp/6._wp*dy*dx*diag%flux_density_div(ji,jk-1,jl)
          else
            diag%theta_v_u(ji,jk,jl) = diag%theta_v_u(ji,jk,jl) - 1._wp/6._wp*dy*dx*diag%flux_density_div(ji,jk,jl)
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
          diag%theta_v_v(ji,jk,jl) = 0.5_wp*(diag%scalar_placeholder(ji,jk,jl) + diag%scalar_placeholder(ji-1,jk,jl))
          if (state%wind_v(ji,jk,jl)>=0._wp) then
            diag%theta_v_v(ji,jk,jl) = diag%theta_v_v(ji,jk,jl) - 1._wp/6._wp*dy*dx*diag%flux_density_div(ji,jk,jl)
          else
            diag%theta_v_v(ji,jk,jl) = diag%theta_v_v(ji,jk,jl) - 1._wp/6._wp*dy*dx*diag%flux_density_div(ji-1,jk,jl)
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
          diag%theta_v_u(ji,1,jl) = 0.5_wp*(diag%scalar_placeholder(ji,nx,jl) + diag%scalar_placeholder(ji,1,jl))
          if (state%wind_u(ji,1,jl)>=0._wp) then
            diag%theta_v_u(ji,1,jl) = diag%theta_v_u(ji,1,jl) - 1._wp/6._wp*dy*dx*diag%flux_density_div(ji,nx,jl)
          else
            diag%theta_v_u(ji,1,jl) = diag%theta_v_u(ji,1,jl) - 1._wp/6._wp*dy*dx*diag%flux_density_div(ji,1,jl)
          endif
        enddo
      enddo
      !$omp end parallel do
      
      !$omp parallel workshare
      diag%theta_v_u(:,nx+1,:) = diag%theta_v_u(:,1,:)
      !$omp end parallel workshare
      
      !$omp parallel do private(jk,jl)
      do jk=1,nx
        do jl=1,n_layers
          diag%theta_v_v(1,jk,jl) = 0.5_wp*(diag%scalar_placeholder(1,jk,jl) + diag%scalar_placeholder(ny,jk,jl))
          if (state%wind_v(1,jk,jl)>=0._wp) then
            diag%theta_v_v(1,jk,jl) = diag%theta_v_v(1,jk,jl) - 1._wp/6._wp*dy*dx*diag%flux_density_div(1,jk,jl)
          else
            diag%theta_v_v(1,jk,jl) = diag%theta_v_v(1,jk,jl) - 1._wp/6._wp*dy*dx*diag%flux_density_div(ny,jk,jl)
          endif
        enddo
      enddo
      !$omp end parallel do
      
      !$omp parallel workshare
      diag%theta_v_v(ny+1,:,:) = diag%theta_v_v(1,:,:)
      !$omp end parallel workshare
    endif
    
  end subroutine theta_v_adv_3rd_order
  
  subroutine w_free_slip(state,grid)
    
    ! This subroutine calculates the vertical velocity at the surface under kinematic boundary conditions.
    
    type(t_state), intent(inout) :: state ! state variables
    type(t_grid),  intent(in)    :: grid  ! grid properties
    
    ! local variables
    integer :: ji ! horizontal index
    integer :: jk ! horizontal index
    
    !$omp parallel do private(ji,jk)
    do jk=1,nx
      do ji=1,ny
        state%wind_w(ji,jk,n_levels) = &
        grid%inner_product_weights(1,ji,jk,n_layers)*state%wind_u(ji,jk+1,n_layers)*grid%slope_x(ji,jk+1,n_layers) &
        + grid%inner_product_weights(2,ji,jk,n_layers)*state%wind_v(ji,jk,n_layers)*grid%slope_y(ji,jk,n_layers) &
        + grid%inner_product_weights(3,ji,jk,n_layers)*state%wind_u(ji,jk,n_layers)*grid%slope_x(ji,jk,n_layers) &
        + grid%inner_product_weights(4,ji,jk,n_layers)*state%wind_v(ji+1,jk,n_layers)*grid%slope_y(ji+1,jk,n_layers)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine w_free_slip

end module mo_inner_product






