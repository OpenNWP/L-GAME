! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_multiplications
  
  ! This module is a collection of various multiplications of vector and/or scalar fields.
  
  use mo_definitions,          only: wp,t_state,t_diag,t_grid
  use mo_run_nml,              only: ny,nx,n_layers,dx,dy
  use mo_bc_nml,               only: lperiodic
  use mo_gradient_operators,   only: grad_hor,grad_vert
  use mo_divergence_operators, only: div_h
  
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
    integer :: ji ! horizontal index
    integer :: jk ! horizontal index
    
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
    integer :: ji ! horizontal index
    integer :: jk ! horizontal index
    
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
            diag%theta_v_u(ji,jk,jl) = diag%theta_v_u(ji,jk,jl) - 1._wp/8._wp*dy*dx*diag%flux_density_div(ji,jk-1,jl)
          else
            diag%theta_v_u(ji,jk,jl) = diag%theta_v_u(ji,jk,jl) - 1._wp/8._wp*dy*dx*diag%flux_density_div(ji,jk,jl)
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
            diag%theta_v_v(ji,jk,jl) = diag%theta_v_v(ji,jk,jl) - 1._wp/8._wp*dy*dx*diag%flux_density_div(ji,jk,jl)
          else
            diag%theta_v_v(ji,jk,jl) = diag%theta_v_v(ji,jk,jl) - 1._wp/8._wp*dy*dx*diag%flux_density_div(ji-1,jk,jl)
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
            diag%theta_v_u(ji,1,jl) = diag%theta_v_u(ji,1,jl) - 1._wp/8._wp*dy*dx*diag%flux_density_div(ji,nx,jl)
          else
            diag%theta_v_u(ji,1,jl) = diag%theta_v_u(ji,1,jl) - 1._wp/8._wp*dy*dx*diag%flux_density_div(ji,1,jl)
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
            diag%theta_v_v(1,jk,jl) = diag%theta_v_v(1,jk,jl) - 1._wp/8._wp*dy*dx*diag%flux_density_div(1,jk,jl)
          else
            diag%theta_v_v(1,jk,jl) = diag%theta_v_v(1,jk,jl) - 1._wp/8._wp*dy*dx*diag%flux_density_div(ny,jk,jl)
          endif
        enddo
      enddo
      !$omp end parallel do
      
      !$omp parallel workshare
      diag%theta_v_v(ny+1,:,:) = diag%theta_v_v(1,:,:)
      !$omp end parallel workshare
    endif
    
  end subroutine theta_v_adv_3rd_order
  
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














