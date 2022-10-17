! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_vorticity_flux

  ! This module computes the vorticity flux term.
  
  use mo_definitions, only: t_grid,t_diag,wp
  use mo_run_nml,     only: ny,nx,n_layers
  use mo_bc_nml,      only: lperiodic
  
  implicit none
  
  contains
  
  subroutine calc_vorticity_flux_term(diag,grid)
  
    ! This module computes the vorticity flux.
    
    type(t_diag), intent(inout) :: diag ! diagnostic quantities (the vorticity flux term is included here)
    type(t_grid), intent(in)    :: grid ! model grid
    
    ! local variables
    integer :: ji,jk,jl ! spatial indices
    
    ! horizontal velocity tendency due to vertical vorticity and horizontal wind (TRSK)
    ! u
    !$omp parallel do private(ji,jk)
    do ji=1,ny
      do jk=2,nx
        diag%pot_vort_tend_x(ji,jk,:) = &
        grid%trsk_weights_u(ji,1)*diag%v_placeholder(ji,jk-1,:)*0.25_wp* &
        (diag%eta_z(ji,jk-1,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji+1,jk,:)) &
        + grid%trsk_weights_u(ji,2)*diag%u_placeholder(ji,jk-1,:)*0.25_wp* &
        (diag%eta_z(ji,jk-1,:)+diag%eta_z(ji+1,jk-1,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji+1,jk,:)) &
        + grid%trsk_weights_u(ji,3)*diag%v_placeholder(ji+1,jk-1,:)*0.25_wp* &
        (diag%eta_z(ji+1,jk-1,:)+diag%eta_z(ji+1,jk,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji+1,jk,:)) &
        + grid%trsk_weights_u(ji,4)*diag%v_placeholder(ji+1,jk,:)*0.25_wp* &
        (diag%eta_z(ji+1,jk,:)+diag%eta_z(ji+1,jk+1,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji+1,jk,:)) &
        + grid%trsk_weights_u(ji,5)*diag%u_placeholder(ji,jk+1,:)*0.25_wp* &
        (diag%eta_z(ji+1,jk+1,:)+diag%eta_z(ji,jk+1,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji+1,jk,:)) &
        + grid%trsk_weights_u(ji,6)*diag%v_placeholder(ji,jk,:)*0.25_wp* &
        (diag%eta_z(ji,jk,:)+diag%eta_z(ji,jk+1,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji+1,jk,:))
      enddo
      
      ! boundary conditions
      if (lperiodic) then
        diag%pot_vort_tend_x(ji,1,:) = &
        grid%trsk_weights_u(ji,1)*diag%v_placeholder(ji,nx,:)*0.25_wp* &
        (diag%eta_z(ji,nx,:)+diag%eta_z(ji,1,:)+diag%eta_z(ji,1,:)+diag%eta_z(ji+1,1,:)) &
        + grid%trsk_weights_u(ji,2)*diag%u_placeholder(ji,nx,:)*0.25_wp* &
        (diag%eta_z(ji,nx,:)+diag%eta_z(ji+1,nx,:)+diag%eta_z(ji,1,:)+diag%eta_z(ji+1,1,:)) &
        + grid%trsk_weights_u(ji,3)*diag%v_placeholder(ji+1,nx,:)*0.25_wp* &
        (diag%eta_z(ji+1,nx,:)+diag%eta_z(ji+1,1,:)+diag%eta_z(ji,1,:)+diag%eta_z(ji+1,1,:)) &
        + grid%trsk_weights_u(ji,4)*diag%v_placeholder(ji+1,1,:)*0.25_wp* &
        (diag%eta_z(ji+1,1,:)+diag%eta_z(ji+1,2,:)+diag%eta_z(ji,1,:)+diag%eta_z(ji+1,1,:)) &
        + grid%trsk_weights_u(ji,5)*diag%u_placeholder(ji,2,:)*0.25_wp* &
        (diag%eta_z(ji+1,2,:)+diag%eta_z(ji,2,:)+diag%eta_z(ji,1,:)+diag%eta_z(ji+1,1,:)) &
        + grid%trsk_weights_u(ji,6)*diag%v_placeholder(ji,1,:)*0.25_wp* &
        (diag%eta_z(ji,1,:)+diag%eta_z(ji,2,:)+diag%eta_z(ji,1,:)+diag%eta_z(ji+1,1,:))
         diag%pot_vort_tend_x(ji,nx+1,:) = diag%pot_vort_tend_x(ji,1,:)
      endif
      
    enddo
    !$omp end parallel do
    ! v
    !$omp parallel do private(ji,jk)
    do jk=1,nx
      do ji=2,ny
        diag%pot_vort_tend_y(ji,jk,:) = &
        grid%trsk_weights_v(ji,1)*diag%u_placeholder(ji,jk,:)*0.25_wp* &
        (diag%eta_z(ji,jk,:)+diag%eta_z(ji+1,jk,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji,jk+1,:)) &
        + grid%trsk_weights_v(ji,2)*diag%u_placeholder(ji,jk+1,:)*0.25_wp* &
        (diag%eta_z(ji,jk+1,:)+diag%eta_z(ji+1,jk+1,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji,jk+1,:)) &
        + grid%trsk_weights_v(ji,3)*diag%u_placeholder(ji-1,jk+1,:)*0.25_wp* &
        (diag%eta_z(ji-1,jk+1,:)+diag%eta_z(ji,jk+1,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji,jk+1,:)) &
        + grid%trsk_weights_v(ji,4)*diag%u_placeholder(ji-1,jk,:)*0.25_wp* &
        (diag%eta_z(ji-1,jk,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji,jk,:)+diag%eta_z(ji,jk+1,:))
      enddo
      
      ! boundary conditions
      if (lperiodic) then
        diag%pot_vort_tend_y(1,jk,:) = &
        grid%trsk_weights_v(1,1)*diag%u_placeholder(1,jk,:)*0.25_wp* &
        (diag%eta_z(1,jk,:)+diag%eta_z(2,jk,:)+diag%eta_z(1,jk,:)+diag%eta_z(1,jk+1,:)) &
        + grid%trsk_weights_v(1,2)*diag%u_placeholder(1,jk+1,:)*0.25_wp* &
        (diag%eta_z(1,jk+1,:)+diag%eta_z(2,jk+1,:)+diag%eta_z(1,jk,:)+diag%eta_z(1,jk+1,:)) &
        + grid%trsk_weights_v(1,3)*diag%u_placeholder(ny,jk+1,:)*0.25_wp* &
        (diag%eta_z(ny,jk+1,:)+diag%eta_z(1,jk+1,:)+diag%eta_z(1,jk,:)+diag%eta_z(1,jk+1,:)) &
        + grid%trsk_weights_v(1,4)*diag%u_placeholder(ny,jk,:)*0.25_wp* &
        (diag%eta_z(ny,jk,:)+diag%eta_z(1,jk,:)+diag%eta_z(1,jk,:)+diag%eta_z(1,jk+1,:))
        diag%pot_vort_tend_y(ny+1,jk,:) = diag%pot_vort_tend_y(1,jk,:)
      endif
      
    enddo
    !$omp end parallel do
    
    ! horizontal velocity tendency due to horizontal vorticity and vertical wind
    ! u
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do ji=1,ny
        do jk=2,nx
          diag%pot_vort_tend_x(ji,jk,jl) = diag%pot_vort_tend_x(ji,jk,jl) &
          - 0.5_wp*grid%inner_product_weights(5,ji,jk-1,jl)*diag%w_placeholder(ji,jk-1,jl)*diag%eta_y(ji,jk,jl) &
          - 0.5_wp*grid%inner_product_weights(6,ji,jk-1,jl)*diag%w_placeholder(ji,jk-1,jl+1)*diag%eta_y(ji,jk,jl+1) &
          - 0.5_wp*grid%inner_product_weights(5,ji,jk,jl)*diag%w_placeholder(ji,jk,jl)*diag%eta_y(ji,jk,jl) &
          - 0.5_wp*grid%inner_product_weights(6,ji,jk,jl)*diag%w_placeholder(ji,jk,jl+1)*diag%eta_y(ji,jk,jl+1)
        enddo
        
        ! boundary conditions
        if (lperiodic) then
          diag%pot_vort_tend_x(ji,1,jl) = diag%pot_vort_tend_x(ji,1,jl) &
          - 0.5_wp*grid%inner_product_weights(5,ji,nx,jl)*diag%w_placeholder(ji,nx,jl)*diag%eta_y(ji,1,jl) &
          - 0.5_wp*grid%inner_product_weights(6,ji,nx,jl)*diag%w_placeholder(ji,nx,jl+1)*diag%eta_y(ji,1,jl+1) &
          - 0.5_wp*grid%inner_product_weights(5,ji,1,jl)*diag%w_placeholder(ji,1,jl)*diag%eta_y(ji,1,jl) &
          - 0.5_wp*grid%inner_product_weights(6,ji,1,jl)*diag%w_placeholder(ji,1,jl+1)*diag%eta_y(ji,1,jl+1)
          diag%pot_vort_tend_x(ji,nx+1,jl) = diag%pot_vort_tend_x(ji,1,jl)
        endif
        
      enddo
    enddo
    !$omp end parallel do
    ! v
    !$omp parallel do private(ji,jk,jl)
    do jl=1,n_layers
      do jk=1,nx
        do ji=2,ny
          diag%pot_vort_tend_y(ji,jk,jl) = diag%pot_vort_tend_y(ji,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(5,ji-1,jk,jl)*diag%w_placeholder(ji-1,jk,jl)*diag%eta_x(ji,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(6,ji-1,jk,jl)*diag%w_placeholder(ji-1,jk,jl+1)*diag%eta_x(ji,jk,jl+1) &
          + 0.5_wp*grid%inner_product_weights(5,ji,jk,jl)*diag%w_placeholder(ji,jk,jl)*diag%eta_x(ji,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(6,ji,jk,jl)*diag%w_placeholder(ji,jk,jl+1)*diag%eta_x(ji,jk,jl+1)
        enddo
        
        ! boundary conditions
        if (lperiodic) then
          diag%pot_vort_tend_y(1,jk,jl) = diag%pot_vort_tend_y(1,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(5,ny,jk,jl)*diag%w_placeholder(ny,jk,jl)*diag%eta_x(1,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(6,ny,jk,jl)*diag%w_placeholder(ny,jk,jl+1)*diag%eta_x(1,jk,jl+1) &
          + 0.5_wp*grid%inner_product_weights(5,1,jk,jl)*diag%w_placeholder(1,jk,jl)*diag%eta_x(1,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(6,1,jk,jl)*diag%w_placeholder(1,jk,jl+1)*diag%eta_x(1,jk,jl+1)
          diag%pot_vort_tend_y(ny+1,jk,jl) = diag%pot_vort_tend_y(1,jk,jl)
        endif
        
      enddo
    enddo
    !$omp end parallel do
    
    ! vertical velocity tendency due to horizontal vorticity and horizontal wind
    !$omp parallel do private(ji,jk,jl)
    do jl=2,n_layers
      do jk=1,nx
        do ji=1,ny
          diag%pot_vort_tend_z(ji,jk,jl) = 0.5_wp*( &
          grid%inner_product_weights(1,ji,jk,jl-1)*diag%u_placeholder(ji,jk+1,jl-1)*diag%eta_y(ji,jk+1,jl) &
          - grid%inner_product_weights(2,ji,jk,jl-1)*diag%v_placeholder(ji,jk,jl-1)*diag%eta_x(ji,jk,jl) &
          + grid%inner_product_weights(3,ji,jk,jl-1)*diag%u_placeholder(ji,jk,jl-1)*diag%eta_y(ji,jk,jl) &
          - grid%inner_product_weights(4,ji,jk,jl-1)*diag%v_placeholder(ji+1,jk,jl-1)*diag%eta_x(ji+1,jk,jl) &
          + grid%inner_product_weights(1,ji,jk,jl)*diag%u_placeholder(ji,jk+1,jl)*diag%eta_y(ji,jk+1,jl) &
          - grid%inner_product_weights(2,ji,jk,jl)*diag%v_placeholder(ji,jk,jl)*diag%eta_x(ji,jk,jl) &
          + grid%inner_product_weights(3,ji,jk,jl)*diag%u_placeholder(ji,jk,jl)*diag%eta_y(ji,jk,jl) &
          - grid%inner_product_weights(4,ji,jk,jl)*diag%v_placeholder(ji+1,jk,jl)*diag%eta_x(ji+1,jk,jl))
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine calc_vorticity_flux_term

end module mo_vorticity_flux








