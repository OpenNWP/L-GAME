! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module vorticity_flux

  ! This module computes the voriticity flux term.
  
  use definitions, only: t_grid,t_diag,wp
  use run_nml,     only: nlins,ncols,nlays
  
  implicit none
  
  private
  
  public :: calc_vorticity_flux_term
  
  contains
  
  subroutine calc_vorticity_flux_term(diag,grid)
  
    ! This module computes the voriticity flux.

    ! input arguments and output
    type(t_diag), intent(inout) :: diag ! diagnostic quantities
    type(t_grid), intent(in)    :: grid ! model grid
    
    ! local variables
    integer :: ji,jk,jl ! loop indices
    
    ! horizontal velocity tendency due to vertical vorticity and horizontal wind (TRSK)
    ! u
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=1,nlins
      do jk=2,ncols
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
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! v
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk)
    do ji=2,nlins
      do jk=1,ncols
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
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! horizontal velocity tendency due to horizontal vorticity and vertical wind
    ! u
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols-1
        do jl=1,nlays
          diag%pot_vort_tend_x(ji,jk,jl) = diag%pot_vort_tend_x(ji,jk,jl) &
          - 0.5_wp*grid%inner_product_weights(ji,jk,jl,5)*diag%w_placeholder(ji,jk,jl)*diag%eta_y(ji,jk+1,jl) &
          - 0.5_wp*grid%inner_product_weights(ji,jk,jl,6)*diag%w_placeholder(ji,jk,jl+1)*diag%eta_y(ji,jk+1,jl+1) &
          - 0.5_wp*grid%inner_product_weights(ji,jk+1,jl,5)*diag%w_placeholder(ji,jk+1,jl)*diag%eta_y(ji,jk+1,jl) &
          - 0.5_wp*grid%inner_product_weights(ji,jk+1,jl,6)*diag%w_placeholder(ji,jk+1,jl+1)*diag%eta_y(ji,jk+1,jl+1)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    ! v
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins-1
      do jk=1,ncols
        do jl=1,nlays
          diag%pot_vort_tend_y(ji,jk,jl) = diag%pot_vort_tend_y(ji,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(ji,jk,jl,5)*diag%w_placeholder(ji,jk,jl)*diag%eta_x(ji+1,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(ji,jk,jl,6)*diag%w_placeholder(ji,jk,jl+1)*diag%eta_x(ji+1,jk,jl+1) &
          + 0.5_wp*grid%inner_product_weights(ji+1,jk,jl,5)*diag%w_placeholder(ji+1,jk,jl)*diag%eta_x(ji+1,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(ji+1,jk,jl,6)*diag%w_placeholder(ji+1,jk,jl+1)*diag%eta_x(ji+1,jk,jl+1)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! vertical velocity tendency due to horizontal vorticity and horizontal wind
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        do jl=2,nlays
          diag%pot_vort_tend_z(ji,jk,jl) = 0.5_wp*( &
          + grid%inner_product_weights(ji,jk,jl-1,1)*diag%u_placeholder(ji,jk+1,jl-1)*diag%eta_y(ji,jk+1,jl) &
          - grid%inner_product_weights(ji,jk,jl-1,2)*diag%v_placeholder(ji,jk,jl-1)*diag%eta_x(ji,jk,jl) &
          + grid%inner_product_weights(ji,jk,jl-1,3)*diag%u_placeholder(ji,jk,jl-1)*diag%eta_y(ji,jk,jl) &
          - grid%inner_product_weights(ji,jk,jl-1,4)*diag%v_placeholder(ji+1,jk,jl-1)*diag%eta_x(ji+1,jk,jl) &
          + grid%inner_product_weights(ji,jk,jl,1)*diag%u_placeholder(ji,jk+1,jl)*diag%eta_y(ji,jk+1,jl) &
          - grid%inner_product_weights(ji,jk,jl,2)*diag%v_placeholder(ji,jk,jl)*diag%eta_x(ji+1,jk,jl) &
          + grid%inner_product_weights(ji,jk,jl,3)*diag%u_placeholder(ji,jk,jl)*diag%eta_y(ji,jk,jl) &
          - grid%inner_product_weights(ji,jk,jl,4)*diag%v_placeholder(ji+1,jk,jl)*diag%eta_x(ji+1,jk,jl))
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine calc_vorticity_flux_term

end module vorticity_flux








