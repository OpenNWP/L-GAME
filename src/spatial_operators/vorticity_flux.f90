! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module vorticity_flux

  ! This module computes the voriticity flux term.
  
  use definitions, only: t_grid,t_diag,wp
  use run_nml,     only: nlins,ncols,nlays
  
  implicit none
  
  private
  
  public :: calc_vorticity_flux_term
  
  contains
  
  subroutine calc_vorticity_flux_term(diag,grid)

    type(t_diag), intent(inout) :: diag ! diagnostic quantities
    type(t_grid), intent(in)    :: grid ! model grid
    
    ! local variables
    integer                     :: ji,jk,jl ! loop indices
    
    ! horizontal velocity tendency due to vertical vorticity and horizontal wind (TRSK)
    ! u
    do ji=1,nlins
      do jk=1,ncols-1
        diag%pot_vort_tend_x(ji,jk,:) = &
        grid%trsk_weights_u(ji,jk,1)*diag%v_placeholder(ji+1,jk+1,:)*0.25_wp* &
        (diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk,:)) &
        + grid%trsk_weights_u(ji,jk,2)*diag%u_placeholder(ji+2,jk,:)*0.25_wp* &
        (diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji,jk,:)) &
        + grid%trsk_weights_u(ji,jk,3)*diag%v_placeholder(ji,jk+1,:)*0.25_wp* &
        (diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji,jk,:)+diag%z_eta_z(ji,jk+1,:)) &
        + grid%trsk_weights_u(ji,jk,4)*diag%v_placeholder(ji,jk+2,:)*0.25_wp* &
        (diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji,jk+2,:)) &
        + grid%trsk_weights_u(ji,jk,5)*diag%u_placeholder(ji+1,jk+2,:)*0.25_wp* &
        (diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji,jk+2,:)+diag%z_eta_z(ji+1,jk+2,:)) &
        + grid%trsk_weights_u(ji,jk,6)*diag%v_placeholder(ji+1,jk+2,:)*0.25_wp* &
        (diag%z_eta_z(ji,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk+2,:)+diag%z_eta_z(ji+1,jk+1,:))
      enddo
    enddo
    ! v
    do ji=1,nlins-1
      do jk=1,ncols
        diag%pot_vort_tend_y(ji,jk,:) = &
        grid%trsk_weights_v(ji,jk,1)*diag%u_placeholder(ji+1,jk,:)*0.25_wp* &
        (diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji,jk,:)) &
        + grid%trsk_weights_v(ji,jk,2)*diag%u_placeholder(ji+1,jk+1,:)*0.25_wp* &
        (diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji,jk+1,:)) &
        + grid%trsk_weights_v(ji,jk,3)*diag%u_placeholder(ji+2,jk+1,:)*0.25_wp* &
        (diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+2,jk+1,:)) &
        + grid%trsk_weights_v(ji,jk,4)*diag%u_placeholder(ji+2,jk,:)*0.25_wp* &
        (diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+1,jk+1,:)+diag%z_eta_z(ji+1,jk,:)+diag%z_eta_z(ji+2,jk,:))
      enddo
    enddo
    
    ! horizontal velocity tendency due to horizontal vorticity and vertical wind
    ! u
    do ji=1,nlins
      do jk=1,ncols-1
        do jl=1,nlays
          diag%pot_vort_tend_x(ji,jk,jl) = diag%pot_vort_tend_x(ji,jk,jl) &
          - 0.5_wp*grid%inner_product_weights(ji,jk  ,jl,5)*diag%w_placeholder(ji,jk  ,jl  )*diag%z_eta_y(ji,jk+1,jl  ) &
          - 0.5_wp*grid%inner_product_weights(ji,jk  ,jl,6)*diag%w_placeholder(ji,jk  ,jl+1)*diag%z_eta_y(ji,jk+1,jl+1) &
          - 0.5_wp*grid%inner_product_weights(ji,jk+1,jl,5)*diag%w_placeholder(ji,jk+1,jl  )*diag%z_eta_y(ji,jk+1,jl  ) &
          - 0.5_wp*grid%inner_product_weights(ji,jk+1,jl,6)*diag%w_placeholder(ji,jk+1,jl+1)*diag%z_eta_y(ji,jk+1,jl+1)
        enddo
      enddo
    enddo
    ! v
    do ji=1,nlins-1
      do jk=1,ncols
        do jl=1,nlays
          diag%pot_vort_tend_y(ji,jk,jl) = diag%pot_vort_tend_y(ji,jk,jl) &
          + 0.5_wp*grid%inner_product_weights(ji  ,jk,jl,5)*diag%w_placeholder(ji  ,jk,jl  )*diag%z_eta_x(ji+1,jk,jl  ) &
          + 0.5_wp*grid%inner_product_weights(ji  ,jk,jl,6)*diag%w_placeholder(ji  ,jk,jl+1)*diag%z_eta_x(ji+1,jk,jl+1) &
          + 0.5_wp*grid%inner_product_weights(ji+1,jk,jl,5)*diag%w_placeholder(ji+1,jk,jl  )*diag%z_eta_x(ji+1,jk,jl  ) &
          + 0.5_wp*grid%inner_product_weights(ji+1,jk,jl,6)*diag%w_placeholder(ji+1,jk,jl+1)*diag%z_eta_x(ji+1,jk,jl+1)
        enddo
      enddo
    enddo
    
    ! vertical velocity tendency due to horizontal vorticity and horizontal wind
    do ji = 1,nlins
      do jk=1,ncols
        do jl=2,nlays-1
          diag%pot_vort_tend_z(ji,jk,jl) = 0.5_wp*( &
          + grid%inner_product_weights(ji,jk,jl-1,1)*diag%u_placeholder(ji+1,jk+1,jl-1)*diag%z_eta_y(ji  ,jk+1,jl) &
          - grid%inner_product_weights(ji,jk,jl-1,2)*diag%v_placeholder(ji+1,jk+1,jl-1)*diag%z_eta_x(ji+1,jk  ,jl) &
          + grid%inner_product_weights(ji,jk,jl-1,3)*diag%u_placeholder(ji+1,jk  ,jl-1)*diag%z_eta_y(ji  ,jk  ,jl) &
          - grid%inner_product_weights(ji,jk,jl-1,4)*diag%v_placeholder(ji  ,jk+1,jl-1)*diag%z_eta_x(ji  ,jk  ,jl) &
          + grid%inner_product_weights(ji,jk,jl  ,1)*diag%u_placeholder(ji+1,jk+1,jl  )*diag%z_eta_y(ji  ,jk+1,jl) &
          - grid%inner_product_weights(ji,jk,jl  ,2)*diag%v_placeholder(ji+1,jk+1,jl  )*diag%z_eta_x(ji+1,jk  ,jl) &
          + grid%inner_product_weights(ji,jk,jl  ,3)*diag%u_placeholder(ji+1,jk  ,jl  )*diag%z_eta_y(ji  ,jk  ,jl) &
          - grid%inner_product_weights(ji,jk,jl  ,4)*diag%v_placeholder(ji  ,jk+1,jl  )*diag%z_eta_x(ji  ,jk  ,jl))
        enddo
      enddo
    enddo
  
  end subroutine calc_vorticity_flux_term

end module vorticity_flux








