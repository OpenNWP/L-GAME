! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_divergence_operators
  
  ! The calculation of the horizontal divergence operator is executed in this module.
  
  use mo_definitions, only: wp,t_grid
  use mo_run_nml,     only: ny,nx,n_layers,n_oro_layers,dtime,n_flat_layers
  use mo_averaging,   only: vertical_contravariant_corr
  
  implicit none
  
  contains

  subroutine div_h(vector_field_x,vector_field_y,result_field,grid)
    
    ! This subroutine computes the divergence of a vector field.
    
    real(wp),     intent(in)  :: vector_field_x(:,:,:) ! x-component of horizontal vector field of which to calculate the divergence
    real(wp),     intent(in)  :: vector_field_y(:,:,:) ! y-component of horizontal vector field of which to calculate the divergence
    real(wp),     intent(out) :: result_field(:,:,:)   ! resulting scalar field
    type(t_grid), intent(in)  :: grid                  ! the grid properties
    
    ! local variables
    integer  :: ji           ! horizontal index
    integer  :: jk           ! horizontal index
    integer  :: jl           ! layer index
    real(wp) :: comp_h       ! horizontal component of divergence
    real(wp) :: comp_v       ! vertical component of divergence
    real(wp) :: contra_upper ! contravariant mass flux density resulting from the horizontal vector components through the upper area
    real(wp) :: contra_lower ! contravariant mass flux density resulting from the horizontal vector components through the lower area

    !$omp parallel do private(ji,jk,jl,contra_upper,contra_lower,comp_h,comp_v)
    do jl=1,n_layers
      do jk=1,nx
        do ji=1,ny
          ! the horizontal component
          comp_h = &
          vector_field_x(ji,jk+1,jl)*grid%area_x(ji,jk+1,jl) &
          + vector_field_y(ji,jk,jl)*grid%area_y(ji,jk,jl) &
          - vector_field_x(ji,jk,jl)*grid%area_x(ji,jk,jl) &
          - vector_field_y(ji+1,jk,jl)*grid%area_y(ji+1,jk,jl)
          
          ! the vertical component
          comp_v = 0._wp
          if (jl==n_flat_layers) then
            contra_lower = vertical_contravariant_corr(vector_field_x,vector_field_y,ji,jk,jl+1,grid)
            comp_v = -contra_lower*grid%area_z(ji,jk,jl+1)
          elseif (jl==n_layers) then
            contra_upper = vertical_contravariant_corr(vector_field_x,vector_field_y,ji,jk,jl,grid)
            comp_v = contra_upper*grid%area_z(ji,jk,jl)
          elseif (jl>n_flat_layers) then
            contra_upper = vertical_contravariant_corr(vector_field_x,vector_field_y,ji,jk,jl,grid)
            contra_lower = vertical_contravariant_corr(vector_field_x,vector_field_y,ji,jk,jl+1,grid)
            comp_v &
            = contra_upper*grid%area_z(ji,jk,jl) &
            - contra_lower*grid%area_z(ji,jk,jl+1)
          endif
          
          ! adding the horizontal and the vertical component and dividing by the volume
          result_field(ji,jk,jl) = (comp_h + comp_v)/grid%volume(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine div_h

  subroutine div_h_tracers(vector_field_x,vector_field_y,density_field,wind_field_x,wind_field_y,result_field,grid)

    ! This subroutine computes the divergence of a vector field for tracers.
    
    real(wp),     intent(in)  :: vector_field_x(:,:,:) ! x-component of horizontal vector field of which to calculate the divergence
    real(wp),     intent(in)  :: vector_field_y(:,:,:) ! y-component of horizontal vector field of which to calculate the divergence
    real(wp),     intent(out) :: result_field(:,:,:)   ! resulting scalar field
    real(wp),     intent(in)  :: density_field(:,:,:)  ! density field at the current time step
    real(wp),     intent(in)  :: wind_field_x(:,:,:)   ! x-component of horizontal wind field
    real(wp),     intent(in)  :: wind_field_y(:,:,:)   ! y-component of horizontal wind field
    type(t_grid), intent(in)  :: grid                  ! the grid properties
    
    ! local variables
    integer  :: ji            ! horizontal index
    integer  :: jk            ! horizontal index
    integer  :: jl            ! layer index
    real(wp) :: comp_h        ! horizontal component of divergence
    real(wp) :: comp_v        ! vertical component of divergence
    real(wp) :: contra_upper  ! contravariant mass flux density resulting from the horizontal vector components through the upper area
    real(wp) :: contra_lower  ! contravariant mass flux density resulting from the horizontal vector components through the lower area
    real(wp) :: density_upper ! density at the upper interface
    real(wp) :: density_lower ! density at the lower interface

    !$omp parallel do private(ji,jk,jl,contra_upper,contra_lower,comp_h,comp_v,density_upper,density_lower)
    do jl=1,n_layers
      do jk=1,nx
        do ji=1,ny
          ! the horizontal component
          comp_h = &
          vector_field_x(ji,jk+1,jl)*grid%area_x(ji,jk+1,jl) &
          + vector_field_y(ji,jk,jl)*grid%area_y(ji,jk,jl) &
          - vector_field_x(ji,jk,jl)*grid%area_x(ji,jk,jl) &
          - vector_field_y(ji+1,jk,jl)*grid%area_y(ji+1,jk,jl)
          
          ! the vertical component
          comp_v = 0._wp
          if (jl==n_flat_layers) then
            contra_lower = vertical_contravariant_corr(wind_field_x,wind_field_y,ji,jk,jl+1,grid)
            if (contra_lower<=0._wp) then
              density_lower = density_field(ji,jk,jl)
            else
              density_lower = density_field(ji,jk,jl+1)
            endif
            comp_v = -density_lower*contra_lower*grid%area_z(ji,jk,jl+1)
          elseif (jl==n_layers) then
            contra_upper = vertical_contravariant_corr(wind_field_x,wind_field_y,ji,jk,jl,grid)
            if (contra_upper<=0._wp) then
              density_upper = density_field(ji,jk,jl-1)
            else
              density_upper = density_field(ji,jk,jl)
            endif
            comp_v = density_upper*contra_upper*grid%area_z(ji,jk,jl)
          elseif (jl>n_flat_layers) then
            contra_upper = vertical_contravariant_corr(wind_field_x,wind_field_y,ji,jk,jl,grid)
            if (contra_upper<=0._wp) then
              density_upper = density_field(ji,jk,jl-1)
            else
              density_upper = density_field(ji,jk,jl)
            endif
            contra_lower = vertical_contravariant_corr(wind_field_x,wind_field_y,ji,jk,jl+1,grid)
            if (contra_lower<=0._wp) then
              density_lower = density_field(ji,jk,jl)
            else
              density_lower = density_field(ji,jk,jl+1)
            endif
            comp_v &
            = density_upper*contra_upper*grid%area_z(ji,jk,jl) &
            - density_lower*contra_lower*grid%area_z(ji,jk,jl+1)
          endif
          
          ! adding the horizontal and the vertical component and dividing by the volume
          result_field(ji,jk,jl) = (comp_h + comp_v)/grid%volume(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do

  end subroutine div_h_tracers
  
  subroutine add_vertical_div(in_field,out_field,grid)

    ! This subroutine adds the divergence of the vertical component of a vector field to the input scalar field.

    real(wp),     intent(in)    :: in_field(:,:,:)  ! input vertical vector field
    real(wp),     intent(inout) :: out_field(:,:,:) ! scalar field to which the vertical divergence will be added
    type(t_grid), intent(in)    :: grid             ! grid properties

    ! local variables
    integer  :: ji        ! horizontal index
    integer  :: jk        ! horizontal index
    integer  :: jl        ! layer index
    real(wp) :: cov_upper ! covariant upper component of the vector field
    real(wp) :: cov_lower ! covariant lower component of the vector field
    real(wp) :: comp_v    ! vertical divergence value to add
    
    !$omp parallel do private(ji,jk,jl,cov_upper,cov_lower,comp_v)
    do jl=1,n_layers
      do jk=1,nx
        do ji=1,ny
          if (jl==1) then
            cov_upper = 0._wp
            cov_lower = in_field(ji,jk,jl+1)
          elseif (jl==n_layers) then
            cov_upper = in_field(ji,jk,jl)
            cov_lower = 0._wp
          else
            cov_upper = in_field(ji,jk,jl)
            cov_lower = in_field(ji,jk,jl+1)
          endif
          comp_v = cov_upper*grid%area_z(ji,jk,jl) - cov_lower*grid%area_z(ji,jk,jl+1)
          out_field(ji,jk,jl) = out_field(ji,jk,jl) + comp_v/grid%volume(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine add_vertical_div

end module mo_divergence_operators










