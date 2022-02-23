! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module divergence_operators

  use definitions, only: wp,t_grid
  use run_nml,     only: nlins,ncols,nlays,nlays_oro
  use averaging,   only: vertical_contravariant_corr
  
  implicit none
  
  private
  
  public :: divv_h
  
  contains

  subroutine divv_h(vector_field_x,vector_field_y,result_field,grid)

    ! This subroutine computes the gradient of a scalar field.
    real(wp),     intent(in)    :: vector_field_x(:,:,:) ! x-component of horizontal vector field of which to calculate the divergence
    real(wp),     intent(in)    :: vector_field_y(:,:,:) ! y-component of horizontal vector field of which to calculate the divergence
    real(wp),     intent(inout) :: result_field(:,:,:)   ! resulting scalar field
    type(t_grid), intent(in)    :: grid                  ! the grid properties
    ! local variables
    integer                     :: ji,jk,jl              ! loop variables
    real(wp)                    :: comp_h                ! horizontal component of divergence
    real(wp)                    :: comp_v                ! horizontal component of divergence
    real(wp)                    :: contra_upper          ! contravariant mass flux density resulting 
                                                         ! from the horizontal vector components through the upper area
    real(wp)                    :: contra_lower          ! contravariant mass flux density resulting
                                                         ! from the horizontal vector components through the lower area

    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,comp_h,comp_v)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
          ! the horizontal component
          comp_h = &
          vector_field_x(ji,jk+1,jl)*grid%area_x(ji,jk+1,jl) &
          + vector_field_y(ji,jk+1,jl)*grid%area_y(ji,jk+1,jl) &
          - vector_field_x(ji,jk,jl)*grid%area_x(ji,jk,jl) &
          - vector_field_y(ji,jk,jl)*grid%area_y(ji,jk,jl)
          
          ! the vertical component
          comp_v = 0._wp
          if (jl == nlays-nlays_oro) then
            contra_lower = vertical_contravariant_corr(vector_field_x,vector_field_y,ji,jk,jl+1,grid)
            comp_v = -contra_lower*grid%area_z(ji,jk,jl+1)
          elseif (jl == nlays) then
            contra_upper = vertical_contravariant_corr(vector_field_x,vector_field_y,ji,jk,jl,grid)
            comp_v = contra_upper*grid%area_z(ji,jk,jl)
          elseif (jl > nlays-nlays_oro-1) then
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
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine divv_h

end module divergence_operators










