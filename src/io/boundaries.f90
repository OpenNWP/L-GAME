! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module boundaries

  ! This module handles everything dealing with boundary conditions.

  use definitions, only: t_bc,wp
  use run_nml,     only: nlins,ncols
  use bc_nml,      only: n_swamp

  implicit none
  
  private
  
  public :: update_boundaries
  public :: rescale_tend
  public :: setup_bc_factor
  
  contains
  
  subroutine update_boundaries(bc)
  
    ! updates the boundary conditions
    
    ! input arguments and output
    type(t_bc) :: bc
    
  end subroutine update_boundaries
  
  subroutine rescale_tend
  
    ! rescales the tendencies to account for boundary conditions
  
  end subroutine rescale_tend
  
  subroutine setup_bc_factor(bc)
  
    ! This subroutine calculates the boundary conditions rescale factors.
  
    ! argument and output
    type(t_bc), intent(inout) :: bc
    
    ! local variables
    integer  :: ji,jk,jl           ! loop indices
    real(wp) :: dist_from_boundary ! index distance from the boundaary of the domain
    
    ! rescale factor for scalar fields
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols
        dist_from_boundary = min(ji-1,jk-1,nlins-ji,ncols-jk)
        bc%scalar_bc_factor(ji,jk) = max(1._wp - dist_from_boundary/n_swamp,0._wp)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! u rescale factor
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins
      do jk=1,ncols+1
        dist_from_boundary = min(ji-1,jk-1,nlins-ji,ncols+1-jk)
        bc%u_bc_factor(ji,jk) = max(1._wp - dist_from_boundary/n_swamp,0._wp)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! v rescale factor
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlins+1
      do jk=1,ncols
        dist_from_boundary = min(ji-1,jk-1,nlins+1-ji,ncols-jk)
        bc%v_bc_factor(ji,jk) = max(1._wp - dist_from_boundary/n_swamp,0._wp)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine setup_bc_factor

end module boundaries








