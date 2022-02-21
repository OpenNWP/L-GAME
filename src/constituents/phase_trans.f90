! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module phase_trans

  ! In this module, phase transition rates are being calculated.
  
  use run_nml, only: nlins,ncols,nlays
  
  implicit none
  
  private
  
  public :: calc_h2otracers_source_rates
  
  contains
  
  subroutine calc_h2otracers_source_rates()
  
    ! This subroutine calculates the phase transition rates.
  
    ! local variables
    integer :: ji,jk,jl
  
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl)
    do ji=1,nlays
      do jk=1,ncols
        do jl=1,nlays
          
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine calc_h2otracers_source_rates

end module phase_trans











