! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module io_nml

  use definitions, only: wp
  
  implicit none
  
  integer           :: dt_write_min       ! output interval in minutes
  real(wp)          :: dt_write           ! output interval in seconds
  
  namelist /io/dt_write_min

  contains

  subroutine io_nml_setup
  
    ! local variables
    integer :: fileunit
    
    dt_write_min    = 60
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=io, unit=fileunit)
        
    close(fileunit)
    
    ! calculating the output timestep in seconds
    dt_write        = 60._wp*dt_write_min
  
  end subroutine io_nml_setup
  
end module io_nml









