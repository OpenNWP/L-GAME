! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module hetero_nml

  use definitions, only: wp

  implicit none
  
  integer :: no_of_gaseous_constituents ! number of constituents of the gas phase
  
  namelist /hetero/no_of_gaseous_constituents

  contains

  subroutine hetero_nml_setup
  
    no_of_gaseous_constituents = 1
    
  end subroutine hetero_nml_setup
  
end module hetero_nml






