! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module constituents_nml

  use definitions, only: wp

  implicit none
  
  integer :: no_of_gaseous_constituents   ! number of constituents of the gas phase
  integer :: no_of_condensed_constituents ! number of condensed constituents
  integer :: no_of_constituents           ! the total number of constituents
  logical :: lassume_lte                  ! switch for the local thermodynamic equilibrium option
  
  public :: no_of_constituents
    
  namelist /constituents/no_of_condensed_constituents,no_of_gaseous_constituents,lassume_lte

  contains

  subroutine constituents_nml_setup()
  
    ! local variables
    integer        :: fileunit
    
    ! default values
    no_of_condensed_constituents = 4
    no_of_gaseous_constituents = 2
    lassume_lte = .true.
    
    ! open and read namelist file
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=constituents, unit=fileunit)
    
    close(fileunit)
    
    no_of_constituents = no_of_condensed_constituents + no_of_gaseous_constituents
    
  end subroutine constituents_nml_setup
  
end module constituents_nml






