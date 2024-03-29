! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_constituents_nml
  
  ! This namelist defines the constituents of the model atmosphere.
  
  use mo_definitions, only: wp
  
  implicit none
  
  logical  :: lmoist                   ! moisture switch
  integer  :: n_gaseous_constituents   ! number of constituents of the gas phase
  integer  :: n_condensed_constituents ! number of condensed constituents
  integer  :: n_constituents           ! the total number of constituents
  
  namelist /constituents/lmoist
  
  contains
  
  subroutine constituents_nml_setup()
    
    ! local variables
    integer :: fileunit ! file unit of the namelist file
    
    lmoist = .true.
    n_condensed_constituents = 5
    n_gaseous_constituents = 2
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=constituents,unit=fileunit)
    
    close(fileunit)
    
    ! the dry case
    if (.not. lmoist) then
      n_condensed_constituents = 0
      n_gaseous_constituents = 1
    endif
    
    ! calculating the total number of constituents
    n_constituents = n_condensed_constituents + n_gaseous_constituents
    
  end subroutine constituents_nml_setup
  
end module mo_constituents_nml






