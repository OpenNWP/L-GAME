! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module surface_nml

  ! This namelist defines the surface properties.

  use definitions, only: wp
  
  implicit none
  
  integer :: orography_id ! identifies which orography to use
  logical :: lsoil        ! soil switch
  integer :: nsoillays    ! number of soil layers
  
  namelist /surface/orography_id,lsoil,nsoillays
  
  contains
  
  subroutine surface_nml_setup

    ! local variables
    integer :: fileunit
    
    ! default values
    orography_id = 1
    lsoil = .true.
    nsoillays = 5
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=surface, unit=fileunit)
        
    close(fileunit)
  
  end subroutine surface_nml_setup

end module surface_nml
