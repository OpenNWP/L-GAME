! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module surface_nml

  ! This namelist defines the surface properties.

  use definitions, only: wp
  use diff_nml,    only: lmom_diff_h
  
  implicit none
  
  integer :: orography_id            ! identifies which orography to use
  logical :: lsoil_heat_conduction   ! soil heat conduction switch
  integer :: nsoillays               ! number of soil layers
  logical :: lsfc_sensible_heat_flux ! surface sensible heat flux switch
  logical :: lsfc_phase_trans        ! surface phase transitions switch
  logical :: lpbl                    ! planetary boundary layer switch
  
  namelist /surface/orography_id,lsoil_heat_conduction,nsoillays,lsfc_phase_trans,lpbl
  
  contains
  
  subroutine surface_nml_setup

    ! local variables
    integer :: fileunit
    
    ! default values
    orography_id = 1
    lsoil_heat_conduction = .true.
    nsoillays = 5
    lsfc_sensible_heat_flux = .true.
    lsfc_phase_trans = .true.
    lpbl = .true.
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=surface, unit=fileunit)
        
    close(fileunit)
  
    ! sanity checks
    if (lpbl .and. .not. lmom_diff_h) then
      write(*,*) "Error: lmom_diff_h must be true is lpbl is true."
      call exit(1)
    endif
  
  end subroutine surface_nml_setup

end module surface_nml
