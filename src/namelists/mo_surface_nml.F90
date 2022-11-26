! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_surface_nml
  
  ! This namelist defines the surface properties.
  
  use mo_definitions, only: wp
  use mo_diff_nml,    only: lmom_diff_h
  
  implicit none
  
  integer :: orography_id            ! identifies which orography to use
  logical :: lprog_soil_temp         ! switch for prognostic soil temperature
  integer :: nsoillays               ! number of soil layers
  logical :: lsfc_sensible_heat_flux ! surface sensible heat flux switch
  logical :: lsfc_phase_trans        ! surface phase transitions switch
  logical :: lpbl                    ! planetary boundary layer switch
  logical :: lsleve                  ! SLEVE vertical coordinate switch
  
  namelist /surface/orography_id,lprog_soil_temp,nsoillays,lsfc_sensible_heat_flux,lsfc_phase_trans,lpbl,lsleve
  
  contains
  
  subroutine surface_nml_setup
    
    ! local variables
    integer :: fileunit ! file unit of the namelist file
    
    ! default values
    orography_id = 1
    lprog_soil_temp = .true.
    nsoillays = 5
    lsfc_sensible_heat_flux = .true.
    lsfc_phase_trans = .true.
    lpbl = .true.
    lsleve = .false.
    
    ! Open and read namelist file.
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=surface,unit=fileunit)
    
    close(fileunit)
    
    ! sanity checks
    if (lpbl .and. .not. lmom_diff_h) then
      write(*,*) "Error: lmom_diff_h must be true is lpbl is true."
      call exit(1)
    endif
    
  end subroutine surface_nml_setup
  
end module mo_surface_nml





