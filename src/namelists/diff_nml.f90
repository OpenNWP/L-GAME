! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module diff_nml
  
  ! In this namelist, the diffusion properties are configured.
  
  use definitions, only: wp
  
  implicit none
  
  logical  :: lklemp           ! turns the Klemp damping layer on or off
  real(wp) :: klemp_damp_max   ! the maximum Klemp damping coefficient
  real(wp) :: klemp_begin_rel  ! lower boundary of the Klemp damping layer in relation to the TOA
  logical  :: lmom_diff_h      ! switch for horizontal momentum diffusion
  logical  :: lmom_diff_v      ! switch for vertical momentum diffusion
  real(wp) :: diff_h_smag_div  ! horizontal diffusion Smagorinsky factor acting on divergent movements
  real(wp) :: diff_h_smag_rot  ! horizontal diffusion Smagorinsky factor acting on vortical movements
  
  namelist /diff/lklemp,klemp_damp_max,klemp_begin_rel,lmom_diff_h,lmom_diff_v,diff_h_smag_div,diff_h_smag_rot
  
  contains
  
  subroutine diff_nml_setup

    ! local variables
    integer :: fileunit
    
    ! default values
    lklemp          = .true.
    klemp_damp_max  = 0.25_wp
    klemp_begin_rel = 0.53_wp
    lmom_diff_h     = .true.
    lmom_diff_v     = .true.
    diff_h_smag_div = 0.01_wp
    diff_h_smag_rot = 0.01_wp
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=diff, unit=fileunit)
        
    close(fileunit)
  
  end subroutine diff_nml_setup
  
end module diff_nml
