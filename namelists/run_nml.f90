! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module run_nml

  use definitions, only: wp
  
  implicit none

  integer           :: nlins               ! number of lines
  integer           :: ncols               ! number of columns
  integer           :: nlays               ! number of levels
  integer           :: nlays_oro           ! number of levels following the orography
  real(wp)          :: dy                  ! mesh size in y direction at sea level
  real(wp)          :: dx                  ! mesh size in x direction at sea level at the equator
  real(wp)          :: dtime               ! time step
  real(wp)          :: toa                 ! top of atmosphere
  real(wp)          :: sigma               ! vertical grid stretching parameter
  integer           :: run_span_hr         ! run span in hours
  real              :: t_init              ! epoch time stamp of the initialization
  integer           :: adv_sound_ratio     ! ratio of advective to sound time step
  real(wp)          :: semimajor           ! large halfaxis of the Earth
  real(wp)          :: semiminor           ! small halfaxis of the Earth
  real(wp)          :: re                  ! Earth radius
  logical           :: lrestart            ! switch for restart runs
  logical           :: lcorio              ! switch for the Coriolis force
  logical           :: lideal              ! switch for analytic test cases
  logical           :: l3dvar              ! switch for 3d-Var
  logical           :: l4dvar              ! switch for 4d-Var
  character(len=64) :: scenario            ! scenario for ideal runs
  character(len=64) :: run_id              ! ID of the model run
  real(wp)          :: p_0                 ! reference pressure
  real(wp)          :: omega               ! angular frequency of Earth rotation
  real(wp)          :: lapse_rate          ! lapse_rate within the troposphere
  real(wp)          :: surface_temp        ! the temperature at the surface
  real(wp)          :: tropo_height        ! the tropopause height
  real(wp)          :: inv_height          ! height where the temperature inversion begins
  real(wp)          :: t_grad_inv          ! temperature gradient above the inversion
  real(wp)          :: p_0_standard        ! reference pressure of the standard atmosphere
  real(wp)          :: gravity             ! average surface gravity value
  logical           :: llinear             ! switch to turn momentum advection on or off
  real(wp)          :: impl_weight         ! implicit weight of the pressure gradient
  real(wp)          :: partial_impl_weight ! partial derivatives new time step weight
  
  namelist /run/nlins,ncols,nlays,dy,dx,run_span_hr, &
  adv_sound_ratio,toa,scenario,llinear,run_id,lcorio

  contains

  subroutine run_nml_setup
  
    ! local variables
    integer :: fileunit
    
    nlins               = 101
    ncols               = 101
    nlays               = 80
    dy                  = 800._wp
    dx                  = 850._wplcorio
    run_span_hr         = 63
    t_init              = 0._wp
    adv_sound_ratio     = 1
    toa                 = 40000._wp
    sigma               = 1.3_wp
    semiminor           = 6356752.314_wp
    semimajor           = 6378137.0_wp
    re                  = (semimajor*semimajor*semiminor)**(1._wp/3._wp)
    lrestart            = .false.
    lideal              = .true.
    l3dvar              = .false.
    l4dvar              = .false.
    scenario            = "standard"
    run_id              = "ideal"
    p_0                 = 100000._wp
    omega               = 7.292115e-5
    lapse_rate          = 0.0065_wp
    surface_temp        = 288.15_wp
    tropo_height        = 12000._wp
    inv_height          = 20000._wp
    t_grad_inv          = 0.001_wp
    p_0_standard        = 101325._wp
    gravity             = 9.80616_wp
    llinear             = .false.
    lcorio              = .true.
    impl_weight         = 0.75_wp
    partial_impl_weight = 0.5_wp
    
    ! Open and read Namelist file.
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=run, unit=fileunit)
        
    close(fileunit)
    
    ! this calculates the time step using the CFL criterion
    dtime           = 0.7_wp*dy/350._wp
    ! number of layers following orography
    nlays_oro       = int(0.66*nlays)
        
    ! checking input data for correctness
    if (mod(nlins, 2) == 0) then
      write(*,*) "Error: nlins must be odd. Aborting."
      call exit(1)
    endif
    if (mod(ncols, 2) == 0) then
      write(*,*) "Error: ncols must be odd. Aborting."
      call exit(1)
    endif
    if (toa <= 0) then
      write(*,*) "Error: TOA must be positive."
      call exit(1)
    endif
    
    write(*,*) "Time step: ", dtime, " s."
    
  end subroutine run_nml_setup
  
end module run_nml









