! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module run_nml

  use definitions, only: wp
  
  implicit none

  integer           :: nlins               ! number of lines
  integer           :: ncols               ! number of columns
  integer           :: nlays               ! number of layers
  integer           :: nlays_oro           ! number of levels following the orography
  real(wp)          :: dy                  ! mesh size in y direction at sea level
  real(wp)          :: dx                  ! mesh size in x direction at sea level at the equator
  real(wp)          :: dtime               ! time step
  real(wp)          :: toa                 ! top of atmosphere
  real(wp)          :: sigma               ! vertical grid stretching parameter
  integer           :: run_span_hr         ! run span in hours
  real              :: t_init              ! epoch time stamp of the initialization
  logical           :: lrestart            ! switch for restart runs
  logical           :: lcorio              ! switch for the Coriolis force
  logical           :: lideal              ! switch for analytic test cases
  character(len=64) :: scenario            ! scenario for ideal runs
  character(len=64) :: run_id              ! ID of the model run
  logical           :: llinear             ! switch to turn momentum advection on or off
  real(wp)          :: impl_weight         ! implicit weight of the pressure gradient
  real(wp)          :: partial_impl_weight ! partial derivatives new time step weight
  real(wp)          :: PRANDTL_HEIGHT      ! height of the Prandtl layer
  real(wp)          :: lat_center          ! latitude of the center of the model domain
  real(wp)          :: lon_center          ! longitude of the center of the model domain
  real(wp)          :: x_dir_deg           ! direction of the x-axis of the model in degrees
  
  namelist /run/nlins,ncols,nlays,dy,dx,run_span_hr,sigma, &
  toa,scenario,llinear,run_id,lcorio,nlays_oro,lat_center,lon_center,x_dir_deg

  contains

  subroutine run_nml_setup
  
    ! local variables
    integer :: fileunit
    
    nlins               = 51
    ncols               = 51
    nlays               = 40
    nlays_oro           = 30
    dy                  = 48000._wp
    dx                  = 52000._wp
    run_span_hr         = 63
    t_init              = 0._wp
    toa                 = 40000._wp
    sigma               = 1.3_wp
    lrestart            = .false.
    lideal              = .true.
    scenario            = "standard"
    run_id              = "ideal"
    llinear             = .false.
    lcorio              = .true.
    impl_weight         = 0.75_wp
    partial_impl_weight = 0.5_wp
    PRANDTL_HEIGHT      = 100._wp
    lat_center          = 0._wp
    lon_center          = 0._wp
    x_dir_deg           = 90._wp
    
    ! open and read Namelist file
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=run, unit=fileunit)
        
    close(fileunit)
    
    ! this calculates the time step using the CFL criterion
    dtime           = 1.5_wp*dy/1000._wp
        
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









