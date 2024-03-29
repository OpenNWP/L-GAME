! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_run_nml
  
  ! This is the namelist that configures the basic run properties of a model integration.
  
  use mo_definitions, only: wp
  use mo_constants,   only: M_PI
  
  implicit none
  
  integer           :: ny                   ! number of gridpoints in y-direction
  integer           :: nx                   ! number of gridpoints in x-direction
  integer           :: n_layers             ! number of layers
  integer           :: n_levels             ! number of levels
  integer           :: n_oro_layers         ! number of layers following the orography
  integer           :: n_flat_layers        ! number of flat layers
  real(wp)          :: dy                   ! mesh size in y direction at sea level
  real(wp)          :: dx                   ! mesh size in x direction at sea level at the equator
  real(wp)          :: eff_hor_res          ! effective horizontal resolution
  real(wp)          :: dtime                ! time step
  real(wp)          :: toa                  ! top of atmosphere
  real(wp)          :: stretching_parameter ! vertical grid stretching parameter
  integer           :: run_span_min         ! run span in minutes
  real(wp)          :: t_init               ! epoch timestamp of the initialization
  integer           :: start_year           ! year of the model run beginning
  integer           :: start_month          ! month of the model run beginning
  integer           :: start_day            ! day of the model run beginning
  integer           :: start_hour           ! hour of the model run beginning
  integer           :: start_minute         ! minute of the model run beginning
  logical           :: lrestart             ! switch for restart runs
  logical           :: lcorio               ! switch for the Coriolis force
  logical           :: lideal               ! switch for analytic test cases
  logical           :: lplane               ! plane geometry switch
  character(len=64) :: scenario             ! scenario for ideal runs
  character(len=64) :: run_id               ! ID of the model run
  logical           :: llinear              ! switch to turn momentum advection on or off
  real(wp)          :: lat_center           ! latitude of the center of the model domain
  real(wp)          :: lon_center           ! longitude of the center of the model domain
  integer           :: theta_adv_order      ! theta advection order (2 or 3)
  logical           :: luse_bg_state        ! switch for using the hydrostatic background state
  
  namelist /run/ny,nx,n_layers,dy,dx,run_span_min,stretching_parameter,theta_adv_order, &
                toa,scenario,llinear,run_id,lcorio,n_oro_layers,lat_center,lon_center, &
                start_year,start_month,start_day,start_hour,start_minute,lplane,luse_bg_state
  
  contains
  
  subroutine run_nml_setup()
    
    ! local variables
    integer :: fileunit ! file unit of the namelist file
    
    ny = 35
    nx = 37
    n_layers = 50
    n_oro_layers = 40
    dy = 25e3_wp
    dx = 25e3_wp
    run_span_min = 1440
    start_year = 2000
    start_month = 1
    start_day = 1
    start_hour = 0
    start_minute = 0
    t_init = 0._wp
    toa = 40000._wp
    stretching_parameter = 1.3_wp
    lrestart = .false.
    lideal = .true.
    lplane = .false.
    scenario = "standard"
    run_id = "ideal"
    llinear = .false.
    lcorio = .true.
    lat_center = 0._wp
    lon_center = 0._wp
    theta_adv_order = 2
    luse_bg_state = .true.
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=run,unit=fileunit)
    
    close(fileunit)
    
    ! derived quantities
    n_levels = n_layers+1
    n_flat_layers = n_layers-n_oro_layers
    
    ! this calculates the time step using the CFL criterion
    eff_hor_res = sqrt(dx*dy)
    dtime = 1.61_wp*eff_hor_res/1000._wp
    
    ! calculating the Unix time of the model start
    t_init = (start_year-1970)*365*24*3600 + leap_year_correction(start_year)*24*3600 &
    + month_day_vector(start_month,start_year)*24*3600 + &
    (start_day-1)*24*3600 + start_hour*3600 + start_minute*60 &
    ! these are the leap seconds
    + 27
    
    ! sanity checks
    if (ny<3) then
      write(*,*) "Error: ny must be larger or equal than three. Aborting."
      call exit(1)
    endif
    if (nx<3) then
      write(*,*) "Error: nx must be larger or equal than three. Aborting."
      call exit(1)
    endif
    if (mod(ny, 2)==0) then
      write(*,*) "Error: ny must be odd. Aborting."
      call exit(1)
    endif
    if (mod(nx, 2)==0) then
      write(*,*) "Error: nx must be odd. Aborting."
      call exit(1)
    endif
    if (toa <= 0) then
      write(*,*) "Error: TOA must be positive."
      call exit(1)
    endif
    if (lat_center<-0.5_wp*M_PI .or. lat_center>0.5_wp*M_PI) then
      write(*,*) "Error: lat_center must be between -90 and 90."
      call exit(1)
    endif
    if (lon_center<-M_PI .or. lon_center>M_PI) then
      write(*,*) "Error: lon_center must be between -180 and 180."
      call exit(1)
    endif
    
    if (n_oro_layers>=n_layers) then
      write(*,*) "Error: it must be n_oro_layers<n_layers."
      call exit(1)
    endif
    
    if (theta_adv_order/=2 .and. theta_adv_order/=3) then
      write(*,*) "Error: it must be theta_adv_order=2 or theta_adv_order=3."
      call exit(1)
    endif
    
  end subroutine run_nml_setup
  
  function leap_year_correction(year)
    
    ! This is a helper function for calculating the Unix time.
    ! It returns the number of 29th of Februaries since 1970.
    
    integer, intent(in) :: year                 ! the year for which we want to calculate the leap year correction
    integer             :: leap_year_correction ! result
    
    leap_year_correction = (year - 1969)/4
    if (year>2000) then
      leap_year_correction = leap_year_correction - 1
    endif
    
  end function leap_year_correction
  
  function month_day_vector(month,year)
    
    ! This is a helper function for calculating the Unix time.
    ! It returns the amount of days in the wanted year in the previous months.
    
    integer, intent(in) :: month
    integer, intent(in) :: year
    integer             :: month_day_vector ! result
    
    ! local variables
    integer :: month_days(12)
    
    month_days(1) = 31
    month_days(2) = 28
    month_days(3) = 31
    month_days(4) = 30
    month_days(5) = 31
    month_days(6) = 30
    month_days(7) = 31
    month_days(8) = 31
    month_days(9) = 30
    month_days(10) = 31
    month_days(11) = 30
    month_days(12) = 31
    
    ! leap years
    if (year/4==0 .and. year/=2000 .and. month>2) then
      month_days(2) = 29
    endif
    
    month_day_vector = 0
    
    if (month>1) then
      month_day_vector = sum(month_days(1:(month-1)))
    endif
    
  end function month_day_vector
  
end module mo_run_nml












