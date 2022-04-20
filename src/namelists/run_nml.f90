! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module run_nml

  ! This is the namelists the configures the basic run properties of a model integration.
  
  use iso_c_binding
  use definitions, only: wp
  use constants,   only: M_PI
  
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
  integer           :: start_year          ! year when to begin the model run
  integer           :: start_month         ! month when to begin the model run
  integer           :: start_day           ! day when to begin the model run
  integer           :: start_hour          ! hour when to begin the model run
  integer           :: start_minute        ! minute when to begin the model run
  logical           :: lrestart            ! switch for restart runs
  logical           :: lcorio              ! switch for the Coriolis force
  logical           :: lideal              ! switch for analytic test cases
  logical           :: lplane              ! plane geometry switch
  character(len=64) :: scenario            ! scenario for ideal runs
  character(len=64) :: run_id              ! ID of the model run
  logical           :: llinear             ! switch to turn momentum advection on or off
  real(wp)          :: impl_weight         ! implicit weight of the pressure gradient
  real(wp)          :: partial_impl_weight ! partial derivatives new time step weight
  real(wp)          :: lat_center          ! latitude of the center of the model domain
  real(wp)          :: lon_center          ! longitude of the center of the model domain
  
  namelist /run/nlins,ncols,nlays,dy,dx,run_span_hr,sigma, &
  toa,scenario,llinear,run_id,lcorio,nlays_oro,lat_center,lon_center, &
  start_year,start_month,start_day,start_hour,start_minute,lplane

  contains

  subroutine run_nml_setup
  
    ! local variables
    integer :: fileunit
    
    nlins = 25
    ncols = 25
    nlays = 50
    nlays_oro = 40
    dy = 25e3_wp
    dx = 25e3_wp
    run_span_hr = 84
    start_year = 2000
    start_month = 1
    start_day = 1
    start_hour = 0
    start_minute = 0
    t_init = 0._wp
    toa = 40000._wp
    sigma = 1.3_wp
    lrestart = .false.
    lideal = .true.
    lplane = .false.
    scenario = "standard"
    run_id = "ideal"
    llinear = .false.
    lcorio = .true.
    impl_weight = 0.75_wp
    partial_impl_weight = 0.5_wp
    lat_center = 0._wp
    lon_center = 0._wp
    
    ! open and read Namelist file
    open(action="read", file="namelist.nml", newunit=fileunit)
    read(nml=run, unit=fileunit)
        
    close(fileunit)
    
    ! this calculates the time step using the CFL criterion
    dtime           = 1.61_wp*sqrt(dx*dy)/1000._wp
    
    ! calculating the Unix time of the model start
    t_init = (start_year-1970)*365*24*3600 + leap_year_correction(start_year)*24*3600 &
    + month_day_vector(start_month,start_year)*24*3600 + &
    (start_day-1)*24*3600 + start_hour*3600 + start_minute*60 &
    ! these are the leap seconds
    + 27
    
    ! sanity checks
    if (nlins<3) then
      write(*,*) "Error: nlins must be larger or equal than three. Aborting."
      call exit(1)
    endif
    if (ncols<3) then
      write(*,*) "Error: ncols must be larger or equal than three. Aborting."
      call exit(1)
    endif
    if (mod(nlins, 2)==0) then
      write(*,*) "Error: nlins must be odd. Aborting."
      call exit(1)
    endif
    if (mod(ncols, 2)==0) then
      write(*,*) "Error: ncols must be odd. Aborting."
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
    
    write(*,*) "Time step: ", dtime, " s."
    
  end subroutine run_nml_setup
  
  function leap_year_correction(year)
  
    ! This is a helper function for calculating the Unix time.
    ! It returns the number of 29th of Februaries since 1970.
    
    ! input
    integer, intent(in) :: year ! the year for which we want to calculate the leap year correction
    ! output
    integer             :: leap_year_correction
    
    leap_year_correction = (year - 1969)/4
    if (year>2000) then
      leap_year_correction = leap_year_correction - 1
    endif
  
  end function leap_year_correction
  
  function month_day_vector(month,year)
  
    ! This is a helper function for calculating the Unix time.
    ! It returns the amount of days in the wanted year in the previous months.
  
    ! input
    integer, intent(in) :: month
    integer, intent(in) :: year
    ! output
    integer             :: month_day_vector
    
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
  
end module run_nml












