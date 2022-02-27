! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module set_initial_state

  ! This module handles everything dealing with IO.

  use definitions,      only: t_state,wp,t_diag,t_grid
  use netcdf
  use run_nml,          only: nlins,ncols,nlays,scenario,run_id
  use constituents_nml, only: no_of_condensed_constituents,no_of_constituents
  use dictionary,       only: specific_gas_constants,spec_heat_capacities_p_gas
  use constants,        only: tropo_height,surface_temp,lapse_rate,inv_height,p_0, &
                              gravity,p_0_standard,re,t_grad_inv
  use io_nml,           only: restart_filename

  implicit none
  
  private
  
  public :: ideal
  public :: bg_temp
  public :: bg_pres
  public :: geopot
  public :: restart
  public :: nc_check
  
  contains
  
  subroutine ideal(state,diag,grid)
  
    ! sets the initial state of the model calculation in terms of analytic functions
    
    type(t_state), intent(inout) :: state ! state to write the initial state to
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! model grid
    
    ! local variables
    real(wp) :: pres_lowest_layer(nlins,ncols) ! pressure in the lowest layer
    real(wp) :: n_squared                      ! Brunt-Väisälä frequency for the Schär test case
    real(wp) :: gravity                        ! gravity acceleration
    real(wp) :: delta_z                        ! delta z
    real(wp) :: T_0                            ! MSLP temperature variable
    real(wp) :: r_d                            ! specific gas constant of dry air
    real(wp) :: c_p                            ! specific heat capacity at const. pressure of dry air
    integer  :: ji,jk,jl                       ! loop indices
    
    r_d = specific_gas_constants(0)
    c_p = spec_heat_capacities_p_gas(0)
      
    select case (trim(scenario))
    
      case("standard")
      
        ! This test case is the standard atmosphere.
      
        ! There is no horizontal wind in this case.
        state%wind_u = 0._wp
        state%wind_v = 0._wp
        
        !$OMP PARALLEL
        !$OMP DO PRIVATE(ji,jk,jl)
        do ji=1,nlins
          do jk=1,ncols
            do jl=1,nlays
              diag%scalar_placeholder(ji,jk,jl) = bg_temp(grid%z_scalar(ji,jk,jl))
            enddo
            pres_lowest_layer(ji,jk) = bg_pres(grid%z_scalar(ji,jk,nlays))
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

      case("schaer")
      
        ! Schär et al. (2001): A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models
        
        state%wind_u = 18.71_wp
        state%wind_v = 0._wp
       
        n_squared = (0.01871_wp)**2
        T_0 = 273.16_wp
       
        !$OMP PARALLEL
        !$OMP DO PRIVATE(ji,jk,jl,delta_z,gravity)
        do ji=1,nlins
          do jk=1,ncols
            ! potential temperature in the lowest layer
            ! calculating delta_z
            delta_z = grid%z_scalar(ji,jk,nlays-1)-grid%z_scalar(ji,jk,nlays)
            ! calculating the gravity
            gravity = (geopot(grid%z_scalar(ji,jk,nlays-1))-geopot(grid%z_scalar(ji,jk,nlays)))/delta_z
            ! result
            state%theta_pert(ji,jk,nlays) &
            ! value at MSL
            = T_0 &
            + T_0/gravity*n_squared*grid%z_scalar(ji,jk,nlays) &
            ! substracting the background state
            - grid%theta_bg(ji,jk,nlays)
            state%exner_pert(ji,jk,nlays)=(exp(-grid%z_scalar(ji,jk,nlays)/8000._wp))**(r_d/c_p) &
            - grid%exner_bg(ji,jk,nlays)
            ! stacking the potential temperature
            do jl=nlays-1,1,-1
              ! calculating delta_z
              delta_z = grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1)
              ! calculating the gravity
              gravity = (geopot(grid%z_scalar(ji,jk,jl))-geopot(grid%z_scalar(ji,jk,jl+1)))/delta_z
              ! result
              state%theta_pert(ji,jk,jl) &
              ! value in the lower layer
              = state%theta_pert(ji,jk,jl+1)+grid%theta_bg(ji,jk,jl+1) &
              + (state%theta_pert(ji,jk,jl+1)+grid%theta_bg(ji,jk,jl+1))/gravity*n_squared*delta_z &
              ! substracting the background state
              - grid%theta_bg(ji,jk,jl)
            enddo
            ! stacking the Exner pressure
            do jl=nlays-1,1,-1
              state%exner_pert(ji,jk,jl)=state%exner_pert(ji,jk,jl+1) &
              - (state%theta_pert(ji,jk,jl)+state%theta_pert(ji,jk,jl+1))/(grid%theta_bg(ji,jk,jl)+grid%theta_bg(ji,jk,jl+1) &
              +state%theta_pert(ji,jk,jl)+state%theta_pert(ji,jk,jl+1))*(grid%exner_bg(ji,jk,jl)-grid%exner_bg(ji,jk,jl+1))
            enddo
            ! filling up what's needed for the unessential_init routine
            ! temperature
            do jl=1,nlays
              diag%scalar_placeholder(ji,jk,jl)=(grid%exner_bg(ji,jk,jl)+state%exner_pert(ji,jk,jl)) &
              *(grid%theta_bg(ji,jk,jl)+state%theta_pert(ji,jk,jl))
            enddo
            ! pressure in the lowest layer
            pres_lowest_layer(ji,jk)=p_0*(grid%exner_bg(ji,jk,nlays)+state%exner_pert(ji,jk,nlays))**(c_p/r_d)
          enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    
    endselect
    
    ! humidity
    state%rho(:,:,:,no_of_constituents) = 0._wp
    
    ! condensates
    state%rho(:,:,:,1:no_of_condensed_constituents) = 0._wp
    state%condensed_rho_t = 0._wp
    
    
    do ji=1,nlins
      do jk=1,ncols
        state%temperature_soil(ji,jk,:) = 280._wp
      enddo
    enddo
    
    call unessential_init(state,diag,grid,pres_lowest_layer)
    
  end subroutine ideal
  
  subroutine restart(state,diag,grid)
  
    ! reads the initial state of the model calculation from a netcdf file
    
    ! input arguments and output
    type(t_state), intent(inout) :: state ! state to write the initial state to
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! model grid
    
    ! local variables
    integer           :: ncid                           ! ID of the NetCDF file
    character(len=64) :: filename                       ! filename containing the initial state
    integer           :: varid_rho                      ! variable ID of the densities
    integer           :: varid_t                        ! variable ID of the temperature
    integer           :: varid_u                        ! variable ID of the u-wind
    integer           :: varid_v                        ! variable ID of the v-wind
    real(wp)          :: pres_lowest_layer(nlins,ncols) ! pressure in the lowest layer
    
    filename = "../../real_weather/" // trim(restart_filename)
    
    ! creating the NetCDF file
    call nc_check(nf90_open(trim(restart_filename),NF90_CLOBBER,ncid))
    
    ! reading the variable IDs
    call nc_check(nf90_inq_varid(ncid,"rho",varid_rho))
    call nc_check(nf90_inq_varid(ncid,"T",varid_t))
    call nc_check(nf90_inq_varid(ncid,"u",varid_u))
    call nc_check(nf90_inq_varid(ncid,"v",varid_v))
    
    ! reading the NetCDF fields
    call nc_check(nf90_get_var(ncid,varid_rho,state%rho))
    call nc_check(nf90_get_var(ncid,varid_t,diag%scalar_placeholder))
    call nc_check(nf90_get_var(ncid,varid_u,state%wind_u))
    call nc_check(nf90_get_var(ncid,varid_v,state%wind_v))
    
    ! closing the NetCDF file
    call nc_check(nf90_close(ncid))
    
    ! setting the pressure in the lowest layer
    ! moist case
    if (no_of_constituents>1) then
      pres_lowest_layer &
      = state%rho(:,:,nlays,no_of_condensed_constituents+1)*specific_gas_constants(0)*diag%scalar_placeholder(:,:,nlays) &
      + state%rho(:,:,nlays,no_of_condensed_constituents+2)*specific_gas_constants(1)*diag%scalar_placeholder(:,:,nlays)
    ! dry case
    else
      pres_lowest_layer &
      = state%rho(:,:,nlays,1)*specific_gas_constants(0)*diag%scalar_placeholder(:,:,nlays)
    endif
    
    call unessential_init(state,diag,grid,pres_lowest_layer)
    
  end subroutine restart
  
  subroutine unessential_init(state,diag,grid,pres_lowest_layer)
  
    ! setting the unessential quantities of an initial state
    ! scalar_placeholder is the temperature here
    
    type(t_state), intent(inout) :: state                  ! state to work with
    type(t_diag),  intent(in)    :: diag                   ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid                   ! model grid
    real(wp),      intent(in)    :: pres_lowest_layer(:,:) ! pressure in the lowest layer
    
    ! local variables
    real(wp) :: b,c         ! abbreviations needed for the hydrostatic initialization routine
    real(wp) :: temperature ! single temperature value
    real(wp) :: pressure    ! single pressure value
    real(wp) :: r_d         ! specific gas constant of dry air
    real(wp) :: c_p         ! specific heat capacity at const. pressure of dry air
    integer  :: ji,jk,jl    ! loop indices
    
    r_d = specific_gas_constants(0)
    c_p = spec_heat_capacities_p_gas(0)
    
    ! integrating the hydrostatic initial state according to the given temperature field and pressure in the lowest layer
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,b,c,temperature,pressure)
    do ji=1,nlins
      do jk=1,ncols  
        ! integrating from bottom to top
        do jl=nlays,1,-1
          temperature = diag%scalar_placeholder(ji,jk,jl)
          ! lowest layer
          if (jl == nlays) then
            pressure = pres_lowest_layer(ji,jk)
            state%exner_pert(ji,jk,jl) = (pressure/p_0)**(r_d/c_p)
            state%theta_pert(ji,jk,jl) = temperature/state%exner_pert(ji,jk,jl)
          ! other layers
          else
            ! solving a quadratic equation for the Exner pressure
            b = -0.5_wp*state%exner_pert(ji,jk,jl+1)/bg_temp(grid%z_scalar(ji,jk,jl+1)) &
            *(temperature - bg_temp(grid%z_scalar(ji,jk,jl+1)) + 2.0_wp/ &
            c_p*(geopot(grid%z_scalar(ji,jk,jl)) - geopot(grid%z_scalar(ji,jk,jl+1))))
            c = state%exner_pert(ji,jk,jl+1)**2*temperature/bg_temp(grid%z_scalar(ji,jk,jl+1))
            state%exner_pert(ji,jk,jl) = b+sqrt(b**2+c)
            state%theta_pert(ji,jk,jl) = temperature/state%exner_pert(ji,jk,jl)
          endif
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    
    ! density
    state%rho(:,:,:,no_of_condensed_constituents+1) = p_0*(state%exner_pert)**(c_p/r_d) &
    /(r_d*diag%scalar_placeholder)
    ! potential temperature density
    state%rhotheta = state%rho(:,:,:,no_of_condensed_constituents+1)*state%theta_pert
    
    ! substracting the background state
    state%theta_pert = state%theta_pert - grid%theta_bg
    state%exner_pert = state%exner_pert - grid%exner_bg
    
    ! vertical wind velocity is always set to zero
    state%wind_w = 0._wp
    
  end subroutine unessential_init
  
  function bg_temp(height)
  
    ! This function returns the temperature of the background state.
    
    real(wp), intent(in) :: height
    ! output
    real(wp)             :: bg_temp

    ! troposphere
    if (height < tropo_height) then  
      bg_temp = surface_temp - lapse_rate*height
    ! constant temperature layer
    elseif (height < inv_height) then
      bg_temp = surface_temp - lapse_rate*tropo_height
    ! inversion
    else
      bg_temp = surface_temp - lapse_rate*tropo_height + t_grad_inv*(height - inv_height)
    endif
  
  end function bg_temp

  function bg_pres(height)
  
    ! This function returns the pressure of the background state (only used in the lowest layer during the initialization).
    
    real(wp), intent(in) :: height
    ! output
    real(wp)             :: bg_pres

    if (height<inv_height) then  
      bg_pres = p_0_standard*(1 - lapse_rate*height/surface_temp)**(gravity/(specific_gas_constants(0)*lapse_rate))
    elseif (height < tropo_height) then
      bg_pres = p_0_standard*(1 - lapse_rate*tropo_height/surface_temp)**(gravity/(specific_gas_constants(0)*lapse_rate)) &
      *exp(-gravity*(height - tropo_height)/(specific_gas_constants(0)*(surface_temp - lapse_rate*tropo_height)))
    else
      write(*,*) "Argument of bg_pres is above the inversion height. This is unrealistic in the lowest layer."
      write(*,*) "Aborting."
      call exit(1)
    endif
  
  end function bg_pres
  
  function geopot(height)
  
    ! This function returns the geopotential as a function of the geometrical height.
  
    ! input
    real(wp), intent(in) :: height
    ! output
    real(wp)             :: geopot
    
    geopot = -gravity*re**2/(re+height)+gravity*re
  
  end function geopot
  
  subroutine nc_check(i_status)
  
    ! This checks wether a NetCDF function threw an error.
  
    integer, intent(in) :: i_status

    if(i_status /= nf90_noerr) then 
      print *, trim(nf90_strerror(i_status))
      stop "Netcdf threw an error."
    end if
  end subroutine nc_check  

end module set_initial_state








