! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module set_initial_state

  ! This module handles everything dealing with IO.

  use definitions,      only: t_state,wp,t_diag,t_grid
  use netcdf
  use run_nml,          only: ny,nx,nlays,scenario,run_id,lplane,dx
  use constituents_nml, only: n_condensed_constituents,n_constituents
  use constants,        only: tropo_height,surface_temp,lapse_rate,inv_height,p_0, &
                              gravity,p_0_standard,r_e,t_grad_inv,M_PI,r_d,c_d_p,c_d_v
  use io_nml,           only: restart_filename

  implicit none
  
  contains
  
  subroutine ideal_init(state,diag,grid)
  
    ! sets the initial state of the model calculation in terms of analytic functions
    
    type(t_state), intent(inout) :: state ! state to write the initial state to
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! model grid
    
    ! local variables
    real(wp) :: pres_lowest_layer(ny,nx)              ! pressure in the lowest layer
    real(wp) :: n_squared                                   ! Brunt-Väisälä frequency for the Schär test case
    real(wp) :: gravity_local                               ! gravity acceleration
    real(wp) :: delta_z                                     ! delta z
    real(wp) :: T_0                                         ! MSLP temperature variable
    integer  :: ji,jk,jl                                    ! loop indices
    real(wp) :: u_0,z_1,z_2,r,rho_0,x_0,z_0,A_x,A_Z,x_coord ! variables needed for the advection test
      
    select case (trim(scenario))
    
      case("standard","seabreeze")
      
        ! This test case is the standard atmosphere.
      
        ! There is no horizontal wind in this case
        !$omp parallel workshare
        state%wind_u = 0._wp
        state%wind_v = 0._wp
        !$omp end parallel workshare
        
        !$omp parallel do private(ji,jk,jl)
        do ji=1,ny
          do jk=1,nx
            do jl=1,nlays
              diag%scalar_placeholder(ji,jk,jl) = bg_temp(grid%z_scalar(ji,jk,jl))
            enddo
            pres_lowest_layer(ji,jk) = bg_pres(grid%z_scalar(ji,jk,nlays))
          enddo
        enddo
        !$omp end parallel do
    
        !$omp parallel workshare
        ! humidity
        state%rho(:,:,:,n_constituents) = 0._wp
        ! condensates
        state%rho(:,:,:,1:n_condensed_constituents) = 0._wp
        !$omp end parallel workshare

      case("advection")
      
        ! This test case is the advection test from the SLEVE paper
        ! (Schär et al. (2001): A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models).
        
        u_0 = 10._wp
        z_1 = 4e3_wp
        z_2 = 5e3_wp
        A_x = 25e3_wp
        A_z = 3e3_wp
        x_0 = -50e3_wp
        z_0 = 9e3_wp
        rho_0 = 1e-3_wp
        
        ! horizontal wind
        !$omp parallel do private(ji,jk,jl)
        do ji=1,ny
          do jk=1,nx+1
            do jl=1,nlays
              if (grid%z_u(ji,jk,jl)>=z_2) then
                state%wind_u(ji,jk,jl) = u_0
              elseif (z_1<=grid%z_u(ji,jk,jl) .and. grid%z_u(ji,jk,jl)<=z_2) then
                state%wind_u(ji,jk,jl) = u_0*sin(0.5_wp*M_PI*(grid%z_u(ji,jk,jl)-z_1)/(z_2-z_1))**2
              elseif (grid%z_u(ji,jk,jl)<=z_1) then
                state%wind_u(ji,jk,jl) = 0._wp
              endif
            enddo
          enddo
        enddo
        !$omp end parallel do
        
        !$omp parallel workshare
        state%wind_v = 0._wp
        !$omp end parallel workshare
        
        !$omp parallel do private(ji,jk,jl)
        do ji=1,ny
          do jk=1,nx
            do jl=1,nlays
              diag%scalar_placeholder(ji,jk,jl) = bg_temp(grid%z_scalar(ji,jk,jl))
            enddo
            pres_lowest_layer(ji,jk) = bg_pres(grid%z_scalar(ji,jk,nlays))
          enddo
        enddo
        !$omp end parallel do
        
        ! density anomaly
        !$omp parallel do private(ji,jk,jl,r,x_coord)
        do ji=1,ny
          do jk=1,nx
            do jl=1,nlays
              x_coord = dx*jk - (nx/2 + 1)*dx
              r = (((x_coord-x_0)/A_x)**2 + ((grid%z_scalar(ji,jk,jl)-z_0)/A_z)**2)**0.5_wp
              if (r<=1._wp) then
                state%rho(ji,jk,jl,n_condensed_constituents+2) = rho_0*cos(0.5_wp*M_PI*r)**2
              else
                state%rho(ji,jk,jl,n_condensed_constituents+2) = 0._wp
              endif
            enddo
          enddo
        enddo
        !$omp end parallel do
    
        !$omp parallel workshare
        ! condensates
        state%rho(:,:,:,1:n_condensed_constituents) = 0._wp
        !$omp end parallel workshare

      case("schaer")
      
        ! Schär et al. (2001): A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models
        
        !$omp parallel workshare
        state%wind_u = 10._wp
        state%wind_v = 0._wp
        !$omp end parallel workshare
       
        n_squared = (0.01_wp)**2
        T_0 = 288._wp
       
        ! background state not yet substracted here
       
        !$omp parallel do private(ji,jk,jl,delta_z,gravity_local)
        do ji=1,ny
          do jk=1,nx
            ! calculating delta_z
            delta_z = grid%z_scalar(ji,jk,nlays)
            ! calculating the gravity
            gravity_local = geopot(grid%z_scalar(ji,jk,nlays))/delta_z
            
            ! virtual potential temperature in the lowest layer
            state%theta_v_pert(ji,jk,nlays) = T_0 &
            *(1._wp + n_squared*delta_z/(2._wp*gravity_local)) &
            /(1._wp - n_squared*delta_z/(2._wp*gravity_local))
            ! Exner pressure in the lowest layer
            state%exner_pert(ji,jk,nlays) = 1._wp &
            - 2._wp*geopot(grid%z_scalar(ji,jk,nlays))/(c_d_p*(state%theta_v_pert(ji,jk,nlays) + T_0))
            
            ! stacking the virtual potential temperature
            do jl=nlays-1,1,-1
              ! calculating delta_z
              delta_z = grid%z_scalar(ji,jk,jl)-grid%z_scalar(ji,jk,jl+1)
              ! calculating the gravity
              gravity_local = (geopot(grid%z_scalar(ji,jk,jl))-geopot(grid%z_scalar(ji,jk,jl+1)))/delta_z
              
              ! result
              state%theta_v_pert(ji,jk,jl) = state%theta_v_pert(ji,jk,jl+1) &
              *(1._wp + n_squared*delta_z/(2._wp*gravity_local)) &
              /(1._wp - n_squared*delta_z/(2._wp*gravity_local))
            enddo
            
            ! stacking the Exner pressure
            do jl=nlays-1,1,-1
              ! result
              state%exner_pert(ji,jk,jl) = state%exner_pert(ji,jk,jl+1) &
              + 2._wp*(geopot(grid%z_scalar(ji,jk,jl+1)) - geopot(grid%z_scalar(ji,jk,jl))) &
              /(c_d_p*(state%theta_v_pert(ji,jk,jl) + state%theta_v_pert(ji,jk,jl+1)))
            enddo
            
            ! filling up what's needed for the unessential_ideal_init routine
            ! temperature
            do jl=1,nlays
              diag%scalar_placeholder(ji,jk,jl)=state%exner_pert(ji,jk,jl)*state%theta_v_pert(ji,jk,jl)
            enddo
            ! pressure in the lowest layer
            pres_lowest_layer(ji,jk)=p_0*state%exner_pert(ji,jk,nlays)**(c_d_p/r_d)
          enddo
        enddo
        !$omp end parallel do
    
        !$omp parallel workshare
        ! humidity
        state%rho(:,:,:,n_constituents) = 0._wp
        ! condensates
        state%rho(:,:,:,1:n_condensed_constituents) = 0._wp
        !$omp end parallel workshare
        
    endselect
    
    !$omp parallel do private(ji,jk)
    do ji=1,ny
      do jk=1,nx
        state%temperature_soil(ji,jk,:) = 280._wp
      enddo
    enddo
    !$omp end parallel do
    
    call unessential_ideal_init(state,diag,grid,pres_lowest_layer)
    
  end subroutine ideal_init
  
  subroutine restart(state,grid)
  
    ! This subroutine sets the initial state of a NWP run.
    
    type(t_state), intent(inout) :: state ! state to write the arrays to
    type(t_grid), intent(in)     :: grid  ! grid properties
  
    ! local variables
    character(len=64) :: filename ! file to read the initial state from
    
    filename = "../../real_weather/" // trim(restart_filename)
    
    call read_from_nc(state%rho,state%rhotheta_v,state%wind_u,state%wind_v,state%wind_w,filename)
    
    ! setting the virtual potential temperature perturbation
    !$omp parallel workshare
    state%theta_v_pert = state%rhotheta_v/state%rho(:,:,:,n_condensed_constituents+1) - grid%theta_v_bg
    state%exner_pert = (r_d*state%rhotheta_v/p_0)**(r_d/c_d_v) - grid%exner_bg
    !$omp end parallel workshare
  
  end subroutine restart
  
  subroutine read_from_nc(rho,rhotheta_v,wind_u,wind_v,wind_w,filename)
  
    ! This subroutine reads a model state from a NetCDF file.
    
    ! input arguments and output
    real(wp),          intent(out) :: rho(:,:,:,:)    ! mass densities
    real(wp),          intent(out) :: rhotheta_v(:,:,:) ! virtual potential temperature density
    real(wp),          intent(out) :: wind_u(:,:,:)   ! u-wind
    real(wp),          intent(out) :: wind_v(:,:,:)   ! v-wind
    real(wp),          intent(out) :: wind_w(:,:,:)   ! w-wind
    character(len=64), intent(in)  :: filename        ! filename to read from
    
    ! local variables
    integer :: ncid           ! ID of the NetCDF file
    integer :: varid_rho      ! variable ID of the densities
    integer :: varid_rhotheta_v ! variable ID of the virtual potential temperature density
    integer :: varid_u        ! variable ID of the u-wind
    integer :: varid_v        ! variable ID of the v-wind
    integer :: varid_w        ! variable ID of the w-wind
    
    ! opening the NetCDF file
    call nc_check(nf90_open(trim(filename),NF90_CLOBBER,ncid))
    
    ! reading the variable IDs
    call nc_check(nf90_inq_varid(ncid,"rho",varid_rho))
    call nc_check(nf90_inq_varid(ncid,"rhotheta_v",varid_rhotheta_v))
    call nc_check(nf90_inq_varid(ncid,"u",varid_u))
    call nc_check(nf90_inq_varid(ncid,"v",varid_v))
    call nc_check(nf90_inq_varid(ncid,"w",varid_w))
    
    ! reading the NetCDF fields
    call nc_check(nf90_get_var(ncid,varid_rho,rho))
    call nc_check(nf90_get_var(ncid,varid_rhotheta_v,rhotheta_v))
    call nc_check(nf90_get_var(ncid,varid_u,wind_u))
    call nc_check(nf90_get_var(ncid,varid_v,wind_v))
    call nc_check(nf90_get_var(ncid,varid_w,wind_w))
    
    ! closing the NetCDF file
    call nc_check(nf90_close(ncid))
    
  end subroutine read_from_nc
  
  subroutine unessential_ideal_init(state,diag,grid,pres_lowest_layer)
  
    ! setting the unessential quantities of an ideal initial state
    ! scalar_placeholder is the temperature here
    
    type(t_state), intent(inout) :: state                  ! state to work with
    type(t_diag),  intent(in)    :: diag                   ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid                   ! model grid
    real(wp),      intent(in)    :: pres_lowest_layer(:,:) ! pressure in the lowest layer
    
    ! local variables
    real(wp) :: b,c         ! abbreviations needed for the hydrostatic initialization routine
    real(wp) :: pressure    ! single pressure value
    integer  :: ji,jk,jl    ! loop indices
    
    ! integrating the hydrostatic initial state according to the given temperature field and pressure in the lowest layer
    !$omp parallel do private(ji,jk,jl,b,c,pressure)
    do ji=1,ny
      do jk=1,nx  
        ! integrating from bottom to top
        do jl=nlays,1,-1
          ! lowest layer
          if (jl==nlays) then
            pressure = pres_lowest_layer(ji,jk)
            state%exner_pert(ji,jk,jl) = (pressure/p_0)**(r_d/c_d_p)
          ! other layers
          else
            ! solving a quadratic equation for the Exner pressure
            b = -0.5_wp*state%exner_pert(ji,jk,jl+1)/diag%scalar_placeholder(ji,jk,jl+1) &
            *(diag%scalar_placeholder(ji,jk,jl) - diag%scalar_placeholder(ji,jk,jl+1) + 2.0_wp/ &
            c_d_p*(geopot(grid%z_scalar(ji,jk,jl)) - geopot(grid%z_scalar(ji,jk,jl+1))))
            c = state%exner_pert(ji,jk,jl+1)**2*diag%scalar_placeholder(ji,jk,jl)/diag%scalar_placeholder(ji,jk,jl+1)
            state%exner_pert(ji,jk,jl) = b+sqrt(b**2+c)
          endif
          ! this is the full virtual potential temperature here
          state%theta_v_pert(ji,jk,jl) = diag%scalar_placeholder(ji,jk,jl)/state%exner_pert(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel workshare
    ! density
    state%rho(:,:,:,n_condensed_constituents+1) = p_0*(state%exner_pert)**(c_d_p/r_d) &
    /(r_d*diag%scalar_placeholder)
    ! virtual potential temperature density
    state%rhotheta_v = state%rho(:,:,:,n_condensed_constituents+1)*state%theta_v_pert
    ! substracting the background state
    state%theta_v_pert = state%theta_v_pert - grid%theta_v_bg
    state%exner_pert = state%exner_pert - grid%exner_bg
    !$omp end parallel workshare
    
    ! vertical wind velocity is always set to zero
    !$omp parallel workshare
    state%wind_w = 0._wp
    !$omp end parallel workshare
    
  end subroutine unessential_ideal_init
  
  function bg_temp(height)
  
    ! This function returns the temperature of the background state.
    
    real(wp), intent(in) :: height  ! geometric height above mean sea level
    ! output
    real(wp)             :: bg_temp ! the result

    ! troposphere
    if (height<tropo_height) then  
      bg_temp = surface_temp - lapse_rate*height
    ! constant temperature layer
    elseif (height<inv_height) then
      bg_temp = surface_temp - lapse_rate*tropo_height
    ! inversion
    else
      bg_temp = surface_temp - lapse_rate*tropo_height + t_grad_inv*(height - inv_height)
    endif
  
  end function bg_temp

  function bg_pres(height)
  
    ! This function returns the pressure of the background state (only used in the lowest layer during the initialization).
    
    real(wp), intent(in) :: height  ! geomteric height above mean sea level
    ! output
    real(wp)             :: bg_pres ! the result

    if (height<inv_height) then  
      bg_pres = p_0_standard*(1 - lapse_rate*height/surface_temp)**(gravity/(r_d*lapse_rate))
    elseif (height<tropo_height) then
      bg_pres = p_0_standard*(1 - lapse_rate*tropo_height/surface_temp)**(gravity/(r_d*lapse_rate)) &
      *exp(-gravity*(height - tropo_height)/(r_d*(surface_temp - lapse_rate*tropo_height)))
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
    
    geopot = -gravity*r_e**2/(r_e+height)+gravity*r_e
    
    ! Schaer test case
    if (lplane) then
      geopot = gravity*height
    endif
  
  end function geopot
  
  subroutine nc_check(i_status)
  
    ! This checks wether a NetCDF function threw an error.
  
    integer, intent(in) :: i_status

    if(i_status/=nf90_noerr) then 
      print *, trim(nf90_strerror(i_status))
      stop "Netcdf threw an error."
    end if
  end subroutine nc_check  

end module set_initial_state








