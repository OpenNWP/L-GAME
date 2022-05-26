! This source file is part of the Limited-area GAME version (L-GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module phase_trans

  ! In this module, phase transition rates are being calculated.
  
  use run_nml,          only: nlins,ncols,nlays,dtime,wp
  use surface_nml,      only: lsfc_phase_trans
  use constants,        only: T_0,EPSILON_SECURITY
  use definitions,      only: t_state,t_diag,t_irrev,t_grid
  use dictionary,       only: saturation_pressure_over_ice,saturation_pressure_over_water, &
                              enhancement_factor_over_water,enhancement_factor_over_ice, &
                              specific_gas_constants,phase_trans_heat,rel_humidity
  
  implicit none
  
  private
  
  public :: calc_h2otracers_source_rates
  
  contains
  
  subroutine calc_h2otracers_source_rates(state,diag,irrev,grid)
  
    ! This subroutine calculates the phase transition rates.
  
    ! input arguments, output
    type(t_state), intent(in)    :: state ! the state with which to calculate the phase transition rates
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_irrev), intent(inout) :: irrev ! irreversible quantities
    type(t_grid),  intent(in)    :: grid  ! grid properties
  
    ! local variables
    real(wp) :: dtime_sat               ! time step for the saturation adjustment
    real(wp) :: solid_temperature       ! temperature of the cloud ice
    real(wp) :: liquid_temperature      ! temperature of the cloud water
    real(wp) :: max_cloud_water_content ! maximum cloud water content in (kg cloud)/(kg dry air).
    real(wp) :: diff_density            ! amount of water the air can still take
    real(wp) :: phase_trans_density     ! density of water that changes its phase
    real(wp) :: saturation_pressure     ! saturation pressure of the water vapour
    real(wp) :: water_vapour_pressure   ! actual water vapour pressure
    real(wp) :: dry_pressure            ! pressure of dry air
    real(wp) :: air_pressure            ! total air pressure
    real(wp) :: layer_thickness         ! thickness of the lowest layer
    real(wp) :: diff_density_sfc        ! density of water the air can still take above water surfaces
    real(wp) :: saturation_pressure_sfc ! saturation pressure at the surface
    integer  :: ji,jk,jl                ! loop variables
    
    max_cloud_water_content = 0.2e-3
  
    ! calculating the time step for the saturation adjustment
    dtime_sat = 2._wp*dtime
  
    !$OMP PARALLEL
    !$OMP DO PRIVATE(ji,jk,jl,solid_temperature,liquid_temperature,diff_density,phase_trans_density, &
    !$OMP saturation_pressure,water_vapour_pressure,dry_pressure,air_pressure,layer_thickness, &
    !$OMP diff_density_sfc,saturation_pressure_sfc)
    do ji=1,nlins
      do jk=1,ncols
        do jl=1,nlays
        
          ! Preparations
          ! ------------
          
          ! determining the temperature of the cloud ice
          if (state%rho(ji,jk,jl,3)<EPSILON_SECURITY) then
            solid_temperature = T_0
          else
            solid_temperature = diag%temperature_gas(ji,jk,jl)
          endif

          ! determining the temperature of the liquid cloud water
          if (state%rho(ji,jk,jl,4)<EPSILON_SECURITY) then
            liquid_temperature = T_0
          else
            liquid_temperature = diag%temperature_gas(ji,jk,jl)
          endif

          ! determining the saturation pressure
          ! "positive" temperatures (the saturation pressure is different over water compared to over ice)
          if (diag%temperature_gas(ji,jk,jl)>=T_0) then
            saturation_pressure = saturation_pressure_over_water(diag%temperature_gas(ji,jk,jl))
          ! "negative" temperatures
          else
            saturation_pressure = saturation_pressure_over_ice(diag%temperature_gas(ji,jk,jl))
          endif

          ! determining the water vapour pressure (using the EOS)
          water_vapour_pressure = state%rho(ji,jk,jl,6)*specific_gas_constants(1)*diag%temperature_gas(ji,jk,jl)
          
          ! determining the dry air pressure
          dry_pressure = state%rho(ji,jk,jl,5)*specific_gas_constants(0)*diag%temperature_gas(ji,jk,jl)
          
          ! calculating the total air pressure
          air_pressure = water_vapour_pressure + dry_pressure
          
          ! multiplying the saturation pressure by the enhancement factor
          if (diag%temperature_gas(ji,jk,jl)>=T_0) then
            saturation_pressure = enhancement_factor_over_water(air_pressure)*saturation_pressure
          ! "negative" temperatures
          else
            saturation_pressure = enhancement_factor_over_ice(air_pressure)*saturation_pressure
          endif

          ! the amount of water vapour that the air can still take 
          diff_density = (saturation_pressure - water_vapour_pressure)/(specific_gas_constants(1)*diag%temperature_gas(ji,jk,jl))

          ! Clouds
          ! ------
          
          ! the case where the air is not over-saturated
          if (diff_density>=0._wp) then
            ! temperature >= 0 °C
            if (diag%temperature_gas(ji,jk,jl)>=T_0) then
              ! It is assumed that the still present ice vanishes within one time step.
              irrev%mass_source_rates(ji,jk,jl,3) = -state%rho(ji,jk,jl,3)/dtime_sat

              !The amount of liquid water per volume that will evaporate.
              !In case the air cannot take all the water,not everything will evaporate.

              phase_trans_density = min(state%rho(ji,jk,jl,4),diff_density)

              ! The source rate for the liquid water consists of two terms:
              ! 1.) the evaporation
              ! 2.) the melting of ice

              irrev%mass_source_rates(ji,jk,jl,4) = (state%rho(ji,jk,jl,3) - phase_trans_density)/dtime_sat

              ! the tendency for the water vapour
              irrev%mass_source_rates(ji,jk,jl,5) = phase_trans_density/dtime_sat

              ! the heat source rates acting on the ice
              irrev%heat_source_rates(ji,jk,jl,3) = irrev%mass_source_rates(ji,jk,jl,3) &
              *phase_trans_heat(2,solid_temperature)

              ! the heat source rates acting on the liquid water
              irrev%heat_source_rates(ji,jk,jl,4) = &
              ! the evaporation
              -phase_trans_density*phase_trans_heat(0,T_0)/dtime_sat
            ! temperature<0 °C
            else
              ! Everything that can sublimate will sublimate.
              phase_trans_density = min(state%rho(ji,jk,jl,3),diff_density)

              ! the tendency for the ice contains two terms:
              ! 1.) the freezing
              ! 2.) the phase transition through sublimation

              irrev%mass_source_rates(ji,jk,jl,3) = (state%rho(ji,jk,jl,4) - phase_trans_density)/dtime_sat

              ! It is assumed that the still present liquid water vanishes within one time step.
              irrev%mass_source_rates(ji,jk,jl,4) = -state%rho(ji,jk,jl,4)/dtime_sat

              ! the tendency for the water vapour
              irrev%mass_source_rates(ji,jk,jl,5) = phase_trans_density/dtime_sat

              ! the heat source rates acting on the ice
              irrev%heat_source_rates(ji,jk,jl,3) = ( &
              ! the freezing
              state%rho(ji,jk,jl,4)*phase_trans_heat(2,solid_temperature) &
              ! the sublimation
              - phase_trans_density*phase_trans_heat(1,solid_temperature))/dtime_sat

              ! the heat source rates acting on the liquid water
              irrev%heat_source_rates(ji,jk,jl,4) = 0._wp
            endif
          ! the case where the air is over-saturated
          else
            ! the vanishing of water vapour through the phase transition
            irrev%mass_source_rates(ji,jk,jl,5) = diff_density/dtime_sat
            ! temperature >= 0._wp °C
            if (diag%temperature_gas(ji,jk,jl)>=T_0) then
              ! It is assumed that the still present ice vanishes within one time step.
              irrev%mass_source_rates(ji,jk,jl,3) = -state%rho(ji,jk,jl,3)/dtime_sat

              ! The source rate for the liquid water consists of two terms:
              ! 1.) the condensation
              ! 2.) the melting of ice

              irrev%mass_source_rates(ji,jk,jl,4) = (-diff_density + state%rho(ji,jk,jl,3))/dtime_sat

              ! the heat source rates acting on the ice
              irrev%heat_source_rates(ji,jk,jl,3) = -state%rho(ji,jk,jl,3)*phase_trans_heat(2,solid_temperature)/dtime_sat

              ! the heat source rates acting on the liquid water
              irrev%heat_source_rates(ji,jk,jl,4) = &
              ! it is only affected by the condensation
              -diff_density*phase_trans_heat(0,liquid_temperature)/dtime_sat
              ! temperature<0 °C
            else
              ! The source rate for the ice consists of two terms:
              ! 1.) the resublimation
              ! 2.) the melting of ice

              irrev%mass_source_rates(ji,jk,jl,3) = (-diff_density + state%rho(ji,jk,jl,4))/dtime_sat

              ! It is assumed that the liquid water disappears within one time step.
              irrev%mass_source_rates(ji,jk,jl,4) = -state%rho(ji,jk,jl,4)/dtime_sat

              ! the heat source rates acting on the ice
              irrev%heat_source_rates(ji,jk,jl,3) = &
              ! the component through the resublimation
              (-diff_density*phase_trans_heat(1,solid_temperature) &
              ! the component through freezing
              + state%rho(ji,jk,jl,4)*phase_trans_heat(2,solid_temperature))/dtime_sat

              ! the heat source rates acting on the liquid water
              irrev%heat_source_rates(ji,jk,jl,4) = 0._wp
            endif
          endif

          ! Precipitation
          ! -------------
          irrev%mass_source_rates(ji,jk,jl,1) = 0._wp
          irrev%mass_source_rates(ji,jk,jl,2) = 0._wp
          ! snow
          ! this only happens if the air is saturated
          if (diag%temperature_gas(ji,jk,jl)<T_0 .and. diff_density<=0._wp) then
            irrev%mass_source_rates(ji,jk,jl,1) = max(state%rho(ji,jk,jl,3) &
            - max_cloud_water_content*state%rho(ji,jk,jl,5),0._wp)/dtime_sat
            ! the snow creation comes at the cost of cloud ice particles
            irrev%mass_source_rates(ji,jk,jl,3) = irrev%mass_source_rates(ji,jk,jl,3) &
            - irrev%mass_source_rates(ji,jk,jl,1)
          endif
          ! rain
          if (diag%temperature_gas(ji,jk,jl)>=T_0 .and. diff_density<=0._wp) then
            ! this only happens if the air is saturated
            irrev%mass_source_rates(ji,jk,jl,2) = max(state%rho(ji,jk,jl,4) &
            - max_cloud_water_content*state%rho(ji,jk,jl,5),0._wp)/dtime_sat
            ! the rain creation comes at the cost of cloud water particles
            irrev%mass_source_rates(ji,jk,jl,4) = irrev%mass_source_rates(ji,jk,jl,4) &
            - irrev%mass_source_rates(ji,jk,jl,2)
          endif

          ! turning of snow to rain
          if (diag%temperature_gas(ji,jk,jl)>=T_0 .and. state%rho(ji,jk,jl,1)>0._wp) then
            irrev%mass_source_rates(ji,jk,jl,1) = -state%rho(ji,jk,jl,1)/dtime_sat
            irrev%mass_source_rates(ji,jk,jl,2) = irrev%mass_source_rates(ji,jk,jl,2) &
            - irrev%mass_source_rates(ji,jk,jl,1)
            irrev%heat_source_rates(ji,jk,jl,1) = irrev%mass_source_rates(ji,jk,jl,1)*phase_trans_heat(2,T_0)
          endif
          ! turning of rain to snow
          if (diag%temperature_gas(ji,jk,jl)<T_0 .and. state%rho(ji,jk,jl,2)>0._wp) then
            irrev%mass_source_rates(ji,jk,jl,2) = -state%rho(ji,jk,jl,2)/dtime_sat
            irrev%mass_source_rates(ji,jk,jl,1) = irrev%mass_source_rates(ji,jk,jl,1) &
            - irrev%mass_source_rates(ji,jk,jl,2)
            irrev%heat_source_rates(ji,jk,jl,1) = -irrev%mass_source_rates(ji,jk,jl,2)*phase_trans_heat(2,T_0)
          endif
          
        enddo
        
        ! Surface effects
        ! ---------------
        
        ! evaporation and latent heat rates
        if (lsfc_phase_trans .and. grid%is_land(ji,jk)==0) then
          ! saturation pressure at surface temperature
          if (state%temperature_soil(ji,jk,1)>=T_0) then
            saturation_pressure_sfc = saturation_pressure_over_water(state%temperature_soil(ji,jk,1))
            saturation_pressure_sfc = enhancement_factor_over_water(air_pressure)*saturation_pressure_sfc
          else
            saturation_pressure_sfc = saturation_pressure_over_ice(state%temperature_soil(ji,jk,1))
            saturation_pressure_sfc = enhancement_factor_over_ice(air_pressure)*saturation_pressure_sfc
          endif
          
          ! difference water vapour density between saturation at ground temperature and actual absolute humidity in the lowest model layer
          diff_density_sfc = saturation_pressure_sfc/(specific_gas_constants(1)*state%temperature_soil(ji,jk,1)) &
          - state%rho(ji,jk,nlays,6)

          ! the thickness of the lowest model layer (we need it as a result of Guass' theorem)
          layer_thickness = grid%z_w(ji,jk,nlays) - grid%z_w(ji,jk,nlays+1)

          ! evporation,sublimation
          irrev%mass_source_rates(ji,jk,nlays,5) = irrev%mass_source_rates(ji,jk,nlays,5) &
          + max(0._wp,diff_density_sfc/diag%scalar_flux_resistance(ji,jk))/layer_thickness

          ! calculating the latent heat flux density affecting the surface
          if (state%temperature_soil(ji,jk,1)>=T_0) then
            diag%power_flux_density_latent(ji,jk) = -phase_trans_heat(0,state%temperature_soil(ji,jk,1)) &
            *max(0._wp,diff_density_sfc/diag%scalar_flux_resistance(ji,jk))
          else
            diag%power_flux_density_latent(ji,jk) = -phase_trans_heat(1,state%temperature_soil(ji,jk,1)) &
            *max(0._wp,diff_density_sfc/diag%scalar_flux_resistance(ji,jk))
          endif
        endif
          
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
  
  end subroutine calc_h2otracers_source_rates

end module phase_trans











