! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module thermodynamics

	! In this file, thermodynamic relations are calculated.

	use hetero_nml,  only: spec_heat_capacities_v_gas_lookup, spec_heat_capacities_p_gas_lookup, &
	                       specific_gas_constants_lookup
	use definitions, only: wp
	
	implicit none
	
	private
	
	public :: spec_heat_cap_diagnostics_v
	public :: spec_heat_cap_diagnostics_p
	public :: gas_constant_diagnostics
	
	contains

	function spec_heat_cap_diagnostics_v(j_constituents)
		
		integer, intent(in) :: j_constituents
		
		! specific heat capacity at constant volume
		real(wp) :: spec_heat_cap_diagnostics_v
		
		spec_heat_cap_diagnostics_v = spec_heat_capacities_p_gas_lookup(1)
		
	end function spec_heat_cap_diagnostics_v

	function spec_heat_cap_diagnostics_p(j_constituents)
		
		integer, intent(in) :: j_constituents
		
		! specific heat capacity at constant pressure
		real(wp) :: spec_heat_cap_diagnostics_p
		
		spec_heat_cap_diagnostics_p = spec_heat_capacities_p_gas_lookup(1)
		
	end function spec_heat_cap_diagnostics_p
	
	function gas_constant_diagnostics(j_constituents)
		
		integer, intent(in) :: j_constituents
		
		! specific heat capacity at constant pressure
		real(wp) :: gas_constant_diagnostics
		
		gas_constant_diagnostics = specific_gas_constants_lookup(1)
		
	end function gas_constant_diagnostics

end module thermodynamics











