! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module hetero_nml

	use definitions, only: wp

	implicit none
	
	integer :: no_of_gaseous_constituents ! number of constituents of the gas phase
	
	namelist /hetero/no_of_gaseous_constituents

	! interface to C functions
	interface
		real(C_DOUBLE) function specific_gas_constants_lookup(gas_number) bind(c, name = "specific_gas_constants_lookup")
			  use, intrinsic::iso_c_binding
			  implicit none
			  integer(C_INT), value :: gas_number
		end function specific_gas_constants_lookup
	end interface
	interface
		real(C_DOUBLE) function spec_heat_capacities_v_gas_lookup(gas_number) bind(c, name = "spec_heat_capacities_v_gas_lookup")
			  use, intrinsic::iso_c_binding
			  implicit none
			  integer(C_INT), value :: gas_number
		end function spec_heat_capacities_v_gas_lookup
	end interface
	interface
		real(C_DOUBLE) function spec_heat_capacities_p_gas_lookup(gas_number) bind(c, name = "spec_heat_capacities_p_gas_lookup")
			  use, intrinsic::iso_c_binding
			  implicit none
			  integer(C_INT), value :: gas_number
		end function spec_heat_capacities_p_gas_lookup
	end interface

	contains

	subroutine hetero_nml_setup
	
		no_of_gaseous_constituents = 1
		
	end subroutine hetero_nml_setup

	function get_gas_contituents_ids(gas_constituent_id)
	
		integer, intent(in)  :: gas_constituent_id
		integer              :: get_gas_contituents_ids
		! local variable
		integer             :: ji
		integer gas_constituent_ids_vector(no_of_gaseous_constituents)
		
		do ji=1,no_of_gaseous_constituents
			gas_constituent_ids_vector(ji) = ji;
		enddo
		
		get_gas_contituents_ids = gas_constituent_ids_vector(gas_constituent_id)
		
	end function get_gas_contituents_ids

	function specific_gas_constants(gas_constituent_id)

		integer, intent(in)   :: gas_constituent_id
		real(wp)              :: specific_gas_constants
		
		specific_gas_constants = specific_gas_constants_lookup(get_gas_contituents_ids(gas_constituent_id))
		
	end function specific_gas_constants
	
end module hetero_nml






