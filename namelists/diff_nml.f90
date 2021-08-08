! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

module diff_nml
	
	use definitions, only: wp
	
	implicit none
	
	logical  :: lklemp           ! turns the Klemp damping layer on or off
	real(wp) :: klemp_damp_max   ! the maximum Klemp damping coefficient
	real(wp) :: klemp_begin_rel  ! lower boundary of the Klemp damping layer in relation to TOA
	
	namelist /diff/lklemp,klemp_damp_max,klemp_begin_rel
	
	contains
	
	subroutine diff_nml_setup

		! local variables
		integer :: fileunit
		
		lklemp          = .true.
		klemp_damp_max  = 0.25_wp
		klemp_begin_rel = 0.53_wp
	
        ! Open and read Namelist file.
        open(action="read", file="namelist.nml", newunit=fileunit)
        read(nml=diff, unit=fileunit)
        
        close(fileunit)
	
	end subroutine diff_nml_setup
	
end module diff_nml
