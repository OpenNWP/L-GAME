! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

program control

	! This controls the model run from the beginning to the end.

	use run_nml,        only: run_nml_setup
	use grid_generator, only: grid_setup
	
	implicit none
	
	logical :: cont_steppint

	! reading in all namelists so that we know what we have to do
	print *, "Reading in run namelist ..."
	call run_nml_setup
	print *, "... run namelist read."

	! firstly, the grid generator needs to be called to calculate the grid properties
	print *, "Setting up the grid ..."
	call grid_setup
	print *, "... grid setup."

	! reading the initial state
	print *, "Reading the initial state..."
	call read_init()
	print *, "... initial state read."

	! the loop over the time steps
	cont_steppint = .true.
	do while (cont_steppint)
		
		
		
	enddo
  
end program control



