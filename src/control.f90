! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

program control

	! This controls the model run from the beginning to the end.

	! reading in all namelists so that we know what we have to do
	call read_all_namelists()

	! firstly, the grid generator needs to be called to calculate the grid properties
	call grid_generator()

	! reading the initial state
	call read_init()

	! the loop over the time steps
	cont_steppint = .true.
	do while (cont_steppint)
		
		
		
	enddo
  
end program control



