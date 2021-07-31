! This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/RFPET

use grad

program control

  ! firstly, the grid generator needs to be called to calculate the grid properties
  call grid_generator()
  
  
  ! the loop over the time steps
  while
  
end program control
