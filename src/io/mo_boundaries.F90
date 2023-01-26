! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_boundaries
  
  ! This module handles everything dealing with boundary conditions.
  
  use mo_definitions,       only: t_state,t_bc,t_grid,wp
  use mo_run_nml,           only: ny,nx,n_layers,n_levels,t_init
  use mo_bc_nml,            only: n_swamp,bc_root_filename,bc_root_filename,dtime_bc,t_latest_bc
  use mo_constants,         only: M_PI,p_0,r_d,c_d_v
  use mo_constituents_nml,  only: n_condensed_constituents
  use mo_set_initial_state, only: read_from_nc
  
  implicit none
  
  contains
  
  subroutine update_boundaries(state,bc,t_now,grid)
    
    ! This subroutine brings the boundary conditions into the model.
    
    type(t_state), intent(inout) :: state ! state of the model (which will be modified)
    type(t_bc),    intent(inout) :: bc    ! boundary conditions
    real(wp),      intent(in)    :: t_now ! model time
    type(t_grid),  intent(in)    :: grid  ! grid properties (containing the background state)
    
    ! local variables
    real(wp) :: old_weight ! old time step interpolation weight
    real(wp) :: new_weight ! new time step interpolation weight
    integer  :: ji         ! horizontal index
    integer  :: jk         ! horizontal index
    integer  :: jl         ! layer index
    
    ! setting the time interpolation weights
    old_weight = 1._wp - (t_now-t_latest_bc)/dtime_bc
    ! reading from a new boundary conditions file if necessary
    if (old_weight<0._wp) then
      call read_boundaries(bc,t_latest_bc+dtime_bc,bc%index_old)
      ! swapping the indices
      bc%index_old = bc%index_new
      bc%index_new = 1
      if (bc%index_old==1) then
        bc%index_new = 2
      endif
      ! updating the latest boundary conditions read time
      t_latest_bc = t_latest_bc + dtime_bc
    endif
    new_weight = 1._wp - old_weight
    
    ! linear combination of the model state and the boundary conditions
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx
        do jl=1,n_layers
          state%rho(ji,jk,jl,:) = bc%scalar_bc_factor(ji,jk)*(old_weight*bc%rho(ji,jk,jl,:,bc%index_old) &
          + new_weight*bc%rho(ji,jk,jl,:,bc%index_new)) + (1._wp - bc%scalar_bc_factor(ji,jk))*state%rho(ji,jk,jl,:)
          state%rhotheta_v(ji,jk,jl) = bc%scalar_bc_factor(ji,jk)*(old_weight*bc%rhotheta_v(ji,jk,jl,bc%index_old) &
          + new_weight*bc%rhotheta_v(ji,jk,jl,bc%index_new)) + (1._wp - bc%scalar_bc_factor(ji,jk))*state%rhotheta_v(ji,jk,jl)
          state%wind_w(ji,jk,jl) = bc%scalar_bc_factor(ji,jk)*(old_weight*bc%wind_w(ji,jk,jl,bc%index_old) &
          + new_weight*bc%wind_w(ji,jk,jl,bc%index_new)) + (1._wp - bc%scalar_bc_factor(ji,jk))*state%wind_w(ji,jk,jl)
        enddo
        state%wind_w(ji,jk,n_levels) = bc%scalar_bc_factor(ji,jk)*(old_weight*bc%wind_w(ji,jk,n_levels,bc%index_old) &
        + new_weight*bc%wind_w(ji,jk,n_levels,bc%index_new)) + (1._wp - bc%scalar_bc_factor(ji,jk))*state%wind_w(ji,jk,n_levels)
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny
      do jk=1,nx+1
        do jl=1,n_layers
          state%wind_u(ji,jk,jl) = bc%u_bc_factor(ji,jk)*(old_weight*bc%wind_u(ji,jk,jl,bc%index_old) &
          + new_weight*bc%wind_u(ji,jk,jl,bc%index_new)) + (1._wp - bc%u_bc_factor(ji,jk))*state%wind_u(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji,jk,jl)
    do ji=1,ny+1
      do jk=1,nx
        do jl=1,n_layers
          state%wind_v(ji,jk,jl) = bc%v_bc_factor(ji,jk)*(old_weight*bc%wind_v(ji,jk,jl,bc%index_old) &
          + new_weight*bc%wind_v(ji,jk,jl,bc%index_new)) + (1._wp - bc%v_bc_factor(ji,jk))*state%wind_v(ji,jk,jl)
        enddo
      enddo
    enddo
    !$omp end parallel do
    
    ! setting the virtual potential temperature perturbation
    !$omp parallel workshare
    state%theta_v_pert = state%rhotheta_v/state%rho(:,:,:,n_condensed_constituents+1) - grid%theta_v_bg
    state%exner_pert = (r_d*state%rhotheta_v/p_0)**(r_d/c_d_v) - grid%exner_bg
    !$omp end parallel workshare
    
  end subroutine update_boundaries
  
  subroutine read_boundaries(bc,t_update,time_step_index)
    
    ! This subroutine reads the boundary conditions from a netCDF file.
    
    type(t_bc), intent(inout) :: bc              ! boundary conditions
    real(wp),   intent(in)    :: t_update        ! valid time of the boundary conditions
    integer,    intent(in)    :: time_step_index ! index of the boundary conditions time step (1 or 2) to which the
                                                 ! data in the netCDF file will be written
    
    ! local variables
    character(len=64) :: filename ! file to read the boundary state from
    
    ! constructing the filename to read the data from
    filename = "../../real_weather/" // trim(bc_root_filename) // "+" // &
    trim(int2string(int(t_update - t_init))) // "s.nc"
    
    write(*,*) "Reading boundary conditions from file", trim(filename), "..."
    
    ! reading the boundary conditions from a the netCDF file
    call read_from_nc(bc%rho(:,:,:,:,time_step_index),bc%rhotheta_v(:,:,:,time_step_index), &
    bc%wind_u(:,:,:,time_step_index),bc%wind_v(:,:,:,time_step_index),bc%wind_w(:,:,:,time_step_index),filename)
    
    write(*,*) "Boundary conditions read."
    
  end subroutine read_boundaries
  
  subroutine setup_bc_factor(bc)
    
    ! This subroutine calculates the boundary conditions rescale factors.
    ! It only needs to be called once (in the beginning).
    
    type(t_bc), intent(inout) :: bc ! boundary conditions type
    
    ! local variables
    real(wp) :: dist_from_boundary ! index distance from the boundary of the domain
    integer  :: ji                 ! horizontal index
    integer  :: jk                 ! horizontal index
    
    ! rescale factor for scalar fields
    !$omp parallel do private(ji,jk)
    do ji=1,ny
      do jk=1,nx
        dist_from_boundary = min(ji-1,jk-1,ny-ji,nx-jk)
        bc%scalar_bc_factor(ji,jk) = max(1._wp - dist_from_boundary/n_swamp,0._wp)
        bc%scalar_bc_factor(ji,jk) = sin(0.5_wp*M_PI*bc%scalar_bc_factor(ji,jk))**2
      enddo
    enddo
    !$omp end parallel do
    
    ! u rescale factor
    !$omp parallel do private(ji,jk)
    do ji=1,ny
      do jk=1,nx+1
        dist_from_boundary = min(ji-1,jk-1,ny-ji,nx+1-jk)
        bc%u_bc_factor(ji,jk) = max(1._wp - dist_from_boundary/n_swamp,0._wp)
        bc%u_bc_factor(ji,jk) = sin(0.5_wp*M_PI*bc%u_bc_factor(ji,jk))**2
      enddo
    enddo
    !$omp end parallel do
    
    ! v rescale factor
    !$omp parallel do private(ji,jk)
    do ji=1,ny+1
      do jk=1,nx
        dist_from_boundary = min(ji-1,jk-1,ny+1-ji,nx-jk)
        bc%v_bc_factor(ji,jk) = max(1._wp - dist_from_boundary/n_swamp,0._wp)
        bc%v_bc_factor(ji,jk) = sin(0.5_wp*M_PI*bc%v_bc_factor(ji,jk))**2
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine setup_bc_factor
  
  character(len=64) function int2string(input)
    
    ! This is a helper function which converts an integer to a string.
    
    integer, intent(in) :: input ! integer which will be converted to a string
    
    write(int2string,*) input
    int2string = adjustl(int2string)
    
  end function int2string

end module mo_boundaries








