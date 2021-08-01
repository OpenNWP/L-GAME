# This source file is part of the Regional Forecasting with Poisson-brackets in Exner-Theta formulation (RFPET), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/RFPET

cmake_minimum_required(VERSION 3.10)
set(version_major 0)
set(version_minor 1)
set(version_patch 1)
set(version ${version_major}.${version_minor}.${version_patch})
project(rfpet VERSION ${version})
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_Fortran_COMPILER mpifort)
enable_language(Fortran)

add_executable(
rfpet
src/control.f90
src/grid_generator/grid_generator.f90
src/waves_coupler/waves_coupler.f90
src/time_stepping/vertical_slice_solvers.f90
src/time_stepping/vector_tendencies_expl.f90
src/time_stepping/scalar_tendencies_expl.f90
src/time_stepping/manage_rkhevi.f90
src/time_stepping/manage_pressure_gradient.f90
src/thermodynamics/thermodynamics.f90
src/subgrid_scale/gravity_waves.f90
src/subgrid_scale/effective_diff_coeffs.f90
src/spatial_operators/vorticity_flux.f90
src/spatial_operators/multiplications.f90
src/spatial_operators/linear_combine_two_states.f90
src/spatial_operators/gradient_operators.f90
src/spatial_operators/averaging.f90
src/spatial_operators/vorticities.f90
src/spatial_operators/momentum_diff_diss.f90
src/spatial_operators/inner_product.f90
src/spatial_operators/divergences.f90
src/soil_coupler/soil_coupler.f90
src/radiation/rterrtmgp_coupler.f90
src/ocean_coupler/ocean_coupler.f90
src/io/io.f90
)

find_package(OpenMP)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O1 -Wall")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS} -lgfortran -O1 -Wall -Wno-tabs -Wno-c-binding-type -Wno-unused-dummy-argument")
                                                # insert output of 'nf-config --fflags --flibs' here
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm")
target_link_libraries(rfpet m netcdf netcdff atmostracers)
install(TARGETS rfpet DESTINATION rfpet)









