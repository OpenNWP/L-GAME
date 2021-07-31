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
src/control.f
src/waves_coupler/waves_coupler.f
src/time_stepping/vertical_slice_solvers.f
src/time_stepping/vector_tendencies_expl.f
src/time_stepping/scalar_tendencies_expl.f
src/time_stepping/manage_rkhevi.f
src/time_stepping/manage_pressure_gradient.f
src/thermodynamics/thermodynamics.f
src/subgrid_scale/gravity_waves.f
src/subgrid_scale/effective_diff_coeffs.f
src/spatial_operators/vorticity_flux.f
src/spatial_operators/multiplications.f
src/spatial_operators/linear_combine_two_states.f
src/spatial_operators/gradient_operators.f
src/spatial_operators/averaging.f
src/spatial_operators/vorticities.f
src/spatial_operators/momentum_diff_diss.f
src/spatial_operators/inner_product.f
src/spatial_operators/divergences.f
src/soil_coupler/soil_coupler.f
src/radiation/rterrtmgp_coupler.f
src/ocean_coupler/ocean_coupler.f
src/io/io.f
)

find_package(OpenMP)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O1 -Wall")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS} -lgfortran -O1 -Wall -Wno-c-binding-type -Wno-unused-dummy-argument")
                                                # insert output of 'nf-config --fflags --flibs' here
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm")
target_link_libraries(rfpet m netcdf netcdff atmostracers)
install(TARGETS rfpet DESTINATION rfpet)









