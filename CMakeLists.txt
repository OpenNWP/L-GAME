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
                 # insert the root directory of RTE+RRTMGP here
set(RTERRTMGPDIR /home/max/code/rte-rrtmgp)
add_executable(
rfpet
src/definitions/definitions.f90
src/control.f90
src/io/io.f90
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
namelists/bc_nml.f90
namelists/diff_nml.f90
namelists/hetero_nml.f90
namelists/init_nml.f90
namelists/rad_nml.f90
namelists/run_nml.f90
namelists/surface_nml.f90
${RTERRTMGPDIR}/rrtmgp/mo_gas_optics.F90
${RTERRTMGPDIR}/rrtmgp/mo_gas_optics_rrtmgp.F90
${RTERRTMGPDIR}/rrtmgp/mo_rrtmgp_constants.F90
${RTERRTMGPDIR}/rrtmgp/mo_rrtmgp_util_string.F90
${RTERRTMGPDIR}/rrtmgp/mo_rrtmgp_util_reorder.F90
${RTERRTMGPDIR}/rrtmgp/mo_gas_concentrations.F90
${RTERRTMGPDIR}/rrtmgp/kernels/mo_gas_optics_kernels.F90
${RTERRTMGPDIR}/rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90
${RTERRTMGPDIR}/rte/mo_rte_lw.F90
${RTERRTMGPDIR}/rte/mo_rte_sw.F90
${RTERRTMGPDIR}/rte/mo_fluxes.F90
${RTERRTMGPDIR}/rte/mo_optical_props.F90
${RTERRTMGPDIR}/rte/mo_rte_config.F90
${RTERRTMGPDIR}/rte/mo_rte_kind.F90
${RTERRTMGPDIR}/rte/mo_source_functions.F90
${RTERRTMGPDIR}/rte/mo_rte_util_array.F90
${RTERRTMGPDIR}/rte/kernels/mo_fluxes_broadband_kernels.F90
${RTERRTMGPDIR}/rte/kernels/mo_optical_props_kernels.F90
${RTERRTMGPDIR}/rte/kernels/mo_rte_solver_kernels.F90
${RTERRTMGPDIR}/extensions/mo_fluxes_byband.F90
${RTERRTMGPDIR}/extensions/mo_fluxes_byband_kernels.F90
${RTERRTMGPDIR}/examples/mo_load_coefficients.F90
${RTERRTMGPDIR}/examples/mo_simple_netcdf.F90)

find_package(OpenMP)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O1 -Wall")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS} -lgfortran -O1 -Wall -Wno-tabs -Wno-c-binding-type -Wno-unused-dummy-argument")
                                                # insert output of 'nf-config --fflags --flibs' here
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdff -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,-z,now -lnetcdf -lnetcdf -ldl -lm")
target_link_libraries(rfpet m netcdf netcdff atmostracers)









