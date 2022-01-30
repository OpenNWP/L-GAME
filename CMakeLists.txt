# This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
# Github repository: https://github.com/OpenNWP/L-GAME

cmake_minimum_required(VERSION 3.10)
set(version_major 0)
set(version_minor 1)
set(version_patch 1)
set(version ${version_major}.${version_minor}.${version_patch})
project(lgame VERSION ${version})
enable_language(Fortran)
add_executable(
lgame
src/definitions.f90
src/control.f90
src/io.f90
src/boundaries.f90
src/grid_generator.f90
src/time_stepping/column_solvers.f90
src/time_stepping/explicit_vector_tendencies.f90
src/time_stepping/explicit_scalar_tendencies.f90
src/time_stepping/manage_rkhevi.f90
src/time_stepping/pressure_gradient.f90
src/thermodynamics.f90
src/subgrid_scale/effective_diff_coeffs.f90
src/spatial_operators/vorticity_flux.f90
src/spatial_operators/multiplications.f90
src/spatial_operators/linear_combine_two_states.f90
src/spatial_operators/gradient_operators.f90
src/spatial_operators/averaging.f90
src/spatial_operators/vorticities.f90
src/spatial_operators/momentum_diff_diss.f90
src/spatial_operators/inner_product.f90
src/spatial_operators/divergence_operators.f90
src/rterrtmgp_coupler.f90
src/namelists/bc_nml.f90
src/namelists/diff_nml.f90
src/namelists/hetero_nml.f90
src/namelists/init_nml.f90
src/namelists/rad_nml.f90
src/namelists/run_nml.f90
src/namelists/io_nml.f90
src/namelists/surface_nml.f90
../rte-rrtmgp/rrtmgp/mo_gas_optics.F90
../rte-rrtmgp/rrtmgp/mo_gas_optics_rrtmgp.F90
../rte-rrtmgp/rrtmgp/mo_rrtmgp_constants.F90
../rte-rrtmgp/rrtmgp/mo_rrtmgp_util_string.F90
../rte-rrtmgp/rrtmgp/mo_rrtmgp_util_reorder.F90
../rte-rrtmgp/rrtmgp/mo_gas_concentrations.F90
../rte-rrtmgp/rrtmgp/kernels/mo_gas_optics_kernels.F90
../rte-rrtmgp/rrtmgp/kernels/mo_rrtmgp_util_reorder_kernels.F90
../rte-rrtmgp/rte/mo_rte_lw.F90
../rte-rrtmgp/rte/mo_rte_sw.F90
../rte-rrtmgp/rte/mo_fluxes.F90
../rte-rrtmgp/rte/mo_optical_props.F90
../rte-rrtmgp/rte/mo_rte_config.F90
../rte-rrtmgp/rte/mo_rte_kind.F90
../rte-rrtmgp/rte/mo_source_functions.F90
../rte-rrtmgp/rte/mo_rte_util_array.F90
../rte-rrtmgp/rte/kernels/mo_fluxes_broadband_kernels.F90
../rte-rrtmgp/rte/kernels/mo_optical_props_kernels.F90
../rte-rrtmgp/rte/kernels/mo_rte_solver_kernels.F90
../rte-rrtmgp/extensions/mo_fluxes_byband.F90
../rte-rrtmgp/extensions/cloud_optics/mo_cloud_optics.F90
../rte-rrtmgp/examples/mo_load_coefficients.F90
../rte-rrtmgp/examples/mo_simple_netcdf.F90
../rte-rrtmgp/examples/all-sky/mo_load_cloud_coefficients.F90
)
find_package(OpenMP)
SET(CMAKE_Fortran_FLAGS "-O2 ${OpenMP_Fortran_FLAGS} -Wall -Wno-c-binding-type -Wno-unused-dummy-argument -I/usr/include -L/usr/lib/x86_64-linux-gnu -lnetcdff -lnetcdf")
target_link_libraries(lgame netcdff atmostracers geos95)









