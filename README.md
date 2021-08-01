## Regional Forecasting with Poisson brackets in Exner-Theta formulation

RFPET is a numerical weather prediction (NWP) model. Properties:

* regional
* non-hydrostatic
* Eulerian
* (rotated) latitude-longitude grid
* C-grid
* using a hybrid of finite volume and finite difference methods
* time stepping: two-time-level Runge-Kutta scheme, modified into a HEVI (horizontally explicit, vertically implicit) and forward-backward scheme for stability, horizontal pressure gradient extrapolated and kept constant
* radiation: coupled to [RTE+RRTMGP](https://github.com/earth-system-radiation/rte-rrtmgp)
* using the Poisson bracket formulation by Gassmann and Herzog (2008) and Gassmann (2013)
* assigning individual mass densities to all tracers and calculating interactions using the [atmostracers](https://github.com/OpenNWP/atmostracers) library

RFPET is meant to be applied at resolutions finer than 7 km. The default configuration is a 800 m resolution implementation nested into ICON-D2. Other nesting options are of course possible but require more work.

## Installation

### Dependencies

Everything is easy and quick to install. These instructions are for Ubuntu.

* cmake: `sudo apt-get install cmake`
* [atmostracers](https://github.com/OpenNWP/atmostracers)
* clone the RTE+RRTMGP repository: `git clone https://github.com/earth-system-radiation/rte-rrtmgp` and replace all `!$omp` with `!!$omp` in the subdirectories `rte`, `rrtmgp` as well as `extensions`

### Installation

```
git clone https://github.com/opennwp/rfpet.git
./setup_directories.sh
./compile.sh
```

### Execution

```
cd run_scripts
./ideal.sh
```

Output will be placed in the directory `output`.

## Fundamental literature

* Thuburn, John. (2008). Numerical wave propagation on the hexagonal C-grid. Journal of Computational Physics. 227. 5836-5858. 10.1016/j.jcp.2008.02.010. 
* Gassmann, Almut & Herzog, Hans-Joachim. (2008). Towards a consistent numerical compressible non‐hydrostatic model using generalized Hamiltonian tools. Quarterly Journal of the Royal Meteorological Society. 134. 1597 - 1613. 10.1002/qj.297.
* Thuburn, John et al. “Numerical representation of geostrophic modes on arbitrarily structured C-grids.” J. Comput. Phys. 228 (2009): 8321-8335.
* Ringler, Todd & Thuburn, John & Klemp, J. & Skamarock, W.C.. (2010). A unified approach to energy conservation and potential vorticity dynamics on arbitrarily structured C-grids. J. Comput. Physics. 229. 3065-3090. 10.1016/j.jcp.2009.12.007.
* Gassmann, A. (2013), A global hexagonal C‐grid non‐hydrostatic dynamical core (ICON‐IAP) designed for energetic consistency. Q.J.R. Meteorol. Soc., 139: 152-175. doi:10.1002/qj.1960






