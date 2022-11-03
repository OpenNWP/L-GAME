# Limited-area GAME (L-GAME)

L-GAME is a numerical weather prediction (NWP) model. It is the application of the theory behind [GAME](https://github.com/opennwp/game) to a regional quadrilateral grid. Properties:

* non-hydrostatic
* Eulerian
* (rotated) latitude-longitude grid
* C-grid
* uses a hybrid of finite volume and finite difference methods
* time stepping: two-time-level predictor-corrector scheme, modified into a HEVI (horizontally explicit, vertically implicit) and forward-backward scheme for stability, horizontal pressure gradient extrapolated and kept constant
* radiation: coupled to [RTE+RRTMGP](https://github.com/earth-system-radiation/rte-rrtmgp)
* uses the Poisson bracket formulation by Gassmann and Herzog (2008) and Gassmann (2013)
* assigns individual mass densities to all tracers and calculates interactions

## Installation

### Dependencies

Everything is easy and quick to install.

	sudo apt-get install gfortran make cmake wget python3-pip libnetcdff-dev

* Clone the RTE+RRTMGP repository: `git clone https://github.com/earth-system-radiation/rte-rrtmgp`
* `pip3 install global-land-mask`

#### For using the plotting routines

The following packages are additionally required if you want to make use of the plotting routines:

* Python visualization library scitools-iris (installation manual: https://scitools-iris.readthedocs.io/en/latest/installing.html#installing-from-source-without-conda-on-debian-based-linux-distros-developers)

### Download and compilation

```
git clone https://github.com/OpenNWP/L-GAME.git
cd L-GAME
./create_directories.sh
./compile.sh
```

## Execution

Modify the variable lgame_home_dir in the run scripts (files in the directory run_scripts). Then you can use these files to execute certain model runs, for example:

```
./run_scripts/schaer.sh
```

Output will be placed in the directory `output`.

## Fundamental literature

* Thuburn, John. (2008). Numerical wave propagation on the hexagonal C-grid. Journal of Computational Physics. 227. 5836-5858. 10.1016/j.jcp.2008.02.010. 
* Gassmann, Almut & Herzog, Hans-Joachim. (2008). Towards a consistent numerical compressible non‐hydrostatic model using generalized Hamiltonian tools. Quarterly Journal of the Royal Meteorological Society. 134. 1597 - 1613. 10.1002/qj.297.
* Thuburn, John et al. “Numerical representation of geostrophic modes on arbitrarily structured C-grids.” J. Comput. Phys. 228 (2009): 8321-8335.
* Ringler, Todd & Thuburn, John & Klemp, J. & Skamarock, W.C.. (2010). A unified approach to energy conservation and potential vorticity dynamics on arbitrarily structured C-grids. J. Comput. Physics. 229. 3065-3090. 10.1016/j.jcp.2009.12.007.
* Gassmann, A. (2013), A global hexagonal C‐grid non‐hydrostatic dynamical core (ICON‐IAP) designed for energetic consistency. Q.J.R. Meteorol. Soc., 139: 152-175. doi:10.1002/qj.1960






