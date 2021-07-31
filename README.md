## Regional Forecasting with Poisson brackets in Exner-Theta formulation.

RFPET is a numerical weather prediction model. Properties:

* regional
* non-hydrostatic
* Eulerian
* quadrilateral mesh on a (rotated) latitude-longitude grid
* using finite volume methods instead of grid point formulations
* time stepping: two-time-level Runge-Kutta scheme, HEVI (horizontally explicit, vertically implicit)
* radiation: coupled to RTE+RRTMGP
* uses the Poisson bracket formulation by Gassmann and Herzog (2008) and Gassmann (2013)

## Fundamental literature

* Thuburn, John. (2008). Numerical wave propagation on the hexagonal C-grid. Journal of Computational Physics. 227. 5836-5858. 10.1016/j.jcp.2008.02.010. 
* Gassmann, Almut & Herzog, Hans-Joachim. (2008). Towards a consistent numerical compressible non‐hydrostatic model using generalized Hamiltonian tools. Quarterly Journal of the Royal Meteorological Society. 134. 1597 - 1613. 10.1002/qj.297.
* Thuburn, John et al. “Numerical representation of geostrophic modes on arbitrarily structured C-grids.” J. Comput. Phys. 228 (2009): 8321-8335.
* Ringler, Todd & Thuburn, John & Klemp, J. & Skamarock, W.C.. (2010). A unified approach to energy conservation and potential vorticity dynamics on arbitrarily structured C-grids. J. Comput. Physics. 229. 3065-3090. 10.1016/j.jcp.2009.12.007.
* Gassmann, A. (2013), A global hexagonal C‐grid non‐hydrostatic dynamical core (ICON‐IAP) designed for energetic consistency. Q.J.R. Meteorol. Soc., 139: 152-175. doi:10.1002/qj.1960

## Installation

```
git clone https://github.com/opennwp/rfpet.git
./setup_directories.sh
```
