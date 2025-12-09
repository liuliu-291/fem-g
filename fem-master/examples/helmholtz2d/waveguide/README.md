# Assessment of Local Absorbing Boundary Conditions for Heterogeneous Helmholtz Problems

Dependency(ies):
* GmshFEM : v 1.0.0
* Gmsh : v 4.8.2

## Problem description

The routines solve a two-dimensional straight waveguide Helmholtz problem
```math
\rho_0^{-1}\partial_x(\rho_0\partial_x) u + \rho_0^{-1}\partial_y(\rho_0\partial_y) u + \omega^2 c_0^{-2} u = 0, \quad \text{in} \; \Omega=[0, L]\times[0, H]
```
with a either
* a longitudinal, linear speed of sound profile $`c_0(x)=ax+b, \, a>0`$ (folder `longitudinal/`),
* a transverse, Gaussian $`c_0(y) ,\, \rho_0(y)`$ variation of the speed of sound and density (folder `transverse/`). In addition, a runcase allows to test a constant piecewise speed of sound profile.

High-order finite elements are used to discretize the weak formulations.
A given mode of the form $`g=\cos(n \pi y / H), n\in\mathbb{N}`$ is imposed as input boundary condition at $`x=0`$.  
Homogeneous Dirichlet or Neumann boundary conditions are available for the upper and lower walls.

Different ABCs are tested on the output boundary.
* In the folder `longitudinal/`, the relative $`L^2`$-error is computed thanks to an analytic solution.
* In the folder `transverse/`, the relative $`L^2`$-error is computed thanks to a PML reference solution.
For the PML reference, Bermúdez et al.’s function is used with the parameter $`\sigma_0=40`$ and $`N_{pml}=40`$ layers in order to ensure the attenuation of all modes including grazing modes. The reference computational domain may also be extended before the application of the PML thanks to the variable $`Lext`$.

## Installation

For the examples in the folder `longitudinal`, GmshFEM must be compiled with the Complex Bessel library.
Clone the GmshFEM repository and run

```
  git clone https://gitlab.onelab.info/gmsh/fem.git
  mkdir build && cd build
  cmake -DENABLE_BESSEL=1 ..
  make
  sudo make install
```

Once GmshFEM is installed, go to the `longitudinal/` or `transverse/` folder and run

```
  mkdir build && cd build
  cmake ..
  make
```

## Usage

Simply run:
```
  ./example [PARAM]
```
with `[PARAM]`:
* `-h [x]` where `[x]` is the mesh size.
* `-L [x]` is the waveguide length and ABC position.
* `-H [x]` is the waveguide height.
* `-FEMorder [x]` is the order of the finite element hierarchical basis.
* `-problem [x]` could be `hard` or `soft` and corresponds to a homogeneous Neumann or Dirichlet boundary conditions for the upper and lower duct walls.
* `-omega [x]` is the angular frequency.
* `-mode [x]` is the mode number for the input boundary condition.
* `-abcName [x]` is the ABC imposed on the outer boundary.
* `-padeOrder [x]` is the Pade order if `-abcName` is Pade-based.
* `-angle [x]` is the rotation angle of the branch-cut if `-abcName` is Pade-based.
* `-runcase [x]` selects a specific configuration for reproducibility.

Each subfolder contains
* the main script (`main.cpp`)
* a loop for testing the ABCs over a given frequency range (`main_omega.cpp`)
* a loop for varying the position of the ABC (`main_length.cpp`)
To run one of the loops just change the script name in `CMakeLists.txt` and recompile

## Reference

> P. Marchner, X. Antoine, C. Geuzaine, H. Bériot,
Construction and Numerical Assessment of Local Absorbing Boundary Conditions for Heterogeneous Time-Harmonic Acoustic Problems, (2021), Preprint [HAL-03196015](https://hal.archives-ouvertes.fr/hal-03196015).

## Results reproductibility
##### Section 3, folder `longitudinal/`
Pressure maps from Figure 2 can be obtained by compiling `main.cpp` and running for example
```
  ./example -omega 20
```
The solution `v_exact.msh` can be open in `Gmsh`.

For Figure 3 use `main_omega.cpp` and run
```
  ./example -runcase [i]
```
with $`i=1`$ for Taylor-based ABCs or $`i=2`$ for Padé-based ABCs.

For Figure 4 use `main_length.cpp`
```
  ./example -runcase [i]
```
with $`i=\{1,2,3\}`$ to obtain the data from Figures 4a, 4b and 4c, respectively.

##### Section 5, folder `transverse/`
For Figures 13a, 13b and 13c, use respectively `main.cpp`, `main_omega.cpp` and `main_length.cpp` with
```
./example -runcase 1
```

Figure 15a and 15b, use respectively `main.cpp` and `main_omega.cpp` with
```
./example -runcase 2
```
