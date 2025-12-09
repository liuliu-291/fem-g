# Assessment of Local Absorbing Boundary Conditions for Heterogeneous Convected Helmholtz Problems

Dependency(ies):
* GmshFEM : v 1.0.0
* Gmsh : v 4.8.2

## Problem description

The routines solve a convected Helmholtz problem in a two-dimensional straight waveguide 
```math
\frac{D_0}{Dt}\left(\frac{1}{c_0^2}\frac{D_0}{Dt}\right) u - \rho_0^{-1}\nabla \cdot (\rho_0 \nabla ) u = 0, \quad \text{in} \; \Omega=[0, L]\times[0, H]
```
with either a uniform (`uniform/`) or non-uniform flow (`nonUniform/`)

Different ABCs are tested on the output boundary.
For the non-uniform flow, the relative $`L^2`$-error is computed thanks to a PML reference solution.

## Installation

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
* `-omega [x]` is the angular frequency for (`uniform/`) or `-kinf [x]` where `[x]` is the reference wavenumber for (`nonUniform/`).
* `-mode [x]` is the mode number for the input boundary condition.
* `-abcName [x]` is the ABC imposed on the outer boundary.
* `-padeOrder [x]` is the Pade order if `-abcName` is Pade-based.
* `-angle [x]` is the rotation angle of the branch-cut if `-abcName` is Pade-based.
* `-runcase [x]` selects a specific configuration for reproducibility.

Each subfolder contains
* the main script (`main.cpp`)
* a loop for testing the ABCs over a given frequency range (`main_omega.cpp`)
To run one of the loops just change the name in `CMakeLists.txt` and recompile

## Reference

> P. Marchner, X. Antoine, C. Geuzaine, H. Bériot,
Construction and Numerical Assessment of Local Absorbing Boundary Conditions for Heterogeneous Time-Harmonic Acoustic Problems, (2021), Preprint [HAL-03196015](https://hal.archives-ouvertes.fr/hal-03196015).

## Results reproductibility
##### Folder `uniform/`

Figure 6 can be obtained with `main_omega.cpp` and
```
  ./example -runcase [i]
```
with $`i=1`$ or $`i=2`$ for respectively Figures 6a and 6b.

##### Folder `nonUniform/`
Pressure maps from Figure 8 can be obtained by compiling `main.cpp` and running
```
  ./example -kinf 40
```
The solution `vpml.msh` can be open in Gmsh.
  
For Figure 9, compile `main_omega.cpp` and run the command
```
./example -runcase 1 -L [i]
```
with $`i=\{0.75,1,1.25,2\}`$ to obtain the 4 subfigures.

For Figure 10, compile `main_omega.cpp` and run the command
```
./example -runcase 2 -L [i]
```
with $`i=\{0.75,1,1.25,2\}`$ to obtain the 4 sub-figures.
