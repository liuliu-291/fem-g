
# Rectangular metal pipe waveguide

Dependency(ies):
* GmshFEM : v 1.0.0
* Gmsh : v 4.8.2

## Problem description

The routines solve a three-dimensional rectangular metal pipe waveguide problem 
```math
\textbf{rot}\textbf{rot}\textbf{E} - k^{2}\textbf{E} = 0, \quad \text{in} \; \Omega=[0, L]\times[0, a]\times[0, b]
```
High-order finite elements are used to discretize the weak formulations.
A given TM21 mode is imposed as input boundary condition at $`x=0`$. 
The physicality of electromagnetic field in the waveguide is enforced through the Silver-Muller absorbing boundary condition at  $`x=L`$.


## Installation

Clone the GmshFEM repository and run

```
  git clone https://gitlab.onelab.info/gmsh/fem.git
  mkdir build && cd build
  cmake ..
  make
  sudo make install
```

Once GmshFEM is installed, go to the `maxwell3d/` folder and run

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
* `-lc [x]` where `[x]` is the mesh size.
* `-k [x]` where `[x]` is the wavenumber in the waveguide.
* `-L [x]` is the waveguide length and ABC position.
* `-a [x]` is the waveguide height.
* `-b [x]` is the waveguide width.
* `-order [x]` is the order of the finite element hierarchical basis.
* `-useLagMult [x]` could be `true` or `false` for using or not Lagrange multipliers to impose the TM excitation 


The folder contains
* the main script (`main.cpp`)

