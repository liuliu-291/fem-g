# Absorbing boundary condtions for free field convected propagation

Dependency(ies):
* GmshFEM : v 1.0.0
* Gmsh

## Problem description

These routines solve a convected Helmholtz propagation problem with a Dirac right-hand-side, either in a circular or square domain. 
Different ABCs are tested on the output boundary, and are compared to the analytical solution.

## Installation and usage 
Simply run 

```
  mkdir build && cd build
  cmake ..
  make
  ./example [PARAM]
```
with `[PARAM]`:
* `-geometry [x]` where `[x]` is the domain geometry (0: disk, 1: square)
* `-k [x]` is the freefield wavenumber, $`k>0`$
* `-R [x]` where `[x]` is the disk radius or square size, $`R>0`$
* `-xs [x]` and `-ys [x]` specify the point source location
* `-M [x]` is the mean flow Mach number, $`0<M<1`$
* `-theta [x]` is the mean flow angle orientation, $\theta \in [0, 2\pi]$
* `-abcName [x]` is the ABC type. The available choices are `ABC-0`, `ABC-2`, `Pade` and `pml`
* `-padeOrder [x]` and `-angle [x]` specify respectively the number of auxiliary function and branch rotation angle for the `Pade` abc.

By default a circular domain with a Padé condition is used, for a centered source at the frequency $k=6\pi$, with a mean flow at $M=0.8$ oriented at $\theta=\pi/4$.

## Reference
> P. Marchner, H. Bériot, S. Le Bras, X. Antoine and C. Geuzaine. A domain decomposition solver for large scale time-harmonic flow acoustics problems. (2023), Preprint [HAL-04254633](https://hal.science/hal-04254633).

## Results reproducibility
We assess the performance of different ABCs in the case of a radiating point source set in a circular domain.
A loop on the Mach number is launched by executing the file `runTests_ABC_circle.sh`, which allows to obtain the data shown at the end of Section 2.

