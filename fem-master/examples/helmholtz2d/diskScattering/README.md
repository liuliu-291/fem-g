# Scattering of a plane wave by a disk
> Royer A.

Dependency(ies):
* GmshFEM : v 1.0.0
* Gmsh : v 4.8.2

## Project description

This code solves the scattering problem of the plane wave `u_inc` by a sound-hard or a sound-soft circular cylinder of radius R centered at the origin and generating a scatted unknown field `u`. 

## Installation


```
  mkdir build && cd build
  cmake -DCMAKE_CXX_FLAGS="-g -march=native" ..
  make
```

## Usage

Simply run:
```
  ./example [PARAM]
```
with `[PARAM]`:
* `-lc [x]` where `[x]` is the mesh size.
* `-order [x]` where `[x]` is the mesh order.
* `-problem [x]` where `[x]` could be `hard` or `soft` and corresponds to a hard or a soft scattering problem.
* `-gauss [x]` where `[x]` means a Gauss quadrature suited for integrating xth order polynomials.
* `-k [x]` where `[x]` is a wavenumber.
* `-theta [x]` where `[x]` is the angle of inscidence of the plan wave.
* `-abcName [x]` where `[x]` is the ABC imposed on the exterior boundary of the domain. It could be `sommerfeld`, `bayliss-turkel`, `pade`or `pml`.
* `-geometry [x]` where `[x]`could be 0 for a cylindrical domain or 1 for a square domain.
* `-saveUExact` to save the exact solution `u_ex`.
* `-saveError` to save the error field `u - u_ex`
* `-computeError` to compute the L_2 error `sqrt(int_\Omega |u - u_ex|^2 / int_\Omega |u_ex|^2)`.
* `-padeOrder [x]` where `[x]` is the Pade order if `abcName = pade`. 
* `-withCorner` to impose the corner traitment if `abcName = pade` and `geometry = 1`.


## References

//
