# Scattering of a plane wave on a plane surface
> Martin B.

Dependency(ies):
* GmshFEM : v 1.0.0
* Gmsh : v 4.8.2

## Project description

This code solves the scattering problem of the plane wave `u_inc` over a plane wall. 100% of the wave is reflected with a similar angle to the normal of the plane. With clever periodic BCs, the exact solution is recovered.

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
* `-k [x]` where `[x]` is a wavenumber.
* `-thetaInc [x]` where `[x]` is the angle of inscidence of the plan wave.
* `-L [x]` where `[x]` is the domain length



## References

//
