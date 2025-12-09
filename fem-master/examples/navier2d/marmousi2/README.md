# Elastodynamic wave propagation in ground (Marmousi2)
> Anthony Royer

Dependency(ies):
* GmshFEM : v 1.0.0
* Gmsh : v 4.9.0

## Problem description

Elastodynamic wave propagation in ground with Lysmer Kuhlemeyer absorbing boundary condition. 
The ground model is called 'Marmousi2' [1]

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
* `-lc [x]` where `[x]` is the mesh size (default = 20).
* `-gauss [x]` where `[x]` is the integration quadrature (default = "Gauss4").
* `-frequency [x]` where `[x]` is the wave frequency (default = 0.005 [Hz]).

## Reference

[1] https://wiki.seg.org/wiki/AGL_Elastic_Marmousi
