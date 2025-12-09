# This directory contains a basis function visualition tool

## Building

* `mkdir build`
* `cd build`
* `cmake ..`
* `make`

## Running

Simply run:
```
  ./demo -form {form} -bf {bf} [-derivative {derivative}] [-element {element}] [-order {order}] [-lc {lc}] [-orientation {orientation}]
```
with the parameters:
* `{form}` the differential form ("form0", "form1", "form2" or "form3");
* `{bf}` the basis function type depending on the chosen differential form.

and the optional parameters:
* `{derivative}` (default `0`) returns the external derivative of the basis function (value "1") or not (value "0");
* `{element}` (default `line`) selects the element type name ("point", "line", "triangle", "quadrangle", "tetrahedron", "hexahedron", "prism", "pyramid");
* `{order}`(default `1`) selects the basis function order (this is meaninless for iso-parametric basis functions);
* `{lc}` (default `0.1`) specifies the size of the mesh on which basis functions are saved;
* `{orientation}` (default `0`) specifies the orientation of the basis functions.


