#!/bin/zsh

demos=("../demos/shell/polarShell"
       "../demos/poisson/coupling"
       "../demos/poisson/dirichlet"
       "../demos/poisson/functions"
       "../demos/poisson/multiRHS"
       "../demos/poisson/neumann"
       "../demos/helmholtz"
       "../demos/navier"
       "../demos/hierarchical/helmholtz"
       )

path=`pwd`
set -e

for i in "${demos[@]}"
do
  echo "### Compiling $i ..."
  cd "$i"
  mkdir build
  cd build
  cmake ..
  make -j16
  echo "### Done"
  echo "### Running $i ..."
  ./demo -maxThreads 16
  echo "### Done"
  cd ..
  rm -r build
  cd "$path"
done

echo "### Everything works fine!"
