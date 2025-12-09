#include <gmsh.h>
#include <math.h>
#include <iostream>

void meshJet(const double sigma, const double SourceMeshRefinement, const int NrOfLayers, const int MeshElemOrder, 
    std::pair<double, double> &xlim, std::pair<double, double> &ylim, double &lc);