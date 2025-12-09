#include <gmshfem/Formulation.h>
#include <gmshfem/GmshFem.h>

using namespace gmshfem;
using namespace gmshfem::function;

TensorFunction< std::complex< double > > getInverseLorentzTensor(
    ScalarFunction< std::complex< double > > Mx, 
    ScalarFunction< std::complex< double > > My, 
    ScalarFunction< std::complex< double > > beta);

TensorFunction< std::complex< double > > getInversePMLJacobian(
    ScalarFunction< std::complex< double > > k, 
    std::pair<double, double> xlim, std::pair<double, double> ylim, double Lpml,
    ScalarFunction< std::complex< double > > beta,
    ScalarFunction< std::complex< double > > &detJ);