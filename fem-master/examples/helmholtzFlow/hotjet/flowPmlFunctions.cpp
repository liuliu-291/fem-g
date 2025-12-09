#include "flowPmlFunctions.h"

TensorFunction< std::complex< double > > getInverseLorentzTensor(
    ScalarFunction< std::complex< double > > Mx, ScalarFunction< std::complex< double > > My
    , ScalarFunction< std::complex< double > > beta)
{
  ScalarFunction< std::complex< double > > Alphax = 1 + Mx * Mx / (beta * (1 + beta));
  ScalarFunction< std::complex< double > > Alphay = 1 + My * My / (beta * (1 + beta));
  ScalarFunction< std::complex< double > > K = Mx * My / (beta * (1 + beta));
  return tensor< std::complex< double > >(beta * Alphay, -beta * K, 0., -beta * K, beta * Alphax, 0., 0., 0., 0.);
}

TensorFunction< std::complex< double > > getInversePMLJacobian(
    ScalarFunction< std::complex< double > > k, std::pair<double, double> xlim, std::pair<double, double> ylim, double Lpml,
    ScalarFunction< std::complex< double > > beta,
    ScalarFunction< std::complex< double > > &detJ)
{
  std::complex< double > im(0., 1.);
  ScalarFunction< std::complex< double > > Sigma0 = beta;
  // x-PML
  double Lxpml1 = xlim.first - Lpml; // PML left limit in x (-)
  double Lxpml2 = xlim.second + Lpml; // PML right limit in x (+)
  ScalarFunction< std::complex< double > > SigmaX = heaviside(-(x< std::complex< double > >() - xlim.first) ) * Sigma0 / ( abs(Lxpml1) - abs(x< std::complex< double > >()) )
    + heaviside( x< std::complex< double > >() - xlim.second ) * Sigma0 / ( abs(Lxpml2) - abs(x< std::complex< double > >()) );
  
  // y-PML
  double Lypml1 = ylim.first - Lpml; // PML lower limit in y (-)
  double Lypml2 = ylim.second + Lpml; // PML upper limit in y (+) 
  ScalarFunction< std::complex< double > > SigmaY = heaviside(-(y< std::complex< double > >() - ylim.first) ) * Sigma0 / ( abs(Lypml1) - abs(y< std::complex< double > >()) )
    + heaviside( y< std::complex< double > >() - ylim.second ) * Sigma0 / ( abs(Lypml2) - abs(y< std::complex< double > >()) );
  // PML stretchings
  ScalarFunction< std::complex< double > > gammaX = 1 - (im / k) * SigmaX;
  ScalarFunction< std::complex< double > > gammaY = 1 - (im / k) * SigmaY;

  detJ = gammaX * gammaY;
  return tensor< std::complex< double > >(1. / gammaX, 0., 0., 0., 1. / gammaY, 0., 0., 0., 0.);
}