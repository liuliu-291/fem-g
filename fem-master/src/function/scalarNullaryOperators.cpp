// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "scalarNullaryOperators.h"

#include "Domain.h"
#include "Message.h"
#include "NullaryNode.h"
#include "instantiate.h"
#include "nullaryOperations.h"

#include <gmsh.h>

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > bilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > &data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree0 > >(BilinearInterpolation< T_Scalar, Degree::Degree0 >(x, y, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , bilinearInterpolation, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > bilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > &&data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree0 > >(BilinearInterpolation< T_Scalar, Degree::Degree0 >(std::move(x), std::move(y), std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , bilinearInterpolation, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > &data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree0 > >(BilinearInterpolation< T_Scalar, Degree::Degree0 >(nx, ny, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , bilinearInterpolation, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > &&data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree0 > >(BilinearInterpolation< T_Scalar, Degree::Degree0 >(nx, ny, pMin, pMax, std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , bilinearInterpolation, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > bilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > *data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree0 > >(BilinearInterpolation< T_Scalar, Degree::Degree0 >(x, y, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , bilinearInterpolation, 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > *data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree0 > >(BilinearInterpolation< T_Scalar, Degree::Degree0 >(nx, ny, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , bilinearInterpolation, 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > linearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > &data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree0 > >(LinearInterpolation< T_Scalar, Degree::Degree0 >(x, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , linearInterpolation, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > linearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > &&data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree0 > >(LinearInterpolation< T_Scalar, Degree::Degree0 >(std::move(x), std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , linearInterpolation, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > &data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree0 > >(LinearInterpolation< T_Scalar, Degree::Degree0 >(nx, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , linearInterpolation, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > &&data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree0 > >(LinearInterpolation< T_Scalar, Degree::Degree0 >(nx, pMin, pMax, std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , linearInterpolation, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > linearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > *data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree0 > >(LinearInterpolation< T_Scalar, Degree::Degree0 >(x, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , linearInterpolation, 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > *data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree0 > >(LinearInterpolation< T_Scalar, Degree::Degree0 >(nx, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , linearInterpolation, 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > probeScalarView(const unsigned int view, const unsigned int step, const int dim, const double distanceMax)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< ProbeView< T_Scalar, Degree::Degree0 > >(ProbeView< T_Scalar, Degree::Degree0 >(view, step, dim, false, distanceMax)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , probeScalarView, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const int, const double))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > phi()
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Phi< T_Scalar > >());
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , phi, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS())


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > r2d()
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< R2D< T_Scalar > >());
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , r2d, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS())


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > r3d()
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< R3D< T_Scalar > >());
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , r3d, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS())


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > theta()
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Theta< T_Scalar > >());
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , theta, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS())
  
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > trilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< scalar::Precision< T_Scalar > > &z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > > &data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree0 > >(TrilinearInterpolation< T_Scalar, Degree::Degree0 >(x, y, z, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , trilinearInterpolation, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > trilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< scalar::Precision< T_Scalar > > &&z, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > > &&data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree0 > >(TrilinearInterpolation< T_Scalar, Degree::Degree0 >(std::move(x), std::move(y), std::move(z), std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , trilinearInterpolation, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > > &data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree0 > >(TrilinearInterpolation< T_Scalar, Degree::Degree0 >(nx, ny, nz, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , trilinearInterpolation, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > > &&data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree0 > >(TrilinearInterpolation< T_Scalar, Degree::Degree0 >(nx, ny, nz, pMin, pMax, std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , trilinearInterpolation, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const  unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > trilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< scalar::Precision< T_Scalar > > *z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > > *data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree0 > >(TrilinearInterpolation< T_Scalar, Degree::Degree0 >(x, y, z, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , trilinearInterpolation, 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > > *data)
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree0 > >(TrilinearInterpolation< T_Scalar, Degree::Degree0 >(nx, ny, nz, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , trilinearInterpolation, 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree0 >::Object > > > *))

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > x()
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< X< T_Scalar > >());
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , x, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS())


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > x(const domain::Domain &domain)
  {
    if(domain.isEmpty()) {
      msg::warning << "Try to get the x-coordinate of an empty domain: return 0" << msg::endl;
      return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Constant< T_Scalar, Degree::Degree0 > >(Constant< T_Scalar, Degree::Degree0 >(0.)));
    }
    std::vector< double > xMinAll, xMaxAll;
    double xMin = 0., yMin = 0., zMin = 0., xMax = 0., yMax = 0., zMax = 0.;
    for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
      gmsh::model::getBoundingBox(it->first, it->second, xMin, yMin, zMin, xMax, yMax, zMax);
      xMinAll.push_back(xMin);
      xMaxAll.push_back(xMax);
    }
    auto min = std::min_element(xMinAll.begin(), xMinAll.end());
    auto max = std::max_element(xMaxAll.begin(), xMaxAll.end());
    if(std::abs(*max - *min) <= scalar::Epsilon< T_Scalar >::value) {
      return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Constant< T_Scalar, Degree::Degree0 > >(Constant< T_Scalar, Degree::Degree0 >((*min + *max) / 2.)));
    }
    return x< T_Scalar >();
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , x, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const domain::Domain &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > y()
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Y< T_Scalar > >());
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , y, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS())


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > y(const domain::Domain &domain)
  {
    if(domain.isEmpty()) {
      msg::warning << "Try to get the y-coordinate of an empty domain: return 0" << msg::endl;
      return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Constant< T_Scalar, Degree::Degree0 > >(Constant< T_Scalar, Degree::Degree0 >(0.)));
    }
    std::vector< double > yMinAll, yMaxAll;
    double xMin = 0., yMin = 0., zMin = 0., xMax = 0., yMax = 0., zMax = 0.;
    for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
      gmsh::model::getBoundingBox(it->first, it->second, xMin, yMin, zMin, xMax, yMax, zMax);
      yMinAll.push_back(yMin);
      yMaxAll.push_back(yMax);
    }
    auto min = std::min_element(yMinAll.begin(), yMinAll.end());
    auto max = std::max_element(yMaxAll.begin(), yMaxAll.end());
    if(std::abs(*max - *min) <= scalar::Epsilon< T_Scalar >::value) {
      return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Constant< T_Scalar, Degree::Degree0 > >(Constant< T_Scalar, Degree::Degree0 >((*min + *max) / 2.)));
    }
    return y< T_Scalar >();
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , y, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const domain::Domain &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > z()
  {
    return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Z< T_Scalar > >());
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , z, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS())


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > z(const domain::Domain &domain)
  {
    if(domain.isEmpty()) {
      msg::warning << "Try to get the z-coordinate of an empty domain: return 0" << msg::endl;
      return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Constant< T_Scalar, Degree::Degree0 > >(Constant< T_Scalar, Degree::Degree0 >(0.)));
    }
    std::vector< double > zMinAll, zMaxAll;
    double xMin = 0., yMin = 0., zMin = 0., xMax = 0., yMax = 0., zMax = 0.;
    for(auto it = domain.cbegin(); it != domain.cend(); ++it) {
      gmsh::model::getBoundingBox(it->first, it->second, xMin, yMin, zMin, xMax, yMax, zMax);
      zMinAll.push_back(zMin);
      zMaxAll.push_back(zMax);
    }
    auto min = std::min_element(zMinAll.begin(), zMinAll.end());
    auto max = std::max_element(zMaxAll.begin(), zMaxAll.end());
    if(std::abs(*max - *min) <= scalar::Epsilon< T_Scalar >::value) {
      return Function< T_Scalar, Degree::Degree0 >(new NullaryNode< Constant< T_Scalar, Degree::Degree0 > >(Constant< T_Scalar, Degree::Degree0 >((*min + *max) / 2.)));
    }
    return z< T_Scalar >();
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree0 >), , z, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const domain::Domain &))


} // namespace gmshfem::function
