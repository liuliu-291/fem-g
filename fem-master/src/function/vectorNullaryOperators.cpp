// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "vectorNullaryOperators.h"

#include "NullaryNode.h"
#include "instantiate.h"
#include "nullaryOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > bilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > &data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree1 > >(BilinearInterpolation< T_Scalar, Degree::Degree1 >(x, y, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , bilinearInterpolation, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > bilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > &&data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree1 > >(BilinearInterpolation< T_Scalar, Degree::Degree1 >(std::move(x), std::move(y), std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , bilinearInterpolation, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > &data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree1 > >(BilinearInterpolation< T_Scalar, Degree::Degree1 >(nx, ny, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , bilinearInterpolation, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > &&data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree1 > >(BilinearInterpolation< T_Scalar, Degree::Degree1 >(nx, ny, pMin, pMax, std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , bilinearInterpolation, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > bilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > *data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree1 > >(BilinearInterpolation< T_Scalar, Degree::Degree1 >(x, y, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , bilinearInterpolation, 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > *data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree1 > >(BilinearInterpolation< T_Scalar, Degree::Degree1 >(nx, ny, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , bilinearInterpolation, 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > linearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > &data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree1 > >(LinearInterpolation< T_Scalar, Degree::Degree1 >(x, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , linearInterpolation, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > linearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > &&data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree1 > >(LinearInterpolation< T_Scalar, Degree::Degree1 >(std::move(x), std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , linearInterpolation, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > &data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree1 > >(LinearInterpolation< T_Scalar, Degree::Degree1 >(nx, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , linearInterpolation, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > &&data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree1 > >(LinearInterpolation< T_Scalar, Degree::Degree1 >(nx, pMin, pMax, std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , linearInterpolation, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > linearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > *data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree1 > >(LinearInterpolation< T_Scalar, Degree::Degree1 >(x, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , linearInterpolation, 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > *data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree1 > >(LinearInterpolation< T_Scalar, Degree::Degree1 >(nx, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , linearInterpolation, 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > probeGradScalarView(const unsigned int view, const unsigned int step, const int dim, const double distanceMax)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< ProbeView< T_Scalar, Degree::Degree1 > >(ProbeView< T_Scalar, Degree::Degree1 >(view, step, dim, true, distanceMax)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , probeGradScalarView, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const int, const double))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > probeVectorView(const unsigned int view, const unsigned int step, const int dim, const double distanceMax)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< ProbeView< T_Scalar, Degree::Degree1 > >(ProbeView< T_Scalar, Degree::Degree1 >(view, step, dim, false, distanceMax)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , probeVectorView, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const int, const double))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > normal()
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< Normal< T_Scalar > >());
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , normal, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS())


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > tangent(const unsigned int component)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< Tangent< T_Scalar > >(Tangent< T_Scalar >(component)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , tangent, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int))
  
  
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > trilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< scalar::Precision< T_Scalar > > &z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > > &data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree1 > >(TrilinearInterpolation< T_Scalar, Degree::Degree1 >(x, y, z, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , trilinearInterpolation, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > trilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< scalar::Precision< T_Scalar > > &&z, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > > &&data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree1 > >(TrilinearInterpolation< T_Scalar, Degree::Degree1 >(std::move(x), std::move(y), std::move(z), std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , trilinearInterpolation, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > > &data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree1 > >(TrilinearInterpolation< T_Scalar, Degree::Degree1 >(nx, ny, nz, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , trilinearInterpolation, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > > &&data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree1 > >(TrilinearInterpolation< T_Scalar, Degree::Degree1 >(nx, ny, nz, pMin, pMax, std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , trilinearInterpolation, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const  unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > trilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< scalar::Precision< T_Scalar > > *z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > > *data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree1 > >(TrilinearInterpolation< T_Scalar, Degree::Degree1 >(x, y, z, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , trilinearInterpolation, 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree1 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > > *data)
  {
    return Function< T_Scalar, Degree::Degree1 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree1 > >(TrilinearInterpolation< T_Scalar, Degree::Degree1 >(nx, ny, nz, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree1 >), , trilinearInterpolation, 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree1 >::Object > > > *))


} // namespace gmshfem::function
