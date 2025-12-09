// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "tensorNullaryOperators.h"

#include "NullaryNode.h"
#include "instantiate.h"
#include "nullaryOperations.h"

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > &data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree2 > >(BilinearInterpolation< T_Scalar, Degree::Degree2 >(x, y, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , bilinearInterpolation, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > &&data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree2 > >(BilinearInterpolation< T_Scalar, Degree::Degree2 >(std::move(x), std::move(y), std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , bilinearInterpolation, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > &data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree2 > >(BilinearInterpolation< T_Scalar, Degree::Degree2 >(nx, ny, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , bilinearInterpolation, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > &&data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree2 > >(BilinearInterpolation< T_Scalar, Degree::Degree2 >(nx, ny, pMin, pMax, std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , bilinearInterpolation, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > *data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree2 > >(BilinearInterpolation< T_Scalar, Degree::Degree2 >(x, y, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , bilinearInterpolation, 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > *data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< BilinearInterpolation< T_Scalar, Degree::Degree2 > >(BilinearInterpolation< T_Scalar, Degree::Degree2 >(nx, ny, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , bilinearInterpolation, 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > *))


  template< class T_Scalar, unsigned int T_Rank >
  Function< T_Scalar, Degree(TensorRank< T_Rank >::value) > identity()
  {
    if constexpr (T_Rank == 2) {
      return Function< T_Scalar,Degree(TensorRank< T_Rank >::value) >(new function::NullaryNode< function::Constant< T_Scalar, Degree::Degree2 > >(function::Constant< T_Scalar, Degree::Degree2 >(Eigen::Matrix3< T_Scalar >::Identity())));
    }
    else if constexpr (T_Rank == 4) {
      typename gmshfem::MathObject< T_Scalar, gmshfem::Degree::Degree4 >::Object tmp;
      for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
          for(int k = 0; k < 3; ++k) {
            for(int l = 0; l < 3; ++l) {
              if(i == k && j == l) {
                tmp(i,j)(k,l) = 1.;
              }
              else {
                tmp(i,j)(k,l) = 0.;
              }
            }
          }
        }
      }
      return Function< T_Scalar, Degree::Degree4 >(new function::NullaryNode< function::Constant< T_Scalar, Degree::Degree4 > >(function::Constant< T_Scalar, Degree::Degree4 >(tmp)));
    }
  }

  INSTANTIATE_FCT_2(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree(TensorRank< TEMPLATE_PARAM_2 >::value) >), , identity, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), 2, unsigned int, TEMPLATE_ARGS(2, 4), TEMPLATE_PARAMS())


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > &data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree2 > >(LinearInterpolation< T_Scalar, Degree::Degree2 >(x, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , linearInterpolation, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > &&data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree2 > >(LinearInterpolation< T_Scalar, Degree::Degree2 >(std::move(x), std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , linearInterpolation, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > &data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree2 > >(LinearInterpolation< T_Scalar, Degree::Degree2 >(nx, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , linearInterpolation, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > &&data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree2 > >(LinearInterpolation< T_Scalar, Degree::Degree2 >(nx, pMin, pMax, std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , linearInterpolation, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > *data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree2 > >(LinearInterpolation< T_Scalar, Degree::Degree2 >(x, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , linearInterpolation, 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > *data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< LinearInterpolation< T_Scalar, Degree::Degree2 > >(LinearInterpolation< T_Scalar, Degree::Degree2 >(nx, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , linearInterpolation, 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const scalar::Precision< TEMPLATE_PARAM_1 > &, const scalar::Precision< TEMPLATE_PARAM_1 > &, const std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > probeGradVectorView(const unsigned int view, const unsigned int step, const int dim, const double distanceMax)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< ProbeView< T_Scalar, Degree::Degree2 > >(ProbeView< T_Scalar, Degree::Degree2 >(view, step, dim, true, distanceMax)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , probeGradVectorView, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const int, const double))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > probeTensorView(const unsigned int view, const unsigned int step, const int dim, const double distanceMax)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< ProbeView< T_Scalar, Degree::Degree2 > >(ProbeView< T_Scalar, Degree::Degree2 >(view, step, dim, false, distanceMax)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , probeTensorView, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const int, const double))
  
  
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< scalar::Precision< T_Scalar > > &z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > &data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree2 > >(TrilinearInterpolation< T_Scalar, Degree::Degree2 >(x, y, z, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , trilinearInterpolation, 0, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< scalar::Precision< T_Scalar > > &&z, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > &&data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree2 > >(TrilinearInterpolation< T_Scalar, Degree::Degree2 >(std::move(x), std::move(y), std::move(z), std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , trilinearInterpolation, 1, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &&, std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > &data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree2 > >(TrilinearInterpolation< T_Scalar, Degree::Degree2 >(nx, ny, nz, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , trilinearInterpolation, 2, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > > &))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > &&data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree2 > >(TrilinearInterpolation< T_Scalar, Degree::Degree2 >(nx, ny, nz, pMin, pMax, std::move(data))));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , trilinearInterpolation, 3, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const  unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > > &&))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< scalar::Precision< T_Scalar > > *z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > *data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree2 > >(TrilinearInterpolation< T_Scalar, Degree::Degree2 >(x, y, z, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , trilinearInterpolation, 4, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > *, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > > *))


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > *data)
  {
    return Function< T_Scalar, Degree::Degree2 >(new NullaryNode< TrilinearInterpolation< T_Scalar, Degree::Degree2 > >(TrilinearInterpolation< T_Scalar, Degree::Degree2 >(nx, ny, nz, pMin, pMax, data)));
  }

  INSTANTIATE_FCT(TEMPLATE_RETURN(Function< TEMPLATE_PARAM_1, Degree::Degree2 >), , trilinearInterpolation, 5, 4, class, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_PARAMS(const unsigned int, const unsigned int, const unsigned int, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< scalar::Precision< TEMPLATE_PARAM_1 > > &, const std::vector< std::vector< std::vector< typename MathObject< TEMPLATE_PARAM_1, Degree::Degree2 >::Object > > > *))


} // namespace gmshfem::function
