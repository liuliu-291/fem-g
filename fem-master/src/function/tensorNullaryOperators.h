// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TENSORNULLARYOPERATORS
#define H_GMSHFEM_TENSORNULLARYOPERATORS

#include "Function.h"

// *************************************
// Nullary operators for tensor function
// *************************************

namespace gmshfem::function
{


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > &data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > &&data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > &data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > &&data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > *data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > bilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > *data);

  // /!\ This function may return tensor or tensor<4> function depending on the template parameter T_Rank
  template< class T_Scalar, unsigned int T_Rank = 2 >
  Function< T_Scalar, Degree(TensorRank< T_Rank >::value) > identity();

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > &data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > &&data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > &data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > &&data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > *data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > linearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > *data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > probeGradVectorView(const unsigned int view, const unsigned int step = 0, const int dim = -1, const double distanceMax = -1);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > probeTensorView(const unsigned int view, const unsigned int step = 0, const int dim = -1, const double distanceMax = -1);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< scalar::Precision< T_Scalar > > &z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > &data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< scalar::Precision< T_Scalar > > &&z, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > &&data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > &data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > &&data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< scalar::Precision< T_Scalar > > *z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > *data);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree2 > trilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > > *data);


} // namespace gmshfem::function

#endif // H_GMSHFEM_TENSORNULLARYOPERATORS
