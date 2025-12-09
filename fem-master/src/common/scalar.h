// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SCALAR
#define H_GMSHFEM_SCALAR

#include "Exception.h"
#include "OmpInterface.h"

namespace Eigen
{
  class eigen_assert_exception : public gmshfem::common::Exception
  {
  };
} // namespace Eigen

#include <Eigen/Eigen>
#include <Eigen/StdVector>
#include <complex>
#include <limits>
#include <string>
#include <vector>

namespace gmshfem::scalar
{


  // scalarPrecision

  template< class T_Scalar >
  struct scalarPrecision {
    typedef T_Scalar T_Precision;
  };

  template< class T_ComplexScalar >
  struct scalarPrecision< std::complex< T_ComplexScalar > > {
    typedef T_ComplexScalar T_Precision;
  };

  template< class T_Scalar >
  using Precision = typename scalarPrecision< T_Scalar >::T_Precision;

  // scalarComplexPrecision

  template< class T_Scalar >
  struct scalarComplexPrecision {
    typedef std::complex< T_Scalar > T_Precision;
  };

  template< class T_ComplexScalar >
  struct scalarComplexPrecision< std::complex< T_ComplexScalar > > {
    typedef std::complex< T_ComplexScalar > T_Precision;
  };

  template< class T_Scalar >
  using ComplexPrecision = typename scalarComplexPrecision< T_Scalar >::T_Precision;

  // scalarPrecisionDigits

  template< class T_PScalar >
  struct scalarPrecisionDigits {
    constexpr static int value = 15;
  };

  template<>
  struct scalarPrecisionDigits< double > {
    constexpr static int value = 15;
  };

  template<>
  struct scalarPrecisionDigits< float > {
    constexpr static int value = 7;
  };

  template< class T_Scalar >
  using PrecisionDigits = scalarPrecisionDigits< typename scalarPrecision< T_Scalar >::T_Precision >;

  // scalarEpsilon

  template< class T_PScalar >
  struct scalarEpsilon {
    constexpr static T_PScalar value = std::numeric_limits< T_PScalar >::epsilon();
  };

  template< class T_Scalar >
  using Epsilon = scalarEpsilon< typename scalarPrecision< T_Scalar >::T_Precision >;

  // scalarIsComplex

  template< class T_Scalar >
  struct scalarIsComplex {
    constexpr static bool value = false;
  };

  template< class T_ComplexScalar >
  struct scalarIsComplex< std::complex< T_ComplexScalar > > {
    constexpr static bool value = true;
  };

  template< class T_Scalar >
  using IsComplex = scalarIsComplex< T_Scalar >;

  // scalarIsReal

  template< class T_Scalar >
  struct scalarIsReal {
    constexpr static bool value = true;
  };

  template< class T_ComplexScalar >
  struct scalarIsReal< std::complex< T_ComplexScalar > > {
    constexpr static bool value = false;
  };

  template< class T_Scalar >
  using IsReal = scalarIsReal< T_Scalar >;

  // scalarName

  template< class T_Scalar >
  struct scalarName {
    constexpr static const char *name = "Unknown scalar type";
  };

  template<>
  struct scalarName< std::complex< double > > {
    constexpr static const char *name = "std::complex< double >";
  };

  template<>
  struct scalarName< double > {
    constexpr static const char *name = "double";
  };

  template<>
  struct scalarName< std::complex< float > > {
    constexpr static const char *name = "std::complex< float >";
  };

  template<>
  struct scalarName< float > {
    constexpr static const char *name = "float";
  };

  template< class T_Scalar >
  using Name = scalarName< T_Scalar >;

  // move

  template< class T_Scalar1, class T_Scalar2 >
  void move(std::vector< T_Scalar1 > &v1, std::vector< T_Scalar2 > &v2)
  {
    v1.resize(v2.size());
    for(auto i = 0ULL; i < v2.size(); ++i) {
      v1[i] = v2[i];
    }
    v2.clear();
    v2.shrink_to_fit();
  }

  template< class T_Scalar1 >
  void move(std::vector< T_Scalar1 > &v1, std::vector< T_Scalar1 > &v2)
  {
    v1 = std::move(v2);
  }

  // copy

  template< class T_Scalar1, class T_Scalar2 >
  void copy(std::vector< T_Scalar1 > &v1, const std::vector< T_Scalar2 > &v2)
  {
    v1.resize(v2.size());
    for(auto i = 0ULL; i < v2.size(); ++i) {
      v1[i] = static_cast< T_Scalar1 >(v2[i]);
    }
  }

  template< class T_Scalar1 >
  void copy(std::vector< T_Scalar1 > &v1, const std::vector< T_Scalar1 > &v2)
  {
    v1.resize(v2.size());
    v1 = v2;
  }


} // namespace gmshfem::scalar

namespace Eigen
{
  template< class T_Scalar >
  using VectorX = Eigen::Matrix< T_Scalar, Eigen::Dynamic, 1 >;

  template< class T_Scalar >
  using MatrixX = Eigen::Matrix< T_Scalar, Eigen::Dynamic, Eigen::Dynamic >;

  template< class T_Scalar >
  using Vector3 = Eigen::Matrix< T_Scalar, 3, 1 >;

  template< class T_Scalar >
  using Matrix3 = Eigen::Matrix< T_Scalar, 3, 3 >;

  template< class T_Scalar >
  using Tensor3 = Eigen::Matrix< Eigen::Matrix< T_Scalar, 3, 3 >, 3, 3 >;
} // namespace Eigen

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector3< std::complex< double > >)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector3< double >)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector3< std::complex< float > >)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Vector3< float >)

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3< std::complex< double > >)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3< double >)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3< std::complex< float > >)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3< float >)

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Tensor3< std::complex< double > >)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Tensor3< double >)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Tensor3< std::complex< float > >)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Tensor3< float >)


#endif // H_GMSHFEM_SCALAR
