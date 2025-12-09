// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MATRIXCRSFAST
#define H_GMSHFEM_MATRIXCRSFAST

#include "Matrix.h"

namespace gmshfem::algebra
{


  //  Matrix format used during the matrix conversions and to efficiently applied const functions such as save, saveSpyPlot

  template< class T_Scalar >
  class MatrixCRSFast final : public Matrix< T_Scalar >
  {
   protected:
    const unsigned long long *_ai;
    const unsigned long long *_aj;
    const T_Scalar *_a;

    virtual void _buildPetsc() override;

   public:
    MatrixCRSFast();
    MatrixCRSFast(const unsigned long long size0, const unsigned long long size1, const std::vector< unsigned long long > &ai, const std::vector< unsigned long long > &aj, const std::vector< T_Scalar > &a);
    MatrixCRSFast(const MatrixCRSFast &other);
    MatrixCRSFast(MatrixCRSFast &&other);
    virtual ~MatrixCRSFast();

    const Eigen::SparseMatrix< T_Scalar, Eigen::RowMajor, long long > getEigen() const;

    MatrixCRSFast &operator=(const MatrixCRSFast &other);
    MatrixCRSFast &operator=(MatrixCRSFast &&other);

    MatrixCRSFast &operator=(const Matrix< T_Scalar > &other) override;

    virtual void extract(const system::MatrixFactory< T_Scalar > *matrix) override;
    virtual void clear() override;

    MatrixFormat format() const override;

    virtual unsigned long long numberOfNonZero() const override;

    virtual bool isSymmetric(const scalar::Precision< T_Scalar > tolerance = 10. * scalar::Epsilon< T_Scalar >::value) const override;
    virtual bool isHermitian(const scalar::Precision< T_Scalar > tolerance = 10. * scalar::Epsilon< T_Scalar >::value) const override;

    virtual void save(const std::string &path) const override;
    virtual void saveSpyPlot(const std::string &path, const unsigned int pointSize = 1, const common::Color &zero = common::Color::white(), const common::Color &noZero = common::Color::black()) const override;

    const unsigned long long *ai() const;
    const unsigned long long *aj() const;
    const T_Scalar *a() const;
  };


} // namespace gmshfem::algebra

#endif // H_GMSHFEM_MATRIXCRSFAST
