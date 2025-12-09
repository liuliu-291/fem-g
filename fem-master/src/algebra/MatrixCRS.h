// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MATRIXCRS
#define H_GMSHFEM_MATRIXCRS

#include "Matrix.h"
#include "MatrixCRSFast.h"

namespace gmshfem::system
{
  template< class T_Scalar >
  class MatrixFactory;
};

namespace gmshfem::algebra
{


  template< class T_Scalar >
  class MatrixCRS final : public Matrix< T_Scalar >
  {
   protected:
    std::vector< unsigned long long > _ai;
    std::vector< unsigned long long > _aj;
    std::vector< T_Scalar > _a;

    virtual void _buildPetsc() override;

   public:
    MatrixCRS();
    MatrixCRS(const MatrixCRS &other);
    MatrixCRS(MatrixCRS &&other);
    MatrixCRS(std::vector< unsigned long long > &&ai, std::vector< unsigned long long > &&aj, std::vector< T_Scalar > &&a);
    virtual ~MatrixCRS();

    Eigen::SparseMatrix< T_Scalar, Eigen::RowMajor, long long > getEigen();

    MatrixCRS &operator=(const MatrixCRS &other);
    MatrixCRS &operator=(MatrixCRS &&other);

    MatrixCRS &operator=(const Matrix< T_Scalar > &other) override;

    virtual void extract(const system::MatrixFactory< T_Scalar > *matrix) override;
    virtual void clear() override;

    MatrixFormat format() const override;

    virtual unsigned long long numberOfNonZero() const override;

    virtual bool isSymmetric(const scalar::Precision< T_Scalar > tolerance = 10. * scalar::Epsilon< T_Scalar >::value) const override;
    virtual bool isHermitian(const scalar::Precision< T_Scalar > tolerance = 10. * scalar::Epsilon< T_Scalar >::value) const override;

    virtual void save(const std::string &path) const override;
    virtual void saveSpyPlot(const std::string &path, const unsigned int pointSize = 1, const common::Color &zero = common::Color::white(), const common::Color &noZero = common::Color::black()) const override;

    const std::vector< unsigned long long > &ai() const;
    const std::vector< unsigned long long > &aj() const;
    const std::vector< T_Scalar > &a() const;
  };


} // namespace gmshfem::algebra

#endif // H_GMSHFEM_MATRIXCRS
