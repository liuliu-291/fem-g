// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MATRIX
#define H_GMSHFEM_MATRIX

#include "AlgebraicObject.h"
#include "Color.h"
#include "gmshfemDefines.h"
#include "scalar.h"

#include <string>
#include <vector>

typedef struct _p_Mat *Mat;

namespace gmshfem::system
{
  template< class T_Scalar >
  class MatrixFactory;
};

namespace gmshfem::algebra
{


  enum class MatrixFormat {
    CRSFast,
    CRS,
    CCS
  };

  template< class T_Scalar >
  class Matrix : public AlgebraicObject< T_Scalar >
  {
   protected:
#ifdef HAVE_PETSC
    Mat _matPetsc;
#endif
    bool _havePetsc;
    std::vector< void * > _shouldBeDestroyedWithPetsc;

    virtual void _buildPetsc() = 0;

   public:
    Matrix();
    Matrix(const unsigned long long size0);
    Matrix(const unsigned long long size0, const unsigned long long size1);
    Matrix(const Matrix &other);
    Matrix(Matrix &&other);
    virtual ~Matrix();

    virtual Matrix &operator=(const Matrix &other) = 0;

    void matrixCopy(const Matrix &other);
    void matrixCopy(Matrix &&other);

    virtual void extract(const system::MatrixFactory< T_Scalar > *matrix) = 0;
    virtual void clear() = 0;

    Mat getPetsc();
    void removePetsc();

    virtual AlgebraicObjectType type() const override;
    virtual MatrixFormat format() const = 0;

    unsigned long long size(const unsigned int dim) const;
    virtual unsigned long long numberOfNonZero() const = 0;
    double sparsity() const;
    bool isSquare() const;

    virtual bool isSymmetric(const scalar::Precision< T_Scalar > tolerance = 10. * scalar::Epsilon< T_Scalar >::value) const = 0;
    virtual bool isHermitian(const scalar::Precision< T_Scalar > tolerance = 10. * scalar::Epsilon< T_Scalar >::value) const = 0;

    virtual void save(const std::string &path) const = 0;
    virtual void saveSpyPlot(const std::string &path, const unsigned int pointSize = 1, const common::Color &zero = common::Color::white(), const common::Color &noZero = common::Color::black()) const = 0;
  };


} // namespace gmshfem::algebra

#include "MatrixCCS.h"
#include "MatrixCRS.h"
#include "MatrixCRSFast.h"

#endif // H_GMSHFEM_MATRIX
