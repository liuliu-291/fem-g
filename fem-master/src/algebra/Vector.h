// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_VECTOR
#define H_GMSHFEM_VECTOR

#include "AlgebraicObject.h"
#include "gmshfemDefines.h"
#include "scalar.h"

#include <string>
#include <vector>

typedef struct _p_Vec *Vec;

namespace gmshfem::system
{
  template< class T_Scalar >
  class VectorFactory;
};

namespace gmshfem::algebra
{


  template< class T_Scalar >
  class Vector : public AlgebraicObject< T_Scalar >
  {
   protected:
    std::vector< T_Scalar > _vector;
#ifdef HAVE_PETSC
    Vec _vecPetsc;
#endif
    bool _havePetsc;
    std::vector< void * > _shouldBeDestroyedWithPetsc;

    virtual void _buildPetsc();

   public:
    Vector();
    Vector(const unsigned long long size);
    Vector(const Vector &other);
    Vector(Vector &&other);
    Vector(std::vector< T_Scalar > &&vector);
    virtual ~Vector();

    Vector &operator=(const Vector &other);
    Vector &operator=(Vector &&other);
    Vector &operator=(std::vector< T_Scalar > &&vector);
    Vector &operator+=(const Vector &other);
    Vector &operator-=(const Vector &other);
    Vector &operator*=(const T_Scalar& alpha);
    Vector &operator/=(const T_Scalar& alpha);
    Vector operator+(const Vector &other) const;
    Vector operator-(const Vector &other) const;
    Vector operator*(const T_Scalar& alpha) const;
    Vector operator/(const T_Scalar& alpha) const;
    T_Scalar operator*(const Vector &other) const;

    void extract(const system::VectorFactory< T_Scalar > *vector);

    T_Scalar &operator[](const unsigned long long pos);
    const T_Scalar &operator[](const unsigned long long pos) const;
    const std::vector< T_Scalar > &getStdVector() const;

    Vec getPetsc();
    void removePetsc();
    Eigen::VectorX< T_Scalar > getEigen();

    virtual AlgebraicObjectType type() const override;

    unsigned long long size() const;
    void resize(const unsigned long long size);
    void clear();

    void save(const std::string &path) const;

    void concatenate(const Vector< T_Scalar > &other);

    T_Scalar min() const;
    T_Scalar max() const;
    scalar::Precision< T_Scalar > norm(unsigned int p = 2) const;

  private:
    template< class V > using T_EigenBs = Eigen::MatrixBase< Eigen::Map< V > >;
    template< class V > using T_EigenOpVVV =
      Eigen::Map< V >&(T_EigenBs< V >::*)(const T_EigenBs< const V >&);
    template< class V > using T_EigenOpSVV =
      T_Scalar(T_EigenBs< const V >::*)(const T_EigenBs< const V >&) const;
    template< class V > using T_EigenOpVVS =
      Eigen::Map< V >&(Eigen::Map< V >::*)(const T_Scalar&);

    static void check(const Vector &lhs, const Vector &rhs, const std::string &name);
    static void eigenop(Vector &lhs, const Vector &rhs, const std::string &name,
                        T_EigenOpVVV< Eigen::VectorX< T_Scalar > > op);
    static T_Scalar eigenop(const Vector &lhs, const Vector &rhs, const std::string &name,
                            T_EigenOpSVV< Eigen::VectorX< T_Scalar > > op);
    static void eigenop(Vector &lhs, const T_Scalar &rhs, const std::string &name,
                        T_EigenOpVVS< Eigen::VectorX< T_Scalar > > op);
  };


} // namespace gmshfem::algebra

#endif // H_GMSHFEM_VECTOR
