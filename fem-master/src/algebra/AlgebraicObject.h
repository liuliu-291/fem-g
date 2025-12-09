// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_ALGEBRAICOBJECT
#define H_GMSHFEM_ALGEBRAICOBJECT

namespace gmshfem::algebra
{


  enum class AlgebraicObjectType {
    Matrix,
    Vector
  };

  template< class T_Scalar >
  class AlgebraicObject
  {
   protected:
    unsigned long long _size[2];

   public:
    AlgebraicObject();
    AlgebraicObject(const unsigned long long size0, const unsigned long long size1);
    AlgebraicObject(const AlgebraicObject &other);
    AlgebraicObject(AlgebraicObject &&other);
    ~AlgebraicObject();

    void copySize(const AlgebraicObject &other);
    void copySize(AlgebraicObject &&other);

    virtual AlgebraicObjectType type() const = 0;
  };


} // namespace gmshfem::algebra

#endif // H_GMSHFEM_ALGEBRAICOBJECT
