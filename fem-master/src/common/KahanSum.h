// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_KAHANSUM
#define H_GMSHFEM_KAHANSUM

namespace gmshfem::common
{


  template< class T_Object >
  class KahanSum
  {
   private:
    T_Object _sum;
    T_Object _c;
    T_Object _y;
    T_Object _t;

   public:
    KahanSum() :
      _sum(), _c(), _y(), _t()
    {
    }

    KahanSum(const T_Object &b) :
      _sum(b), _c(b), _y(b), _t(b)
    {
    }

    ~KahanSum()
    {
    }

    T_Object sum() const
    {
      return _sum;
    }

    KahanSum &operator+=(const T_Object &b)
    {
      _y = b - _c;
      _t = _sum + _y;
      _c = (_t - _sum) - _y;
      _sum = _t;
      return *this;
    }
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_KAHANSUM
