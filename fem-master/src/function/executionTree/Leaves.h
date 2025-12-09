// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_LEAVES
#define H_GMSHFEM_LEAVES

#include "Exception.h"
#include "MathObject.h"

namespace gmshfem::function
{
  template< class T_Scalar, Degree T_Degree >
  class ExecutionTreeWithDegree;

  template< class T_Scalar >
  class ExecutionTree;
} // namespace gmshfem::function

namespace gmshfem::function
{


  template< class T_Scalar, unsigned int N, Degree T_Degree, Degree... TT_LeavesDegrees >
  class Leaves_t
  {
   private:
    Leaves_t< T_Scalar, N + 1, TT_LeavesDegrees... > _next;
    const ExecutionTreeWithDegree< T_Scalar, T_Degree > *_leaf;

   public:
    Leaves_t() :
      _next(), _leaf(nullptr) {}
    ~Leaves_t()
    {
      if(_leaf) {
        delete _leaf;
      }
    }

    unsigned int size() const
    {
      return _next.size();
    }

    const ExecutionTree< T_Scalar > *get(const unsigned int M) const
    {
      if(M == N) {
        return _leaf;
      }
      return _next.get(M);
    }

    void set(const unsigned int M, const ExecutionTree< T_Scalar > *leaf)
    {
      if(M == N) {
        _leaf = static_cast< const ExecutionTreeWithDegree< T_Scalar, T_Degree > * >(leaf);
      }
      _next.set(M, leaf);
    }
  };

  template< class T_Scalar, unsigned int N, Degree T_Degree >
  class Leaves_t< T_Scalar, N, T_Degree >
  {
   private:
    const ExecutionTreeWithDegree< T_Scalar, T_Degree > *_leaf;

   public:
    Leaves_t() :
      _leaf(nullptr) {}
    ~Leaves_t()
    {
      if(_leaf) {
        delete _leaf;
      }
    }

    unsigned int size() const
    {
      return N + 1;
    }

    const ExecutionTree< T_Scalar > *get(const unsigned int M) const
    {
      if(M > N) {
        throw common::Exception("Cannot access leaf " + std::to_string(M) + ": there are only " + std::to_string(N) + " leaf(ves)");
      }
      return _leaf;
    }

    void set(const unsigned int M, const ExecutionTree< T_Scalar > *leaf)
    {
      if(M > N) {
        throw common::Exception("Cannot access leaf " + std::to_string(M) + ": there are only " + std::to_string(N) + " leaf(ves)");
      }
      _leaf = static_cast< const ExecutionTreeWithDegree< T_Scalar, T_Degree > * >(leaf);
    }
  };

  template< class T_Scalar, Degree... TT_LeavesDegree >
  class Leaves
  {
   protected:
    Leaves_t< T_Scalar, 0, TT_LeavesDegree... > _leaves;

   public:
    Leaves() :
      _leaves(){};
    ~Leaves(){};

    unsigned int size() const
    {
      return _leaves.size();
    }

    const ExecutionTree< T_Scalar > *get(const unsigned int M) const
    {
      return _leaves.get(M);
    }

    void set(const unsigned int M, const ExecutionTree< T_Scalar > *leaf)
    {
      return _leaves.set(M, leaf);
    }
  };

  template< class T_Scalar >
  class Leaves< T_Scalar >
  {
   public:
    Leaves(){};
    ~Leaves(){};

    unsigned int size() const
    {
      return 0;
    }

    const ExecutionTree< T_Scalar > *get(const unsigned int M) const
    {
      return nullptr;
    }
  };


} // namespace gmshfem::function

#endif // H_GMSHFEM_LEAVES
