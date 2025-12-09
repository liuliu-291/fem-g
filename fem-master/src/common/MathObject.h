// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MATHOBJECT
#define H_GMSHFEM_MATHOBJECT

#include "scalar.h"

namespace gmshfem
{


  enum class Degree {
    Empty = -1,
    Degree0 = 1, // scalar
    Degree1 = 3, // Vector
    Degree2 = 9, // Tensor
    Degree4 = 81 // Tensor (4)
  };

  // **********
  // MathObject
  // **********
  template< class T_Scalar, Degree T_Degree >
  struct MathObject {
  };

  template< class T >
  static std::string s_tostr(const T &t)
  {
    std::ostringstream s;
    s << t;
    return s.str();
  }

  template< class T_Scalar >
  struct MathObject< T_Scalar, Degree::Degree0 > {
    typedef T_Scalar Object;
    constexpr static const char *name = "scalar";
    static void zero(Object &obj)
    {
      obj = 0.;
    };
    static void copy(Object &obj, const Eigen::MatrixX< T_Scalar > &mat)
    {
      obj = mat(0);
    };
    static std::string to_string(const Object &obj)
    {
      return s_tostr(obj);
    };
  };

  template< class T_Scalar >
  struct MathObject< T_Scalar, Degree::Degree1 > {
    typedef Eigen::Vector3< T_Scalar > Object;
    constexpr static const char *name = "vector";
    static void zero(Object &obj)
    {
      obj = Eigen::Vector3< T_Scalar >::Zero(3);
    };
    static void copy(Object &obj, const Object &mat)
    {
      obj = mat;
    };
    static std::string to_string(const Object &obj)
    {
      return "[" + s_tostr(obj(0)) + ", " + s_tostr(obj(1)) + ", " + s_tostr(obj(2)) + "]";
    };
  };

  template< class T_Scalar >
  struct MathObject< T_Scalar, Degree::Degree2 > {
    typedef Eigen::Matrix3< T_Scalar > Object;
    constexpr static const char *name = "tensor";
    static void zero(Object &obj)
    {
      obj = Eigen::Matrix3< T_Scalar >::Zero(3, 3);
    };
    static void copy(Object &obj, const Object &mat)
    {
      obj = mat;
    };
    static void copy(Object &obj, const Eigen::MatrixX< T_Scalar > &ten)
    {
      obj << ten(0, 0), ten(1, 0), ten(2, 0), ten(3, 0), ten(4, 0), ten(5, 0), ten(6, 0), ten(7, 0), ten(8, 0);
    };
    static std::string to_string(const Object &obj)
    {
      return "[" + s_tostr(obj(0, 0)) + ", " + s_tostr(obj(0, 1)) + ", " + s_tostr(obj(0, 2)) + "; " + s_tostr(obj(1, 0)) + ", " + s_tostr(obj(1, 1)) + ", " + s_tostr(obj(1, 2)) + "; " + s_tostr(obj(2, 0)) + ", " + s_tostr(obj(2, 1)) + ", " + s_tostr(obj(2, 2)) + "]";
    };
  };

  template< class T_Scalar >
  struct MathObject< T_Scalar, Degree::Degree4 > {
    typedef Eigen::Tensor3< T_Scalar > Object;
    constexpr static const char *name = "tensor<4>";
    static void zero(Object &obj)
    {
      for(unsigned int i = 0; i < 3; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
          obj(i, j) = Eigen::Matrix3< T_Scalar >::Zero(3, 3);
        }
      }
    };
    static void copy(Object &obj, const Object &ten)
    {
      obj = ten;
    };
    static std::string to_string(const Object &obj)
    {
      std::string str = "[A, B, C; D, E, F; G, H, I]\n";
      const std::string mat[3][3] = {{"A", "B", "C"}, {"D", "E", "F"}, {"G", "H", "I"}};
      for(unsigned int i = 0; i < 3; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
          str += mat[i][j] + " = " + MathObject< T_Scalar, Degree::Degree2 >::to_string(obj(i, j)) + "\n";
        }
      }
      return str;
    };
  };

  // ***********
  // MathProduct
  // ***********
  template< Degree T_Degree1, Degree T_Degree2 >
  struct MathProduct {
  };

  // scalar
  template< Degree T_Degree2 >
  struct MathProduct< Degree::Degree0, T_Degree2 > {
    constexpr static Degree ProductDegree = T_Degree2;
  };

  // vector
  template< Degree T_Degree2 >
  struct MathProduct< Degree::Degree1, T_Degree2 > {
    constexpr static Degree ProductDegree = Degree::Degree1;
  };

  template<>
  struct MathProduct< Degree::Degree1, Degree::Degree1 > {
    constexpr static Degree ProductDegree = Degree::Degree0;
  };

  // tensor
  template< Degree T_Degree2 >
  struct MathProduct< Degree::Degree2, T_Degree2 > {
    constexpr static Degree ProductDegree = Degree::Degree2;
  };

  template<>
  struct MathProduct< Degree::Degree2, Degree::Degree1 > {
    constexpr static Degree ProductDegree = Degree::Degree1;
  };

  // tensor4
  template<>
  struct MathProduct< Degree::Degree4, Degree::Degree4 > {
    constexpr static Degree ProductDegree = Degree::Degree4;
  };
  
  template<>
  struct MathProduct< Degree::Degree4, Degree::Degree2 > {
    constexpr static Degree ProductDegree = Degree::Degree2;
  };

  template<>
  struct MathProduct< Degree::Degree4, Degree::Degree0 > {
    constexpr static Degree ProductDegree = Degree::Degree4;
  };

  template<>
  struct MathProduct< Degree::Degree0, Degree::Degree4 > {
    constexpr static Degree ProductDegree = Degree::Degree4;
  };


} // namespace gmshfem


#endif // H_GMSHFEM_MATHOBJECT
