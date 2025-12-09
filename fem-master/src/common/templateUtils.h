// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TEMPLATEUTILS
#define H_GMSHFEM_TEMPLATEUTILS

#include <tuple>
#include <vector>

namespace gmshfem::common
{


  template< auto T_Begin, auto T_End, auto T_Inc, class T_F >
  constexpr void constexpr_for(T_F &&f)
  {
    if constexpr(T_Begin < T_End) {
      f(std::integral_constant< decltype(T_Begin), T_Begin >());
      constexpr_for< T_Begin + T_Inc, T_End, T_Inc >(f);
    }
  }

  template< int T_Begin, int T_End, class T_F, class... T_Args >
  constexpr void turn_into_tuple(T_F &&f, const std::vector< bool > &condition, T_Args... args)
  {
    if constexpr(T_Begin < T_End) {
      if(condition[T_Begin]) {
        turn_into_tuple< T_Begin + 1, T_End, T_F, T_Args..., std::true_type >(std::forward< T_F >(f), condition, args..., std::true_type());
      }
      else {
        turn_into_tuple< T_Begin + 1, T_End, T_F, T_Args..., std::false_type >(std::forward< T_F >(f), condition, args..., std::false_type());
      }
    }
    else {
      f(std::tuple(args...));
    }
  }


} // namespace gmshfem::common

#endif // H_GMSHFEM_TEMPLATEUTILS
