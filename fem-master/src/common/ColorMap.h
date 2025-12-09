// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_COLORMAP
#define H_GMSHFEM_COLORMAP

#include "Color.h"

namespace gmshfem::common
{


  enum class ColorMapName {
    Autumn,
    Blueberry,
    Desert,
    Forest,
    Fresh,
    Gmsh,
    Grey,
    Neutral,
    NightMountain,
    Spectrum,
    Summer,
    Thermal,
    Watermelon,
    Twilight
  };

  class ColorMap
  {
   protected:
    const ColorMapName _name;
    const unsigned int _numItem;

   public:
    ColorMap(const ColorMapName name, const unsigned int numItem);
    ~ColorMap();

    unsigned int size() const;
    Color operator()(const unsigned int item) const;
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_COLORMAP
