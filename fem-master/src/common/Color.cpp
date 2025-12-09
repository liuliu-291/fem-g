// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Color.h"

namespace gmshfem::common
{


  Color &Color::operator=(const Color &other)
  {
    _r = other._r;
    _g = other._g;
    _b = other._b;
    return *this;
  }

  Color &Color::operator=(Color &&other)
  {
    _r = std::move(other._r);
    _g = std::move(other._g);
    _b = std::move(other._b);
    return *this;
  }

  bool Color::operator==(const Color &other) const
  {
    if(_r == other._r && _g == other._g && _b == other._b) {
      return true;
    }
    return false;
  }

  bool Color::operator!=(const Color &other) const
  {
    return !(*this == other);
  }

  void Color::r(unsigned char r)
  {
    _r = r;
  }

  void Color::g(unsigned char g)
  {
    _g = g;
  }

  void Color::b(unsigned char b)
  {
    _b = b;
  }

  unsigned char Color::r() const
  {
    return _r;
  }

  unsigned char Color::g() const
  {
    return _g;
  }

  unsigned char Color::b() const
  {
    return _b;
  }

  void Color::extractRGB(unsigned char *ptr) const
  {
    *ptr = _r;
    *(ptr + 1) = _g;
    *(ptr + 2) = _b;
  }

  void Color::extractGrey(unsigned char *ptr) const
  {
    *ptr = (_r + _g + _b) / 3;
  }

  void Color::extractBW(unsigned char *ptr, const unsigned int pos) const
  {
    if((_r + _g + _b) < 381) {
      *ptr |= (0x1 << (8 - pos - 1));
    }
  }


} // namespace gmshfem::common
