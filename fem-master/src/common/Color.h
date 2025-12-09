// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_COLOR
#define H_GMSHFEM_COLOR

#include <utility>
#include <vector>

namespace gmshfem::common
{


  class Color
  {
   private:
    unsigned char _r;
    unsigned char _g;
    unsigned char _b;

   public:
    constexpr Color() :
      _r(0), _g(0), _b(0) {}
    constexpr Color(const unsigned char r, const unsigned char g, const unsigned char b) :
      _r(r), _g(g), _b(b) {}

    constexpr Color(const Color &other) :
      _r(other._r), _g(other._g), _b(other._b) {}
    constexpr Color(Color &&other) :
      _r(std::move(other._r)), _g(std::move(other._g)), _b(std::move(other._b)) {}

    Color &operator=(const Color &other);
    Color &operator=(Color &&other);

    bool operator==(const Color &other) const;
    bool operator!=(const Color &other) const;

    void r(unsigned char r);
    void g(unsigned char g);
    void b(unsigned char b);

    unsigned char r() const;
    unsigned char g() const;
    unsigned char b() const;

    void extractRGB(unsigned char *ptr) const;
    void extractGrey(unsigned char *ptr) const;
    void extractBW(unsigned char *ptr, const unsigned int pos) const;

    constexpr static Color apricot()
    {
      return Color(255, 215, 180);
    }
    constexpr static Color beige()
    {
      return Color(255, 250, 200);
    }
    constexpr static Color black()
    {
      return Color(0, 0, 0);
    }
    constexpr static Color blue()
    {
      return Color(0, 130, 200);
    }
    constexpr static Color brown()
    {
      return Color(170, 110, 40);
    }
    constexpr static Color cyan()
    {
      return Color(70, 240, 240);
    }
    constexpr static Color green()
    {
      return Color(60, 180, 75);
    }
    constexpr static Color grey()
    {
      return Color(128, 128, 128);
    }
    constexpr static Color lavender()
    {
      return Color(230, 190, 255);
    }
    constexpr static Color lime()
    {
      return Color(210, 245, 60);
    }
    constexpr static Color magenta()
    {
      return Color(240, 50, 230);
    }
    constexpr static Color maroon()
    {
      return Color(128, 0, 0);
    }
    constexpr static Color mint()
    {
      return Color(170, 255, 195);
    }
    constexpr static Color navy()
    {
      return Color(0, 0, 128);
    }
    constexpr static Color olive()
    {
      return Color(128, 128, 0);
    }
    constexpr static Color orange()
    {
      return Color(245, 130, 48);
    }
    constexpr static Color pink()
    {
      return Color(250, 190, 190);
    }
    constexpr static Color purple()
    {
      return Color(145, 30, 180);
    }
    constexpr static Color red()
    {
      return Color(230, 25, 75);
    }
    constexpr static Color teal()
    {
      return Color(0, 0, 128);
    }
    constexpr static Color white()
    {
      return Color(255, 255, 255);
    }
    constexpr static Color yellow()
    {
      return Color(255, 255, 25);
    }
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_COLOR
