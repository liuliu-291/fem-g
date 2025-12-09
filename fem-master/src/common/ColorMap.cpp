// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "ColorMap.h"

#include "Color.h"
#include "Message.h"

#include <cmath>

namespace gmshfem::common
{


  //
  // class ColorMap
  //

  ColorMap::ColorMap(const ColorMapName name, const unsigned int numItem) :
    _name(name), _numItem(numItem)
  {
  }

  ColorMap::~ColorMap()
  {
  }

  unsigned int ColorMap::size() const
  {
    return _numItem;
  }

  Color ColorMap::operator()(const unsigned int item) const
  {
    auto cubicLow = [](const double c_r[], const double c_g[], const double c_b[], double &r, double &g, double &b, const double x) {
      for(auto i = 0; i < 4; ++i) {
        r += c_r[i] * std::pow(x, 3 - i);
        g += c_g[i] * std::pow(x, 3 - i);
        b += c_b[i] * std::pow(x, 3 - i);
      }
    };
    auto quarticLow = [](const double c_r[], const double c_g[], const double c_b[], double &r, double &g, double &b, const double x) {
      for(auto i = 0; i < 5; ++i) {
        r += c_r[i] * std::pow(x, 4 - i);
        g += c_g[i] * std::pow(x, 4 - i);
        b += c_b[i] * std::pow(x, 4 - i);
      }
    };
    switch(_name) {
    case ColorMapName::Autumn: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{-328.51, 409.5, 92.002, 46.};
      const double c_g[4]{715.5, -1089., 487.5, 35.};
      const double c_b[4]{-36., 40.499, -3.4997, 2.e-12};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::Blueberry: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{301.48, -179.98, 56.495, 30.};
      const double c_g[4]{251.98, -148.48, 90.495, 31.};
      const double c_b[4]{193.49, -152.98, 170.5, 38.};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::Desert: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{-670.5, 877.5, -113., 70.};
      const double c_g[4]{-220.5, 202.49, 41.003, 33.};
      const double c_b[4]{121.51, -310.51, 195.0, 26.};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::Forest: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{305.99, -310.48, 201.5, 4.};
      const double c_g[4]{148.49, -76.484, 105., 32.};
      const double c_b[4]{202.49, -125.98, 79.496, 44.};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::Fresh: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{2232., -3303., 1174., 76.};
      const double c_g[4]{1147.5, -1525.5, 390., 181.};
      const double c_b[4]{112.51, -234.02, -123.5, 245.};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::Gmsh: {
      const double x = ((double)item / (_numItem - 1)) * 128.;
      const double c_r = (x <= 46. ? 0. : (x >= 111. ? -0.03125 * (x - 111.) + 1. : (x >= 78. ? 1. : 0.03125 * (x - 46.))));
      const double c_g = (x <= 14. || x >= 111. ? 0. : (x >= 79. ? -0.03125 * (x - 111.) : (x <= 46. ? 0.03125 * (x - 14.) : 1.)));
      const double c_b = (x >= 79. ? 0. : (x >= 47. ? -0.03125 * (x - 79.) : (x <= 14. ? 0.03125 * (x - 14.) + 1. : 1.)));
      return Color((unsigned char)(c_r * 255), (unsigned char)(c_g * 255), (unsigned char)(c_b * 255));
    } break;
    case ColorMapName::Grey: {
      const double x = 1. - (double)item / (_numItem);
      return Color((unsigned char)(x * 255.), (unsigned char)(x * 255.), (unsigned char)(x * 255.));
    } break;
    case ColorMapName::Neutral: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{148.51, -558.02, 490.51, 98.};
      const double c_g[4]{351.02, -859.52, 535.51, 109.};
      const double c_b[4]{535.52, -1084.5, 539.01, 113.};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::NightMountain: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{-252.02, 607.52, -201.5, 54.};
      const double c_g[4]{44.992, 18.01, 36.998, 50.};
      const double c_b[4]{436.5, -661.5, 301., 55.};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::Spectrum: {
      const double x = (double)item / (_numItem - 1);
      double c_r = 0., c_g = 0., c_b = 0.;
      if(x >= 0. && x <= 0.2) {
        c_r = 1.;
        c_g = 5. * x;
        c_b = 0.;
      }
      else if(x > 0.2 && x <= 0.4) {
        c_r = -5. * (x - 0.2) + 1.;
        c_g = 1.;
        c_b = 0.;
      }
      else if(x > 0.4 && x <= 0.6) {
        c_r = 0.;
        c_g = 1.;
        c_b = 5 * (x - 0.4);
      }
      else if(x > 0.6 && x <= 0.8) {
        c_r = 0.;
        c_g = -5. * (x - 0.6) + 1.;
        c_b = 1.;
      }
      else if(x > 0.8 && x <= 1.) {
        c_r = 5. * (x - 0.8);
        c_g = 0.;
        c_b = 1.;
      }
      return Color((unsigned char)(c_r * 255), (unsigned char)(c_g * 255), (unsigned char)(c_b * 255));
    } break;
    case ColorMapName::Summer: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{207.02, -621.02, 418.01, 165.};
      const double c_g[4]{567.01, -882.01, 279., 195.};
      const double c_b[4]{580.51, -891.01, 163.5, 207.};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::Thermal: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{0.0054, -360.01, 600.01, 3.e-11};
      const double c_g[4]{-936.02, 1462.5, -302.5, 21.};
      const double c_b[4]{1579.5, -1746., 325.49, 96.};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::Watermelon: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[4]{-472.49, 494.99, -112.5, 245.};
      const double c_g[4]{-801., 1003.5, -47.496, 37.};
      const double c_b[4]{337.51, -715.52, 333., 73.};
      double r = 0., g = 0., b = 0.;
      cubicLow(c_r, c_g, c_b, r, g, b, x);
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    case ColorMapName::Twilight: {
      const double x = 1. - (double)item / (_numItem - 1);
      const double c_r[5]{2016., -5189.3, 4098., -695.67, 7.};
      const double c_g[5]{629.33, -1754.7, 1152.7, 55.667, 57.};
      const double c_b[5]{906.67, -1882.7, 823.33, 121.67, 121.};
      double r = 0., g = 0., b = 0.;
      quarticLow(c_r, c_g, c_b, r, g, b, x);
      if(r < 0.) r = 0.;
      return Color((unsigned char)r, (unsigned char)g, (unsigned char)b);
    } break;
    default:
      break;
    }
    return Color(0, 0, 0);
  }


} // namespace gmshfem::common
