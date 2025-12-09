// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_PGFPLOTSIO
#define H_GMSHFEM_PGFPLOTSIO

#include "Color.h"
#include "Exception.h"
#include "io.h"

#include <fstream>
#include <string>
#include <vector>

namespace gmshfem::common
{


  class CurveInterface
  {
   protected:
    Color _color;
    std::string _mark;
    std::string _name;
    double _lineWidth;
    std::string _lineStyle;

   public:
    CurveInterface();
    virtual ~CurveInterface();

    void setColor(const Color &color);
    void setMark(const std::string &mark);
    void setName(const std::string &name);
    std::string name() const;
    void setLineWidth(const double lw);
    void setLineStyle(const std::string &ls);

    virtual std::string addPlot() const = 0;
  };

  template< class T_Scalar >
  class Curve final : public CurveInterface
  {
   private:
    std::vector< T_Scalar > _x;
    std::vector< T_Scalar > _y;

   public:
    Curve();
    ~Curve();

    void setData(const std::vector< T_Scalar > &x, const std::vector< T_Scalar > &y);
    std::string addPlot() const override;
  };

  class PGFPlotsio
  {
   private:
    std::ofstream _file;
    std::string _xLabel;
    std::string _yLabel;
    std::string _axisType;
    std::vector< std::string > _otherOptions;
    std::vector< CurveInterface * > _curves;
    bool _printLegend;
    std::string _legendPos;
    bool _grid;

   public:
    PGFPlotsio();
    PGFPlotsio(const std::string &path);
    ~PGFPlotsio();

    bool open(const std::string &path);
    bool isOpen() const;
    void close();

    template< class T_Scalar >
    void addCurve(const Curve< T_Scalar > &curve)
    {
      Curve< T_Scalar > *curvePtr = new Curve< T_Scalar >(curve);
      _curves.push_back(curvePtr);
    }

    void xLabel(const std::string &xLabel);
    void yLabel(const std::string &yLabel);
    void axisType(const std::string &axisType);
    void printLegend(const bool printLegend);
    void legendPos(const std::string &legendPos);
    void grid(const bool grid);
    void addOtherOption(const std::string &option);

    void write(const bool readyToCompile = false);
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_PGFPLOTSIO
