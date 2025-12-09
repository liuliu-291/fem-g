// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "PGFPlotsio.h"

#include "Message.h"
#include "instantiate.h"

#include <cstdio>

namespace gmshfem::common
{


  //
  // class CurveInterface
  //

  CurveInterface::CurveInterface() :
    _color(Color::black()), _mark("none"), _name(), _lineWidth(0.5), _lineStyle()
  {
  }

  CurveInterface::~CurveInterface()
  {
  }

  void CurveInterface::setColor(const Color &color)
  {
    _color = color;
  }

  void CurveInterface::setMark(const std::string &mark)
  {
    _mark = mark;
  }

  void CurveInterface::setName(const std::string &name)
  {
    _name = name;
  }

  std::string CurveInterface::name() const
  {
    return _name;
  }

  void CurveInterface::setLineWidth(const double lw)
  {
    _lineWidth = lw;
  }

  void CurveInterface::setLineStyle(const std::string &ls)
  {
    _lineStyle = ls;
  }

  //
  // class Curve
  //

  template< class T_Scalar >
  Curve< T_Scalar >::Curve() :
    _x(), _y()
  {
  }

  template< class T_Scalar >
  Curve< T_Scalar >::~Curve()
  {
  }

  template< class T_Scalar >
  void Curve< T_Scalar >::setData(const std::vector< T_Scalar > &x, const std::vector< T_Scalar > &y)
  {
    _x = x;
    _y = y;
    if(_x.size() != _y.size()) {
      throw common::Exception("PGF Plots: x and y vectors are not of the same size");
    }
  }

  template< class T_Scalar >
  std::string Curve< T_Scalar >::addPlot() const
  {

    std::string str;
    str += "\t\t\\definecolor{color" + std::to_string(static_cast< unsigned int >(_color.r())) + std::to_string(static_cast< unsigned int >(_color.g())) + std::to_string(static_cast< unsigned int >(_color.b())) + "}{RGB}{" + std::to_string(static_cast< unsigned int >(_color.r())) + ", " + std::to_string(static_cast< unsigned int >(_color.g())) + ", " + std::to_string(static_cast< unsigned int >(_color.b())) + "};\n";

    str += "\t\t\\addplot[" + (_lineStyle.empty() ? "" : _lineStyle + ", ") +
           "color=color" + std::to_string(static_cast< unsigned int >(_color.r())) + std::to_string(static_cast< unsigned int >(_color.g())) + std::to_string(static_cast< unsigned int >(_color.b())) + ", mark=" + _mark + ", line width=" + std::to_string(_lineWidth) + "pt] coordinates {\n";
    for(auto i = 0ULL; i < _x.size(); ++i) {
      char x[64];
      char y[64];
      sprintf(x, "%g", _x[i]);
      sprintf(y, "%g", _y[i]);
      str += "\t\t\t\t(" + std::string(x) + ", " + std::string(y) + ")\n";
    }
    str += "\t\t\t\t};";
    return str;
  }

  INSTANTIATE_CLASS(Curve, 2, TEMPLATE_ARGS(double, float))


  //
  // class PGFPlotsio
  //

  PGFPlotsio::PGFPlotsio() :
    _file(), _xLabel(""), _yLabel(""), _axisType("axis"), _otherOptions(), _curves(), _printLegend(false), _legendPos("north east"), _grid(false)
  {
  }

  PGFPlotsio::PGFPlotsio(const std::string &path) :
    _file(), _xLabel(""), _yLabel(""), _axisType("axis"), _otherOptions(), _curves(), _printLegend(false), _legendPos("north east"), _grid(false)
  {
    open(path);
    if(!isOpen()) {
      throw common::Exception("Unable to open PGF Plots file " + path);
    }
  }

  PGFPlotsio::~PGFPlotsio()
  {
    close();
  }

  bool PGFPlotsio::open(const std::string &path)
  {
    close();
    try {
      _file.open(path + ".tex", std::fstream::out);
    }
    catch(...) {
      return false;
    }

    return true;
  }

  bool PGFPlotsio::isOpen() const
  {
    return _file.is_open();
  }

  void PGFPlotsio::close()
  {
    if(_file.is_open()) {
      _file.close();
    }
    for(auto i = 0ULL; i < _curves.size(); ++i) {
      delete _curves[i];
    }
    _curves.clear();
  }

  void PGFPlotsio::xLabel(const std::string &xLabel)
  {
    _xLabel = xLabel;
  }

  void PGFPlotsio::yLabel(const std::string &yLabel)
  {
    _yLabel = yLabel;
  }

  void PGFPlotsio::axisType(const std::string &axisType)
  {
    _axisType = axisType;
  }

  void PGFPlotsio::printLegend(const bool printLegend)
  {
    _printLegend = printLegend;
  }

  void PGFPlotsio::grid(const bool grid)
  {
    _grid = grid;
  }

  void PGFPlotsio::legendPos(const std::string &legendPos)
  {
    _legendPos = legendPos;
  }

  void PGFPlotsio::write(const bool readyToCompile)
  {
    if(!isOpen()) {
      msg::error << "Unable to write PGF Plots file: no file is opened" << msg::endl;
      return;
    }

    std::string indent = "";
    if(readyToCompile) {
      _file << "\\documentclass[12pt]{article}" << std::endl
            << std::endl;
      _file << "\\usepackage[utf8]{inputenc}" << std::endl;
      _file << "\\usepackage[T1]{fontenc}" << std::endl;
      _file << "\\usepackage[english]{babel}" << std::endl;
      _file << "\\usepackage{color}" << std::endl;
      _file << "\\usepackage{graphicx}" << std::endl;
      _file << "\\usepackage{tikz}" << std::endl;
      _file << "\\usepackage{pgfplots}" << std::endl
            << std::endl;

      _file << "\\begin{document}" << std::endl;

      indent = "\t";
    }

    _file << indent << "\\begin{tikzpicture}" << std::endl;

    _file << indent << "\t\\begin{" << _axisType << "}[" << std::endl;
    _file << indent << "\t\t\txlabel={" << _xLabel << "}," << std::endl;
    _file << indent << "\t\t\tylabel={" << _yLabel << "}," << std::endl;
    _file << indent << "\t\t\tlegend pos=" << _legendPos << "," << std::endl;
    _file << indent << "\t\t\tgrid=" << (_grid ? "both " : "none ");

    if(_otherOptions.size()) {
      for(auto iOption = 0ULL; iOption < _otherOptions.size(); ++iOption) {
        _file << "," << std::endl;
        _file << indent << "\t\t\t" << _otherOptions[iOption];
      }
      _file << "]" << std::endl;
    }
    else {
      _file << "]" << std::endl;
    }

    for(auto i = 0ULL; i < _curves.size(); ++i) {
      _file << _curves[i]->addPlot() << std::endl;
    }

    if(_printLegend) {
      _file << indent << "\t\t\\legend{";
      for(auto i = 0ULL; i < _curves.size(); ++i) {
        _file << _curves[i]->name() << (i == _curves.size() - 1 ? "}" : ", ");
      }
      _file << std::endl;
    }

    _file << indent << "\t\\end{" << _axisType << "}" << std::endl;
    _file << indent << "\\end{tikzpicture}" << std::endl;

    if(readyToCompile) {
      _file << "\\end{document}" << std::endl;
    }
  }

  void PGFPlotsio::addOtherOption(const std::string &option)
  {
    _otherOptions.push_back(option);
  }


} // namespace gmshfem::common
