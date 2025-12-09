// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_NULLARYOPERATIONS
#define H_GMSHFEM_NULLARYOPERATIONS

#include "FunctionSpaceInterface.h"
#include "Message.h"
#include "NullaryNode.h"
#include "OmpInterface.h"

#include <gmsh.h>

namespace gmshfem::function
{


  // Available nullary operations:
  //  BilinearInterpolation
  //  Constant
  //  LinearInterpolation
  //  Normal
  //  Phi
  //  ProbeView
  //  R2D
  //  R3D
  //  Tangent
  //  Theta
  //  TrinearInterpolation
  //  X
  //  Y
  //  Z

  //
  // BilinearInterpolation
  //
  /*
  The 2D table has the following elements
  (nx = number of values along x axis)
  (ny = number of values along y axis)

  data(x(0),y(0))  data(x(0),y(1))  data(x(0),y(2))  ... data(x(0),y(ny))
  data(x(1),y(0))  data(x(1),y(1))  data(x(1),y(2))  ... data(x(1),y(ny))
  data(x(2),y(0))  data(x(2),y(1))  data(x(2),y(2))  ... data(x(2),y(ny))
   ...
  data(x(nx),y(0)) data(x(nx),y(1)) data(x(nx),y(2)) ... data(x(nx),y(ny))
  */
  template< class T_Scalar, Degree T_Degree >
  class BilinearInterpolation final : public NullaryOperation< T_Scalar, T_Degree >
  {
   private:
    const unsigned int _nx;
    const unsigned int _ny;
    const std::vector< scalar::Precision< T_Scalar > > _pMin;
    const std::vector< scalar::Precision< T_Scalar > > _pMax;
    const std::vector< scalar::Precision< T_Scalar > > _x;
    const std::vector< scalar::Precision< T_Scalar > > _y;
    const std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > _data;
    const std::vector< scalar::Precision< T_Scalar > > *_px;
    const std::vector< scalar::Precision< T_Scalar > > *_py;
    const std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > *_pdata;
    const bool _method; // 0 = use _x and _y array, 1 = use _pMin and _pMax
    const bool _dataSource; // 0 = use vector, 1 = use pointer

   public:
    BilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > &data) :
      _nx(x.size()), _ny(y.size()), _pMin{x.front(), y.front()}, _pMax{x.back(), y.back()}, _x(x), _y(y), _data(data), _px(nullptr), _py(nullptr), _pdata(nullptr), _method(0), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("The x-vector size(" + std::to_string(_nx) + ") and the data size(" + std::to_string(_data.size()) + ") do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != _data[i].size()) {
          throw gmshfem::common::Exception("The y-vector size(" + std::to_string(_ny) + ") and the data size(" + std::to_string(_data[i].size()) + ") do not match");
        }
      }
    }

    BilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > &&data) :
      _nx(x.size()), _ny(y.size()), _pMin{x.front(), y.front()}, _pMax{x.back(), y.back()}, _x(std::move(x)), _y(std::move(y)), _data(std::move(data)), _px(nullptr), _py(nullptr), _pdata(nullptr), _method(0), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("The x-vector size and the data size do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != _data[i].size()) {
          throw gmshfem::common::Exception("The y-vector size and the data size do not match");
        }
      }
    }

    BilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > &data) :
      _nx(nx), _ny(ny), _pMin(pMin), _pMax(pMax), _x(), _y(), _data(data), _px(nullptr), _py(nullptr), _pdata(nullptr), _method(1), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("'nx' and the data size do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != _data[i].size()) {
          throw gmshfem::common::Exception("'ny' and the data size do not match");
        }
      }
    }

    BilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > &&data) :
      _nx(nx), _ny(ny), _pMin(pMin), _pMax(pMax), _x(), _y(), _data(std::move(data)), _px(nullptr), _py(nullptr), _pdata(nullptr), _method(1), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("'nx' and the data size do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != _data[i].size()) {
          throw gmshfem::common::Exception("'ny' and the data size do not match");
        }
      }
    }

    BilinearInterpolation(const BilinearInterpolation &other) :
      NullaryOperation< T_Scalar, T_Degree >(other), _nx(other._nx), _ny(other._ny), _pMin(other._pMin), _pMax(other._pMax), _x(other._x), _y(other._y), _data(other._data), _px(other._px), _py(other._py), _pdata(other._pdata), _method(other._method), _dataSource(other._dataSource)
    {
    }

    BilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > *data) :
      _nx(x->size()), _ny(y->size()), _pMin{x->front(), y->front()}, _pMax{x->back(), y->back()}, _x(), _y(), _data(), _px(x), _py(y), _pdata(data), _method(0), _dataSource(1)
    {
      if(_nx != _pdata->size()) {
        throw gmshfem::common::Exception("The x-vector size(" + std::to_string(_nx) + ") and the data size(" + std::to_string(_pdata->size()) + ") do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != (*_pdata)[i].size()) {
          throw gmshfem::common::Exception("The y-vector size(" + std::to_string(_ny) + ") and the data size(" + std::to_string((*_pdata)[i].size()) + ") do not match");
        }
      }
    }

    BilinearInterpolation(const unsigned int nx, const unsigned int ny, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > *data) :
      _nx(nx), _ny(ny), _pMin(pMin), _pMax(pMax), _x(), _y(), _data(), _px(nullptr), _py(nullptr), _pdata(data), _method(1), _dataSource(1)
    {
      if(_nx != _pdata->size()) {
        throw gmshfem::common::Exception("'nx' and the data size do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != (*_pdata)[i].size()) {
          throw gmshfem::common::Exception("'ny' and the data size do not match");
        }
      }
    }

    void operator()(OutputVector< T_Scalar, T_Degree > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      bool warningOutside = false;
      const std::vector< scalar::Precision< T_Scalar > > *px = (_dataSource ? _px : &_x);
      const std::vector< scalar::Precision< T_Scalar > > *py = (_dataSource ? _py : &_y);
      const std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > *pdata = (_dataSource ? _pdata : &_data);
      if(_method == 0) { // with _x and _y
        unsigned int i = 0, j = 0;
#pragma omp for
        for(auto n = 0ULL; n < values.size(); ++n) {
          const scalar::Precision< T_Scalar > xp = points[3 * n];
          const scalar::Precision< T_Scalar > yp = points[3 * n + 1];
          if(xp >= _pMin[0] && xp <= _pMax[0] && yp >= _pMin[1] && yp <= _pMax[1]) {

            // Interpolate point (xp,yp) in a regular grid
            // _x[i] < xp <= _x[i+1]
            // _y[j] < yp <= _y[j+1]
            if(!((*px)[i] <= xp && (*px)[i + 1] > xp)) {
              auto it = std::lower_bound(px->begin(), px->end(), xp);
              i = std::distance(px->begin(), it);
              i = i == 0 ? 0 : i - 1;
            }

            if(!((*py)[j] <= yp && (*py)[j + 1] > yp)) {
              auto it = std::lower_bound(py->begin(), py->end(), yp);
              j = std::distance(py->begin(), it);
              j = j == 0 ? 0 : j - 1;
            }

            const typename MathObject< T_Scalar, T_Degree >::Object d11 = (*pdata)[i][j];
            const typename MathObject< T_Scalar, T_Degree >::Object d21 = (*pdata)[i + 1][j];
            const typename MathObject< T_Scalar, T_Degree >::Object d12 = (*pdata)[i][j + 1];
            const typename MathObject< T_Scalar, T_Degree >::Object d22 = (*pdata)[i + 1][j + 1];
            const scalar::Precision< T_Scalar > x1 = (*px)[i];
            const scalar::Precision< T_Scalar > x2 = (*px)[i + 1];
            const scalar::Precision< T_Scalar > y1 = (*py)[j];
            const scalar::Precision< T_Scalar > y2 = (*py)[j + 1];

            values[n] = scalar::Precision< T_Scalar >(1.) / ((x2 - x1) * (y2 - y1)) * (d11 * (x2 - xp) * (y2 - yp) + d12 * (x2 - xp) * (yp - y1) + d21 * (xp - x1) * (y2 - yp) + d22 * (xp - x1) * (yp - y1));
          }
          else {
#pragma omp critical
            if(!warningOutside) {
              msg::debug << "Bilinear interpolation outside the data range" << msg::endl;
              warningOutside = true;
            }
          }
        }
      }
      else if(_method == 1) { // with _pMin and _pMax
        const scalar::Precision< T_Scalar > xStep = (_pMax[0] - _pMin[0]) / (_nx - 1);
        const scalar::Precision< T_Scalar > yStep = (_pMax[1] - _pMin[1]) / (_ny - 1);
        const scalar::Precision< T_Scalar > coeff = 1. / (xStep * yStep);
        const scalar::Precision< T_Scalar > invXStep = 1. / xStep;
        const scalar::Precision< T_Scalar > invYStep = 1. / yStep;
#pragma omp for
        for(auto n = 0ULL; n < values.size(); ++n) {
          const scalar::Precision< T_Scalar > xp = points[3 * n];
          const scalar::Precision< T_Scalar > yp = points[3 * n + 1];
          if(xp >= _pMin[0] && xp <= _pMax[0] && yp >= _pMin[1] && yp <= _pMax[1]) {

            // Interpolate point (xp,yp) in a regular grid
            // _x[i] <= xp < _x[i+1]
            // _y[j] <= yp < _y[j+1]
            unsigned int i = (xp - _pMin[0]) * invXStep;
            i = (i == _nx - 1 ? _nx - 2 : i);
            unsigned int j = (yp - _pMin[1]) * invYStep;
            j = (j == _ny - 1 ? _ny - 2 : j);

            //  d12    d22
            //  ********
            //  *      *
            //  *      *
            //  ********
            //  d11    d21
            const typename MathObject< T_Scalar, T_Degree >::Object d11 = (*pdata)[i][j];
            const typename MathObject< T_Scalar, T_Degree >::Object d21 = (*pdata)[i + 1][j];
            const typename MathObject< T_Scalar, T_Degree >::Object d12 = (*pdata)[i][j + 1];
            const typename MathObject< T_Scalar, T_Degree >::Object d22 = (*pdata)[i + 1][j + 1];
            const scalar::Precision< T_Scalar > x1 = i * xStep;
            const scalar::Precision< T_Scalar > x2 = (i + 1) * xStep;
            const scalar::Precision< T_Scalar > y1 = j * yStep;
            const scalar::Precision< T_Scalar > y2 = (j + 1) * yStep;

            values[n] = coeff * (d11 * (x2 - xp) * (y2 - yp) + d12 * (x2 - xp) * (yp - y1) + d21 * (xp - x1) * (y2 - yp) + d22 * (xp - x1) * (yp - y1));
          }
          else {
#pragma omp critical
            if(!warningOutside) {
              msg::debug << "Bilinear interpolation outside the data range" << msg::endl;
              warningOutside = true;
            }
          }
        }
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "bilinear interpolation";
    }

    bool operator==(const NullaryOperation< T_Scalar, T_Degree > &other) const override
    {
      if(_dataSource == 0) {
        if(_nx == static_cast< const BilinearInterpolation & >(other)._nx && _ny == static_cast< const BilinearInterpolation & >(other)._ny && _pMin == static_cast< const BilinearInterpolation & >(other)._pMin && _pMax == static_cast< const BilinearInterpolation & >(other)._pMax && _x == static_cast< const BilinearInterpolation & >(other)._x && _y == static_cast< const BilinearInterpolation & >(other)._y && _data == static_cast< const BilinearInterpolation & >(other)._data && _method == static_cast< const BilinearInterpolation & >(other)._method) {
          return true;
        }
      }
      else if(_dataSource == 1) {
        if(_nx == static_cast< const BilinearInterpolation & >(other)._nx && _ny == static_cast< const BilinearInterpolation & >(other)._ny && _pMin == static_cast< const BilinearInterpolation & >(other)._pMin && _pMax == static_cast< const BilinearInterpolation & >(other)._pMax && _px == static_cast< const BilinearInterpolation & >(other)._px && _py == static_cast< const BilinearInterpolation & >(other)._py && _pdata == static_cast< const BilinearInterpolation & >(other)._pdata && _method == static_cast< const BilinearInterpolation & >(other)._method) {
          return true;
        }
      }
      return false;
    }
  };

  //
  // Constant
  //

  template< class T_Scalar, Degree T_Degree >
  class Constant final : public NullaryOperation< T_Scalar, T_Degree >
  {
   private:
    const typename MathObject< T_Scalar, T_Degree >::Object _value;

   public:
    Constant(const typename MathObject< T_Scalar, T_Degree >::Object &value) :
      _value(value)
    {
    }

    Constant(const Constant &other) :
      NullaryOperation< T_Scalar, T_Degree >(other), _value(other._value)
    {
    }

    void operator()(OutputVector< T_Scalar, T_Degree > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = _value;
      }
    }

    bool isConstant() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "constant " + MathObject< T_Scalar, T_Degree >::to_string(_value);
    }

    bool operator==(const NullaryOperation< T_Scalar, T_Degree > &other) const override
    {
      if(_value == static_cast< const Constant & >(other)._value) {
        return true;
      }
      return false;
    }
  };

  template< class T_Scalar >
  class Constant< T_Scalar, Degree::Degree4 > final : public NullaryOperation< T_Scalar, Degree::Degree4 >
  {
   private:
    const typename MathObject< T_Scalar, Degree::Degree4 >::Object _value;

   public:
    Constant(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &value) :
      _value(value)
    {
    }

    Constant(const Constant &other) :
      NullaryOperation< T_Scalar, Degree::Degree4 >(other), _value(other._value)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree4 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = _value;
      }
    }

    bool isConstant() const override
    {
      return true;
    }

    std::string name() const override
    {
      return "constant " + MathObject< T_Scalar, Degree::Degree4 >::to_string(_value);
    }

    bool operator==(const NullaryOperation< T_Scalar, Degree::Degree4 > &other) const override
    {
      for(unsigned int i = 0; i < 3; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
          if(_value(i, j) != static_cast< const Constant & >(other)._value(i, j)) {
            return false;
          }
        }
      }
      return true;
    }
  };

  //
  // LinearInterpolation
  //
  /*
  The 1D table has the following elements
  (nx = number of values along x axis)

  data(x(0))  data(x(1))  data(x(2))  ... data(x(nx))
  */
  template< class T_Scalar, Degree T_Degree >
  class LinearInterpolation final : public NullaryOperation< T_Scalar, T_Degree >
  {
   private:
    const unsigned int _nx;
    const scalar::Precision< T_Scalar > _pMin;
    const scalar::Precision< T_Scalar > _pMax;
    const std::vector< scalar::Precision< T_Scalar > > _x;
    const std::vector< typename MathObject< T_Scalar, T_Degree >::Object > _data;
    const std::vector< scalar::Precision< T_Scalar > > *_px;
    const std::vector< typename MathObject< T_Scalar, T_Degree >::Object > *_pdata;
    const bool _method; // 0 = use _x and _y array, 1 = use _pMin and _pMax
    const bool _dataSource; // 0 = use vector, 1 = use pointer

   public:
    LinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< typename MathObject< T_Scalar, T_Degree >::Object > &data) :
      _nx(x.size()), _pMin(x.front()), _pMax(x.back()), _x(x), _data(data), _px(nullptr), _pdata(nullptr), _method(0), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("The x-vector size and the data size do not match");
      }
    }

    LinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< typename MathObject< T_Scalar, T_Degree >::Object > &&data) :
      _nx(x.size()), _pMin(x.front()), _pMax(x.back()), _x(std::move(x)), _data(std::move(data)), _px(nullptr), _pdata(nullptr), _method(0), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("The x-vector size and the data size do not match");
      }
    }

    LinearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, T_Degree >::Object > &data) :
      _nx(nx), _pMin(pMin), _pMax(pMax), _x(), _data(data), _px(nullptr), _pdata(nullptr), _method(1), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("'nx' and the data size do not match");
      }
    }

    LinearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, std::vector< typename MathObject< T_Scalar, T_Degree >::Object > &&data) :
      _nx(nx), _pMin(pMin), _pMax(pMax), _x(), _data(std::move(data)), _px(nullptr), _pdata(nullptr), _method(1), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("'nx' and the data size do not match");
      }
    }

    LinearInterpolation(const LinearInterpolation &other) :
      NullaryOperation< T_Scalar, T_Degree >(other), _nx(other._nx), _pMin(other._pMin), _pMax(other._pMax), _x(other._x), _data(other._data), _method(other._method), _dataSource(other._dataSource)
    {
    }

    LinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< typename MathObject< T_Scalar, T_Degree >::Object > *data) :
      _nx(x->size()), _pMin(x->front()), _pMax(x->back()), _x(), _data(), _px(x), _pdata(data), _method(0), _dataSource(1)
    {
      if(_nx != _pdata->size()) {
        throw gmshfem::common::Exception("The x-vector size and the data size do not match");
      }
    }

    LinearInterpolation(const unsigned int nx, const scalar::Precision< T_Scalar > &pMin, const scalar::Precision< T_Scalar > &pMax, const std::vector< typename MathObject< T_Scalar, T_Degree >::Object > *data) :
      _nx(nx), _pMin(pMin), _pMax(pMax), _x(), _data(), _px(nullptr), _pdata(data), _method(1), _dataSource(1)
    {
      if(_nx != _pdata->size()) {
        throw gmshfem::common::Exception("'nx' and the data size do not match");
      }
    }

    void operator()(OutputVector< T_Scalar, T_Degree > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      bool warningOutside = false;
      const std::vector< scalar::Precision< T_Scalar > > *px = (_dataSource ? _px : &_x);
      const std::vector< typename MathObject< T_Scalar, T_Degree >::Object > *pdata = (_dataSource ? _pdata : &_data);
      if(_method == 0) { // with _x and _y
        unsigned int i = 0;
#pragma omp for
        for(auto n = 0ULL; n < values.size(); ++n) {
          const scalar::Precision< T_Scalar > xp = points[3 * n];
          if(xp >= _pMin && xp <= _pMax) {

            // Interpolate point (xp) in a regular grid
            // _x[i] < xp <= _x[i+1]
            if(!((*px)[i] <= xp && (*px)[i + 1] > xp)) {
              auto it = std::lower_bound(px->begin(), px->end(), xp);
              i = std::distance(px->begin(), it);
              i = i == 0 ? 0 : i - 1;
            }

            const typename MathObject< T_Scalar, T_Degree >::Object d1 = (*pdata)[i];
            const typename MathObject< T_Scalar, T_Degree >::Object d2 = (*pdata)[i + 1];
            const scalar::Precision< T_Scalar > x1 = (*px)[i];
            const scalar::Precision< T_Scalar > x2 = (*px)[i + 1];

            values[n] = scalar::Precision< T_Scalar >(1.) / (x2 - x1) * (d1 * (x2 - xp) + d2 * (xp - x1));
          }
          else {
#pragma omp critical
            if(!warningOutside) {
              msg::debug << "Linear interpolation outside the data range" << msg::endl;
              warningOutside = true;
            }
          }
        }
      }
      else if(_method == 1) { // with _pMin and _pMax
        const scalar::Precision< T_Scalar > xStep = (_pMax - _pMin) / (_nx - 1);
        const scalar::Precision< T_Scalar > coeff = 1. / (xStep);
        const scalar::Precision< T_Scalar > invXStep = 1. / xStep;
#pragma omp for
        for(auto n = 0ULL; n < values.size(); ++n) {
          const scalar::Precision< T_Scalar > xp = points[3 * n];
          if(xp >= _pMin && xp <= _pMax) {

            // Interpolate point (xp) in a regular grid
            // _x[i] <= xp < _x[i+1]
            unsigned int i = (xp - _pMin) * invXStep;
            i = (i == _nx - 1 ? _nx - 2 : i);

            const typename MathObject< T_Scalar, T_Degree >::Object d1 = (*pdata)[i];
            const typename MathObject< T_Scalar, T_Degree >::Object d2 = (*pdata)[i + 1];
            const scalar::Precision< T_Scalar > x1 = i * xStep;
            const scalar::Precision< T_Scalar > x2 = (i + 1) * xStep;

            values[n] = coeff * (d1 * (x2 - xp) + d2 * (xp - x1));
          }
          else {
#pragma omp critical
            if(!warningOutside) {
              msg::debug << "Linear interpolation outside the data range" << msg::endl;
              warningOutside = true;
            }
          }
        }
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "linear interpolation";
    }

    bool operator==(const NullaryOperation< T_Scalar, T_Degree > &other) const override
    {
      if(_dataSource == 0) {
        if(_nx == static_cast< const LinearInterpolation & >(other)._nx && _pMin == static_cast< const LinearInterpolation & >(other)._pMin && _pMax == static_cast< const LinearInterpolation & >(other)._pMax && _x == static_cast< const LinearInterpolation & >(other)._x && _data == static_cast< const LinearInterpolation & >(other)._data && _method == static_cast< const LinearInterpolation & >(other)._method) {
          return true;
        }
      }
      else if(_dataSource == 1) {
        if(_nx == static_cast< const LinearInterpolation & >(other)._nx && _pMin == static_cast< const LinearInterpolation & >(other)._pMin && _pMax == static_cast< const LinearInterpolation & >(other)._pMax && _px == static_cast< const LinearInterpolation & >(other)._px && _pdata == static_cast< const LinearInterpolation & >(other)._pdata && _method == static_cast< const LinearInterpolation & >(other)._method) {
          return true;
        }
      }
      return false;
    }
  };

  //
  // Normal
  //

  template< class T_Scalar >
  class Normal final : public NullaryOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    mutable unsigned int _nbrNodesByElements;
    mutable std::vector< scalar::Precision< T_Scalar > > _nodesCoord;
    mutable std::vector< scalar::Precision< T_Scalar > > _functions;
    mutable std::vector< double > _gmshJacobians, _gmshDeterminants, _gmshPoints;
    mutable std::vector< scalar::Precision< T_Scalar > > _jacobians;

   public:
    Normal()
    {
    }

    Normal(const Normal &other) :
      NullaryOperation< T_Scalar, Degree::Degree1 >(other)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      if(entity.first == 3) {
        msg::warning << "Cannot apply the normal operator on a volume" << msg::endl;
        return;
      }

      const unsigned int nbrOfElements = points.size() / gaussPoints.size();
#pragma omp single
      {
        std::vector< double > gmshNodesParametricCoord;
        std::vector< std::size_t > gmshNodesTags;
        std::vector< double > gmshNodesCoord;
        gmsh::model::mesh::getNodesByElementType(elementType, gmshNodesTags, gmshNodesCoord, gmshNodesParametricCoord, entity.second, false);
        scalar::move(_nodesCoord, gmshNodesCoord);
      }

#pragma omp single
      {
        field::FunctionSpaceLagrange< scalar::Precision< T_Scalar > > functionSpace;
        _nbrNodesByElements = functionSpace.getBasisFunction(true, _functions, gaussPoints, elementType);
      }
      const unsigned int nbrGaussPoints = gaussPoints.size() / 3;

      if(entity.first == 1) {
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            MathObject< T_Scalar, Degree::Degree1 >::zero(values[i * nbrGaussPoints + j]);
            for(auto k = 0U; k < _nbrNodesByElements; ++k) {
              values[i * nbrGaussPoints + j](0) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3] * _functions[j * _nbrNodesByElements * 3 + k * 3];
              values[i * nbrGaussPoints + j](1) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 1] * _functions[j * _nbrNodesByElements * 3 + k * 3];
              values[i * nbrGaussPoints + j](2) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 2] * _functions[j * _nbrNodesByElements * 3 + k * 3];
            }
            values[i * nbrGaussPoints + j].normalize();
            const T_Scalar x = values[i * nbrGaussPoints + j](0);
            const T_Scalar y = values[i * nbrGaussPoints + j](1);
            values[i * nbrGaussPoints + j](0) = y;
            values[i * nbrGaussPoints + j](1) = -x;
          }
        }
      }
      else if(entity.first == 2) {
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            typename MathObject< T_Scalar, Degree::Degree1 >::Object t1;
            MathObject< T_Scalar, Degree::Degree1 >::zero(t1);
            typename MathObject< T_Scalar, Degree::Degree1 >::Object t2;
            MathObject< T_Scalar, Degree::Degree1 >::zero(t2);
            for(auto k = 0U; k < _nbrNodesByElements; ++k) {
              t1(0) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3] * _functions[j * _nbrNodesByElements * 3 + k * 3];
              t2(0) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3] * _functions[j * _nbrNodesByElements * 3 + k * 3 + 1];
              t1(1) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 1] * _functions[j * _nbrNodesByElements * 3 + k * 3];
              t2(1) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 1] * _functions[j * _nbrNodesByElements * 3 + k * 3 + 1];
              t1(2) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 2] * _functions[j * _nbrNodesByElements * 3 + k * 3];
              t2(2) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 2] * _functions[j * _nbrNodesByElements * 3 + k * 3 + 1];
            }
            values[i * nbrGaussPoints + j](0) = t1(1) * t2(2) - t1(2) * t2(1);
            values[i * nbrGaussPoints + j](1) = t1(2) * t2(0) - t1(0) * t2(2);
            values[i * nbrGaussPoints + j](2) = t1(0) * t2(1) - t1(1) * t2(0);
            values[i * nbrGaussPoints + j].normalize();
          }
        }
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "normal";
    }
  };

  //
  // Phi
  //

  template< class T_Scalar >
  class Phi final : public NullaryOperation< T_Scalar, Degree::Degree0 >
  {
   public:
    Phi()
    {
    }

    Phi(const Phi &other) :
      NullaryOperation< T_Scalar, Degree::Degree0 >(other)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      // spherical angular coordinate (mathematics convention (r, theta, phi) with r is the radius, theta [0, 2pi] is the azimuthal angle is the xy plane and phi [0, pi] is the inclination polar angle)
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::atan2(std::sqrt(points[3 * i] * points[3 * i] + points[3 * i + 1] * points[3 * i + 1]), points[3 * i + 2]);
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "phi";
    }
  };

  //
  // ProbeView
  //

  template< class T_Scalar, Degree T_Degree >
  class ProbeView final : public NullaryOperation< T_Scalar, T_Degree >
  {
   private:
    const unsigned int _view;
    const unsigned int _step;
    const bool _grad;
    const int _dim;

   public:
    ProbeView(const unsigned int view, const unsigned int step, const bool grad, const int dim) :
      _view(view), _step(step), _grad(grad), _dim(dim)
    {
    }

    void operator()(OutputVector< T_Scalar, T_Degree > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
    }

    bool isConstant() const override
    {
      return false;
    }
  };

  // Specification for scalars
  template< class T_Scalar >
  class ProbeView< T_Scalar, Degree::Degree0 > final : public NullaryOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const unsigned int _view;
    const unsigned int _step;
    const int _dim;
    const bool _grad;
    const double _distanceMax;

   public:
    ProbeView(const unsigned int view, const unsigned int step, const int dim, const bool grad, const double distanceMax) :
      _view(view), _step(step), _dim(dim), _grad(grad), _distanceMax(distanceMax)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      std::vector< double > viewval;
      double distance;
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        gmsh::view::probe(_view, points[3 * i], points[3 * i + 1], points[3 * i + 2], viewval, distance, _step, -1, _grad, _distanceMax,
                          std::vector< double >(), std::vector< double >(), std::vector< double >(), _dim);
        if(viewval.empty()) {
          msg::warning << "No value found at x = " << points[3 * i] << " y = " << points[3 * i + 1] << " z = " << points[3 * i + 2] << msg::endl;
          values[i] = 0.;
        }
        else {
          values[i] = viewval[0];
        }
      }
    }

    bool isConstant() const override
    {
      return false;
    }
  };

  // Specification for vectors
  template< class T_Scalar >
  class ProbeView< T_Scalar, Degree::Degree1 > final : public NullaryOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    const unsigned int _view;
    const unsigned int _step;
    const int _dim;
    const bool _grad;
    const double _distanceMax;

   public:
    ProbeView(const unsigned int view, const unsigned int step, const int dim, const bool grad, const double distanceMax) :
      _view(view), _step(step), _dim(dim), _grad(grad), _distanceMax(distanceMax)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      std::vector< double > viewval;
      double distance;
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        gmsh::view::probe(_view, points[3 * i], points[3 * i + 1], points[3 * i + 2], viewval, distance, _step, -1, _grad, _distanceMax,
                          std::vector< double >(), std::vector< double >(), std::vector< double >(), _dim);
        if(viewval.empty()) {
          msg::warning << "No value found at x = " << points[3 * i] << " y = " << points[3 * i + 1] << " z = " << points[3 * i + 2] << msg::endl;
          values[i] << 0., 0., 0.;
        }
        else {
          values[i] << viewval[0], viewval[1], viewval[2];
        }
      }
    }

    bool isConstant() const override
    {
      return false;
    }
  };

  // Specification for tensors
  template< class T_Scalar >
  class ProbeView< T_Scalar, Degree::Degree2 > final : public NullaryOperation< T_Scalar, Degree::Degree2 >
  {
   private:
    const unsigned int _view;
    const unsigned int _step;
    const int _dim;
    const bool _grad;
    const double _distanceMax;

   public:
    ProbeView(const unsigned int view, const unsigned int step, const int dim, const bool grad, const double distanceMax) :
      _view(view), _step(step), _dim(dim), _grad(grad), _distanceMax(distanceMax)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree2 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      std::vector< double > viewval;
      double distance;
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        gmsh::view::probe(_view, points[3 * i], points[3 * i + 1], points[3 * i + 2], viewval, distance, _step, -1, _grad, _distanceMax,
                          std::vector< double >(), std::vector< double >(), std::vector< double >(), _dim);
        if(viewval.empty()) {
          msg::warning << "No value found at x = " << points[3 * i] << " y = " << points[3 * i + 1] << " z = " << points[3 * i + 2] << msg::endl;
          values[i] << 0., 0., 0., 0., 0., 0., 0., 0., 0.;
        }
        else {
          values[i] << viewval[0], viewval[1], viewval[2],
            viewval[3], viewval[4], viewval[5],
            viewval[6], viewval[7], viewval[8];
        }
      }
    }

    bool isConstant() const override
    {
      return false;
    }
  };

  //
  // R2D
  //

  template< class T_Scalar >
  class R2D final : public NullaryOperation< T_Scalar, Degree::Degree0 >
  {
   public:
    R2D()
    {
    }

    R2D(const R2D &other) :
      NullaryOperation< T_Scalar, Degree::Degree0 >(other)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::sqrt(points[3 * i] * points[3 * i] + points[3 * i + 1] * points[3 * i + 1]);
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "r2D";
    }
  };

  //
  // R3D
  //

  template< class T_Scalar >
  class R3D final : public NullaryOperation< T_Scalar, Degree::Degree0 >
  {
   public:
    R3D()
    {
    }

    R3D(const R3D &other) :
      NullaryOperation< T_Scalar, Degree::Degree0 >(other)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::sqrt(points[3 * i] * points[3 * i] + points[3 * i + 1] * points[3 * i + 1] + points[3 * i + 2] * points[3 * i + 2]);
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "r3D";
    }
  };

  //
  // Tangent
  //

  template< class T_Scalar >
  class Tangent final : public NullaryOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    const unsigned int _comp;
    mutable unsigned int _nbrNodesByElements;
    mutable std::vector< scalar::Precision< T_Scalar > > _nodesCoord;
    mutable std::vector< scalar::Precision< T_Scalar > > _functions;
    mutable std::vector< double > _gmshJacobians, _gmshDeterminants, _gmshPoints;
    mutable std::vector< scalar::Precision< T_Scalar > > _jacobians;

   public:
    Tangent(const unsigned int comp) :
      _comp(comp)
    {
    }

    Tangent(const Tangent &other) :
      NullaryOperation< T_Scalar, Degree::Degree1 >(other), _comp(other._comp)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      if(entity.first == 3) {
        msg::warning << "Cannot apply the tangent operator on a volume" << msg::endl;
        return;
      }

      const unsigned int nbrOfElements = points.size() / gaussPoints.size();
#pragma omp single
      {
        std::vector< double > gmshNodesParametricCoord;
        std::vector< std::size_t > gmshNodesTags;
        std::vector< double > gmshNodesCoord;
        gmsh::model::mesh::getNodesByElementType(elementType, gmshNodesTags, gmshNodesCoord, gmshNodesParametricCoord, entity.second, false);
        scalar::move(_nodesCoord, gmshNodesCoord);
      }

#pragma omp single
      {
        field::FunctionSpaceLagrange< scalar::Precision< T_Scalar > > functionSpace;
        _nbrNodesByElements = functionSpace.getBasisFunction(true, _functions, gaussPoints, elementType);
      }
      const unsigned int nbrGaussPoints = gaussPoints.size() / 3;

      if(entity.first == 1) {
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            MathObject< T_Scalar, Degree::Degree1 >::zero(values[i * nbrGaussPoints + j]);
            for(auto k = 0U; k < _nbrNodesByElements; ++k) {
              values[i * nbrGaussPoints + j](0) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3] * _functions[j * _nbrNodesByElements * 3 + k * 3];
              values[i * nbrGaussPoints + j](1) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 1] * _functions[j * _nbrNodesByElements * 3 + k * 3];
              values[i * nbrGaussPoints + j](2) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 2] * _functions[j * _nbrNodesByElements * 3 + k * 3];
            }
            values[i * nbrGaussPoints + j].normalize();
          }
        }
      }
      else if(entity.first == 2) {
        const unsigned int numThreads = omp::getNumThreads();
        const unsigned int myThreadID = omp::getThreadNum();
        std::vector< double > gmshGaussPoints;
        scalar::copy(gmshGaussPoints, gaussPoints);
#pragma omp single
        gmsh::model::mesh::preallocateJacobians(elementType, nbrGaussPoints, true, true, false, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second);
        gmsh::model::mesh::getJacobians(elementType, gmshGaussPoints, _gmshJacobians, _gmshDeterminants, _gmshPoints, entity.second, myThreadID, numThreads);
#pragma omp barrier
#pragma omp single
        scalar::move(_jacobians, _gmshJacobians);
#pragma omp for
        for(auto i = 0U; i < nbrOfElements; ++i) {
          for(auto j = 0U; j < nbrGaussPoints; ++j) {
            MathObject< T_Scalar, Degree::Degree1 >::zero(values[i * nbrGaussPoints + j]);
            const Eigen::Matrix3< scalar::Precision< T_Scalar > > Jinv = Eigen::Map< const Eigen::Matrix3< scalar::Precision< T_Scalar > > >(&_jacobians[i * nbrGaussPoints * 9 + j * 9]).transpose().inverse();
            for(auto k = 0U; k < _nbrNodesByElements; ++k) {
              const scalar::Precision< T_Scalar > fct = Jinv(_comp, 0) * _functions[j * _nbrNodesByElements * 3 + k * 3 + 0] + Jinv(_comp, 1) * _functions[j * _nbrNodesByElements * 3 + k * 3 + 1] + Jinv(_comp, 2) * _functions[j * _nbrNodesByElements * 3 + k * 3 + 2];
              values[i * nbrGaussPoints + j](0) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3] * fct;
              values[i * nbrGaussPoints + j](1) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 1] * fct;
              values[i * nbrGaussPoints + j](2) += _nodesCoord[i * _nbrNodesByElements * 3 + k * 3 + 2] * fct;
            }
            values[i * nbrGaussPoints + j].normalize();
          }
        }
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "tangent";
    }

    bool operator==(const NullaryOperation< T_Scalar, Degree::Degree1 > &other) const override
    {
      if(_comp == static_cast< const Tangent & >(other)._comp) {
        return true;
      }
      return false;
    }
  };

  //
  // Theta
  //

  template< class T_Scalar >
  class Theta final : public NullaryOperation< T_Scalar, Degree::Degree0 >
  {
   public:
    Theta()
    {
    }

    Theta(const Theta &other) :
      NullaryOperation< T_Scalar, Degree::Degree0 >(other)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      // polar angular coordinate
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = std::atan2(points[3 * i + 1], points[3 * i]);
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "theta";
    }
  };

  //
  // TrilinearInterpolation
  //
  /*
  The 3D table has the following elements
  (nx = number of values along x axis)
  (ny = number of values along y axis)
  (nz = number of values along z axis)

  data(x(0),y(0),z(0))  data(x(0),y(0),z(1))  data(x(0),y(0),z(2))  ... data(x(0),y(0),z(nz))
  data(x(0),y(1),z(0))  data(x(0),y(1),z(1))  data(x(0),y(1),z(2))  ... data(x(0),y(1),z(nz))
  data(x(0),y(2),z(0))  data(x(0),y(2),z(1))  data(x(0),y(2),z(2))  ... data(x(0),y(2),z(nz))
   ...
  data(x(0),y(ny),z(0))  data(x(0),y(ny),z(1))  data(x(0),y(ny),z(2))  ... data(x(0),y(ny),z(nz))
  data(x(1),y(0),z(0))  data(x(1),y(0),z(1))  data(x(1),y(0),z(2))  ... data(x(1),y(0),z(nz))
  data(x(1),y(1),z(0))  data(x(1),y(1),z(1))  data(x(1),y(1),z(2))  ... data(x(1),y(1),z(nz))
  data(x(1),y(2),z(0))  data(x(1),y(2),z(1))  data(x(1),y(2),z(2))  ... data(x(1),y(2),z(nz))
   ...
  data(x(nx),y(ny),z(0))  data(x(nx),y(ny),z(1))  data(x(nx),y(ny),z(2))  ... data(x(nx),y(ny),z(nz))
  */
  template< class T_Scalar, Degree T_Degree >
  class TrilinearInterpolation final : public NullaryOperation< T_Scalar, T_Degree >
  {
   private:
    const unsigned int _nx;
    const unsigned int _ny;
    const unsigned int _nz;
    const std::vector< scalar::Precision< T_Scalar > > _pMin;
    const std::vector< scalar::Precision< T_Scalar > > _pMax;
    const std::vector< scalar::Precision< T_Scalar > > _x;
    const std::vector< scalar::Precision< T_Scalar > > _y;
    const std::vector< scalar::Precision< T_Scalar > > _z;
    const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > > _data;
    const std::vector< scalar::Precision< T_Scalar > > *_px;
    const std::vector< scalar::Precision< T_Scalar > > *_py;
    const std::vector< scalar::Precision< T_Scalar > > *_pz;
    const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > > *_pdata;
    const bool _method; // 0 = use _x and _y array, 1 = use _pMin and _pMax
    const bool _dataSource; // 0 = use vector, 1 = use pointer

   public:
    TrilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > &x, const std::vector< scalar::Precision< T_Scalar > > &y, const std::vector< scalar::Precision< T_Scalar > > &z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > > &data) :
      _nx(x.size()), _ny(y.size()), _nz(z.size()), _pMin{x.front(), y.front(), z.front()}, _pMax{x.back(), y.back(), z.back()}, _x(x), _y(y), _z(z), _data(data), _px(nullptr), _py(nullptr), _pz(nullptr), _pdata(nullptr), _method(0), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("The x-vector size(" + std::to_string(_nx) + ") and the data size(" + std::to_string(_data.size()) + ") do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != _data[i].size()) {
          throw gmshfem::common::Exception("The y-vector size(" + std::to_string(_ny) + ") and the data size(" + std::to_string(_data[i].size()) + ") do not match");
        }

        for(auto j = 0U; j < _ny; ++j) {
          if(_nz != _data[i][j].size()) {
            throw gmshfem::common::Exception("The z-vector size(" + std::to_string(_nz) + ") and the data size(" + std::to_string(_data[i][j].size()) + ") do not match");
          }
        }
      }
    }

    TrilinearInterpolation(std::vector< scalar::Precision< T_Scalar > > &&x, std::vector< scalar::Precision< T_Scalar > > &&y, std::vector< scalar::Precision< T_Scalar > > &&z, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > > &&data) :
      _nx(x.size()), _ny(y.size()), _nz(z.size()), _pMin{x.front(), y.front(), z.front()}, _pMax{x.back(), y.back(), z.back()}, _x(std::move(x)), _y(std::move(y)), _z(std::move(z)), _data(std::move(data)), _px(nullptr), _py(nullptr), _pz(nullptr), _pdata(nullptr), _method(0), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("The x-vector size and the data size do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != _data[i].size()) {
          throw gmshfem::common::Exception("The y-vector size and the data size do not match");
        }

        for(auto j = 0U; j < _ny; ++j) {
          if(_nz != _data[i][j].size()) {
            throw gmshfem::common::Exception("The z-vector size and the data size do not match");
          }
        }
      }
    }

    TrilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > > &data) :
      _nx(nx), _ny(ny), _nz(nz), _pMin(pMin), _pMax(pMax), _x(), _y(), _z(), _data(data), _px(nullptr), _py(nullptr), _pz(nullptr), _pdata(nullptr), _method(1), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("'nx' and the data size do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != _data[i].size()) {
          throw gmshfem::common::Exception("'ny' and the data size do not match");
        }

        for(auto j = 0U; j < _ny; ++j) {
          if(_nz != _data[i][j].size()) {
            throw gmshfem::common::Exception("'nz' and the data size do not match");
          }
        }
      }
    }

    TrilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, std::vector< std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > > &&data) :
      _nx(nx), _ny(ny), _nz(nz), _pMin(pMin), _pMax(pMax), _x(), _y(), _z(), _data(std::move(data)), _px(nullptr), _py(nullptr), _pz(nullptr), _pdata(nullptr), _method(1), _dataSource(0)
    {
      if(_nx != _data.size()) {
        throw gmshfem::common::Exception("'nx' and the data size do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != _data[i].size()) {
          throw gmshfem::common::Exception("'ny' and the data size do not match");
        }

        for(auto j = 0U; j < _ny; ++j) {
          if(_nz != _data[i][j].size()) {
            throw gmshfem::common::Exception("'nz' and the data size do not match");
          }
        }
      }
    }

    TrilinearInterpolation(const TrilinearInterpolation &other) :
      NullaryOperation< T_Scalar, T_Degree >(other), _nx(other._nx), _ny(other._ny), _nz(other._nz), _pMin(other._pMin), _pMax(other._pMax), _x(other._x), _y(other._y), _z(other._z), _data(other._data), _px(other._px), _py(other._py), _pz(other._pz), _pdata(other._pdata), _method(other._method), _dataSource(other._dataSource)
    {
    }

    TrilinearInterpolation(const std::vector< scalar::Precision< T_Scalar > > *x, const std::vector< scalar::Precision< T_Scalar > > *y, const std::vector< scalar::Precision< T_Scalar > > *z, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > > *data) :
      _nx(x->size()), _ny(y->size()), _nz(z->size()), _pMin{x->front(), y->front(), z->front()}, _pMax{x->back(), y->back(), z->back()}, _x(), _y(), _z(), _data(), _px(x), _py(y), _pz(z), _pdata(data), _method(0), _dataSource(1)
    {
      if(_nx != _pdata->size()) {
        throw gmshfem::common::Exception("The x-vector size(" + std::to_string(_nx) + ") and the data size(" + std::to_string(_pdata->size()) + ") do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != (*_pdata)[i].size()) {
          throw gmshfem::common::Exception("The y-vector size(" + std::to_string(_ny) + ") and the data size(" + std::to_string((*_pdata)[i].size()) + ") do not match");
        }

        for(auto j = 0U; j < _ny; ++j) {
          if(_nz != (*_pdata)[i][j].size()) {
            throw gmshfem::common::Exception("The z-vector size(" + std::to_string(_nz) + ") and the data size(" + std::to_string((*_pdata)[i][j].size()) + ") do not match");
          }
        }
      }
    }

    TrilinearInterpolation(const unsigned int nx, const unsigned int ny, const unsigned int nz, const std::vector< scalar::Precision< T_Scalar > > &pMin, const std::vector< scalar::Precision< T_Scalar > > &pMax, const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > > *data) :
      _nx(nx), _ny(ny), _nz(nz), _pMin(pMin), _pMax(pMax), _x(), _y(), _z(), _data(), _px(nullptr), _py(nullptr), _pz(nullptr), _pdata(data), _method(1), _dataSource(1)
    {
      if(_nx != _pdata->size()) {
        throw gmshfem::common::Exception("'nx' and the data size do not match");
      }

      for(auto i = 0U; i < _nx; ++i) {
        if(_ny != (*_pdata)[i].size()) {
          throw gmshfem::common::Exception("'ny' and the data size do not match");
        }

        for(auto j = 0U; j < _ny; ++j) {
          if(_nz != (*_pdata)[i][j].size()) {
            throw gmshfem::common::Exception("'nz' and the data size do not match");
          }
        }
      }
    }

    void operator()(OutputVector< T_Scalar, T_Degree > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
      bool warningOutside = false;
      const std::vector< scalar::Precision< T_Scalar > > *px = (_dataSource ? _px : &_x);
      const std::vector< scalar::Precision< T_Scalar > > *py = (_dataSource ? _py : &_y);
      const std::vector< scalar::Precision< T_Scalar > > *pz = (_dataSource ? _pz : &_z);
      const std::vector< std::vector< std::vector< typename MathObject< T_Scalar, T_Degree >::Object > > > *pdata = (_dataSource ? _pdata : &_data);
      if(_method == 0) { // with _x _y and _z
        unsigned int i = 0, j = 0, k = 0;
#pragma omp for
        for(auto n = 0ULL; n < values.size(); ++n) {
          const scalar::Precision< T_Scalar > xp = points[3 * n];
          const scalar::Precision< T_Scalar > yp = points[3 * n + 1];
          const scalar::Precision< T_Scalar > zp = points[3 * n + 2];
          if(xp >= _pMin[0] && xp <= _pMax[0] && yp >= _pMin[1] && yp <= _pMax[1] && zp >= _pMin[2] && zp <= _pMax[2]) {

            // Interpolate point (xp,yp,zp) in a regular grid
            // _x[i] < xp <= _x[i+1]
            // _y[j] < yp <= _y[j+1]
            // _z[k] < zp <= _z[k+1]
            if(!((*px)[i] <= xp && (*px)[i + 1] > xp)) {
              auto it = std::lower_bound(px->begin(), px->end(), xp);
              i = std::distance(px->begin(), it);
              i = i == 0 ? 0 : i - 1;
            }

            if(!((*py)[j] <= yp && (*py)[j + 1] > yp)) {
              auto it = std::lower_bound(py->begin(), py->end(), yp);
              j = std::distance(py->begin(), it);
              j = j == 0 ? 0 : j - 1;
            }

            if(!((*pz)[k] <= zp && (*pz)[k + 1] > zp)) {
              auto it = std::lower_bound(pz->begin(), pz->end(), zp);
              k = std::distance(pz->begin(), it);
              k = k == 0 ? 0 : k - 1;
            }

            const typename MathObject< T_Scalar, T_Degree >::Object d111 = (*pdata)[i][j][k];
            const typename MathObject< T_Scalar, T_Degree >::Object d211 = (*pdata)[i + 1][j][k];
            const typename MathObject< T_Scalar, T_Degree >::Object d121 = (*pdata)[i][j + 1][k];
            const typename MathObject< T_Scalar, T_Degree >::Object d221 = (*pdata)[i + 1][j + 1][k];
            const typename MathObject< T_Scalar, T_Degree >::Object d112 = (*pdata)[i][j][k + 1];
            const typename MathObject< T_Scalar, T_Degree >::Object d212 = (*pdata)[i + 1][j][k + 1];
            const typename MathObject< T_Scalar, T_Degree >::Object d122 = (*pdata)[i][j + 1][k + 1];
            const typename MathObject< T_Scalar, T_Degree >::Object d222 = (*pdata)[i + 1][j + 1][k + 1];
            const scalar::Precision< T_Scalar > x1 = (*px)[i];
            const scalar::Precision< T_Scalar > x2 = (*px)[i + 1];
            const scalar::Precision< T_Scalar > y1 = (*py)[j];
            const scalar::Precision< T_Scalar > y2 = (*py)[j + 1];
            const scalar::Precision< T_Scalar > z1 = (*pz)[k];
            const scalar::Precision< T_Scalar > z2 = (*pz)[k + 1];

            values[n] = scalar::Precision< T_Scalar >(1.) / ((x2 - x1) * (y2 - y1) * (z2 - z1)) * (d111 * (x2 - xp) * (y2 - yp) * (z2 - zp) + d121 * (x2 - xp) * (yp - y1) * (z2 - zp) + d211 * (xp - x1) * (y2 - yp) * (z2 - zp) + d221 * (xp - x1) * (yp - y1) * (z2 - zp) + d112 * (x2 - xp) * (y2 - yp) * (zp - z1) + d122 * (x2 - xp) * (yp - y1) * (zp - z1) + d212 * (xp - x1) * (y2 - yp) * (zp - z1) + d222 * (xp - x1) * (yp - y1) * (zp - z1));
          }
          else {
#pragma omp critical
            if(!warningOutside) {
              msg::debug << "Trilinear interpolation outside the data range" << msg::endl;
              warningOutside = true;
            }
          }
        }
      }
      else if(_method == 1) { // with _pMin and _pMax
        const scalar::Precision< T_Scalar > xStep = (_pMax[0] - _pMin[0]) / (_nx - 1);
        const scalar::Precision< T_Scalar > yStep = (_pMax[1] - _pMin[1]) / (_ny - 1);
        const scalar::Precision< T_Scalar > zStep = (_pMax[2] - _pMin[2]) / (_nz - 1);
        const scalar::Precision< T_Scalar > coeff = 1. / (xStep * yStep * zStep);
        const scalar::Precision< T_Scalar > invXStep = 1. / xStep;
        const scalar::Precision< T_Scalar > invYStep = 1. / yStep;
        const scalar::Precision< T_Scalar > invZStep = 1. / zStep;
#pragma omp for
        for(auto n = 0ULL; n < values.size(); ++n) {
          const scalar::Precision< T_Scalar > xp = points[3 * n];
          const scalar::Precision< T_Scalar > yp = points[3 * n + 1];
          const scalar::Precision< T_Scalar > zp = points[3 * n + 2];
          if(xp >= _pMin[0] && xp <= _pMax[0] && yp >= _pMin[1] && yp <= _pMax[1] && zp >= _pMin[2] && zp <= _pMax[2]) {

            // Interpolate point (xp,yp,zp) in a regular grid
            // _x[i] <= xp < _x[i+1]
            // _y[j] <= yp < _y[j+1]
            // _z[k] <= zp < _z[k+1]
            unsigned int i = (xp - _pMin[0]) * invXStep;
            i = (i == _nx - 1 ? _nx - 2 : i);
            unsigned int j = (yp - _pMin[1]) * invYStep;
            j = (j == _ny - 1 ? _ny - 2 : j);
            unsigned int k = (zp - _pMin[2]) * invZStep;
            k = (k == _nz - 1 ? _nz - 2 : k);

            const typename MathObject< T_Scalar, T_Degree >::Object d111 = (*pdata)[i][j][k];
            const typename MathObject< T_Scalar, T_Degree >::Object d211 = (*pdata)[i + 1][j][k];
            const typename MathObject< T_Scalar, T_Degree >::Object d121 = (*pdata)[i][j + 1][k];
            const typename MathObject< T_Scalar, T_Degree >::Object d221 = (*pdata)[i + 1][j + 1][k];
            const typename MathObject< T_Scalar, T_Degree >::Object d112 = (*pdata)[i][j][k + 1];
            const typename MathObject< T_Scalar, T_Degree >::Object d212 = (*pdata)[i + 1][j][k + 1];
            const typename MathObject< T_Scalar, T_Degree >::Object d122 = (*pdata)[i][j + 1][k + 1];
            const typename MathObject< T_Scalar, T_Degree >::Object d222 = (*pdata)[i + 1][j + 1][k + 1];
            const scalar::Precision< T_Scalar > x1 = i * xStep;
            const scalar::Precision< T_Scalar > x2 = (i + 1) * xStep;
            const scalar::Precision< T_Scalar > y1 = j * yStep;
            const scalar::Precision< T_Scalar > y2 = (j + 1) * yStep;
            const scalar::Precision< T_Scalar > z1 = k * zStep;
            const scalar::Precision< T_Scalar > z2 = (k + 1) * zStep;

            values[n] = coeff * (d111 * (x2 - xp) * (y2 - yp) * (z2 - zp) +
                                 d121 * (x2 - xp) * (yp - y1) * (z2 - zp) +
                                 d211 * (xp - x1) * (y2 - yp) * (z2 - zp) +
                                 d221 * (xp - x1) * (yp - y1) * (z2 - zp) +
                                 d112 * (x2 - xp) * (y2 - yp) * (zp - z1) +
                                 d122 * (x2 - xp) * (yp - y1) * (zp - z1) +
                                 d212 * (xp - x1) * (y2 - yp) * (zp - z1) +
                                 d222 * (xp - x1) * (yp - y1) * (zp - z1));
          }
          else {
#pragma omp critical
            if(!warningOutside) {
              msg::debug << "Trilinear interpolation outside the data range" << msg::endl;
              warningOutside = true;
            }
          }
        }
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "trilinear interpolation";
    }

    bool operator==(const NullaryOperation< T_Scalar, T_Degree > &other) const override
    {
      if(_dataSource == 0) {
        if(_nx == static_cast< const TrilinearInterpolation & >(other)._nx &&
           _ny == static_cast< const TrilinearInterpolation & >(other)._ny &&
           _nz == static_cast< const TrilinearInterpolation & >(other)._nz &&
           _pMin == static_cast< const TrilinearInterpolation & >(other)._pMin &&
           _pMax == static_cast< const TrilinearInterpolation & >(other)._pMax &&
           _x == static_cast< const TrilinearInterpolation & >(other)._x &&
           _y == static_cast< const TrilinearInterpolation & >(other)._y &&
           _z == static_cast< const TrilinearInterpolation & >(other)._z &&
           _data == static_cast< const TrilinearInterpolation & >(other)._data &&
           _method == static_cast< const TrilinearInterpolation & >(other)._method) {
          return true;
        }
      }
      else if(_dataSource == 1) {
        if(_nx == static_cast< const TrilinearInterpolation & >(other)._nx &&
           _ny == static_cast< const TrilinearInterpolation & >(other)._ny &&
           _nz == static_cast< const TrilinearInterpolation & >(other)._nz &&
           _pMin == static_cast< const TrilinearInterpolation & >(other)._pMin &&
           _pMax == static_cast< const TrilinearInterpolation & >(other)._pMax &&
           _px == static_cast< const TrilinearInterpolation & >(other)._px &&
           _py == static_cast< const TrilinearInterpolation & >(other)._py &&
           _pz == static_cast< const TrilinearInterpolation & >(other)._pz &&
           _pdata == static_cast< const TrilinearInterpolation & >(other)._pdata &&
           _method == static_cast< const TrilinearInterpolation & >(other)._method) {
          return true;
        }
      }
      return false;
    }
  };

  //
  // X
  //

  template< class T_Scalar >
  class X final : public NullaryOperation< T_Scalar, Degree::Degree0 >
  {
   public:
    X()
    {
    }

    X(const X &other) :
      NullaryOperation< T_Scalar, Degree::Degree0 >(other)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = points[3 * i];
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "x";
    }
  };

  //
  // Y
  //

  template< class T_Scalar >
  class Y final : public NullaryOperation< T_Scalar, Degree::Degree0 >
  {
   public:
    Y()
    {
    }

    Y(const Y &other) :
      NullaryOperation< T_Scalar, Degree::Degree0 >(other)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = points[3 * i + 1];
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "y";
    }
  };

  //
  // Z
  //

  template< class T_Scalar >
  class Z final : public NullaryOperation< T_Scalar, Degree::Degree0 >
  {
   public:
    Z()
    {
    }

    Z(const Z &other) :
      NullaryOperation< T_Scalar, Degree::Degree0 >(other)
    {
    }

    void operator()(OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const
    {
#pragma omp for
      for(auto i = 0ULL; i < values.size(); ++i) {
        values[i] = points[3 * i + 2];
      }
    }

    bool isConstant() const override
    {
      return false;
    }

    std::string name() const override
    {
      return "z";
    }
  };


} // namespace gmshfem::function


#endif // H_GMSHFEM_NULLARYOPERATIONS
