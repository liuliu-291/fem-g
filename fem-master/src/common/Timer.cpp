// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Timer.h"

#include "Message.h"
#include "OmpInterface.h"
#include "instantiate.h"

#include <ctime>
#include <string>

namespace gmshfem::common
{


  Timer::Timer() :
    _t1(0.), _t2(0.), _dt(0.) {}

  Timer::Timer(const Timer &other) :
    _t1(other._t1), _t2(other._t2), _dt(other._dt)
  {
  }

  Timer::Timer(Timer &&other) :
    _t1(other._t1), _t2(other._t2), _dt(other._dt)
  {
    other._t1 = other._t2 = other._dt = 0.;
  }

  Timer::~Timer() {}

  Timer &Timer::operator=(const Timer &other)
  {
    _t1 = other._t1;
    _t2 = other._t2;
    _dt = other._dt;

    return *this;
  }

  Timer &Timer::operator=(Timer &&other)
  {
    _t1 = other._t1;
    _t2 = other._t2;
    _dt = other._dt;
    other._t1 = other._t2 = other._dt = 0.;

    return *this;
  }

  void Timer::reset()
  {
    _t1 = 0.;
    _t2 = 0.;
    _dt = 0.;
  }

  void Timer::tick()
  {
    _t1 = omp::getTime();
  }

  void Timer::tock()
  {
    _t2 = omp::getTime();
    _dt += _t2 - _t1;
    _t1 = _t2;
  }

  double Timer::time() const
  {
    return _dt;
  }

  Timer &Timer::operator+=(const double &other)
  {
    _dt += other;
    return *this;
  }

  Timer &Timer::operator+=(const Timer &other)
  {
    _dt += other._dt;
    return *this;
  }

  Timer &Timer::operator-=(const double &other)
  {
    _dt -= other;
    return *this;
  }

  Timer &Timer::operator-=(const Timer &other)
  {
    _dt -= other._dt;
    return *this;
  }

  Timer &Timer::operator*=(const double &other)
  {
    _dt *= other;
    return *this;
  }

  Timer &Timer::operator*=(const Timer &other)
  {
    _dt *= other._dt;
    return *this;
  }

  Timer &Timer::operator/=(const double &other)
  {
    _dt /= other;
    return *this;
  }

  Timer &Timer::operator/=(const Timer &other)
  {
    _dt /= other._dt;
    return *this;
  }


  Timer operator+(const Timer &a, const Timer &b)
  {
    Timer ret = a;
    ret += b;
    return ret;
  }

  Timer operator-(const Timer &a, const Timer &b)
  {
    Timer ret = a;
    ret -= b;
    return ret;
  }

  template< MessageType T_Type >
  Message< T_Type > &operator<<(Message< T_Type > &message, const Timer &timer)
  {
    return (message << timer.time());
  }

  INSTANTIATE_OPLL(Message< TEMPLATE_PARAM_1 > &, , 0, 5, MessageType, TEMPLATE_ARGS(common::MessageType::None, common::MessageType::Info, common::MessageType::Warning, common::MessageType::Error, common::MessageType::Debug), TEMPLATE_PARAMS(Message< TEMPLATE_PARAM_1 > &, const Timer &))

  std::string day()
  {
    time_t timer;
    time(&timer);
    struct tm *info = localtime(&timer);
    switch(info->tm_mon) {
    case 0: return "January"; break;
    case 1: return "February"; break;
    case 2: return "March"; break;
    case 3: return "April"; break;
    case 4: return "May"; break;
    case 5: return "June"; break;
    case 6: return "July"; break;
    case 7: return "August"; break;
    case 8: return "September"; break;
    case 9: return "October"; break;
    case 10: return "November"; break;
    case 11: return "December"; break;
    default: break;
    }
    return "The end of time";
  }

  std::string today()
  {
    time_t timer;
    time(&timer);
    struct tm *info = localtime(&timer);
    std::string str = day() + " " + std::to_string(info->tm_mday) + ", " + std::to_string(info->tm_year + 1900);
    return str;
  }

  std::string hour()
  {
    time_t timer;
    time(&timer);
    struct tm *info = localtime(&timer);
    std::string str = (info->tm_hour < 10 ? "0" + std::to_string(info->tm_hour) : std::to_string(info->tm_hour)) + ":" + (info->tm_min < 10 ? "0" + std::to_string(info->tm_min) : std::to_string(info->tm_min)) + ":" + (info->tm_sec < 10 ? "0" + std::to_string(info->tm_sec) : std::to_string(info->tm_sec));
    return str;
  }


} // namespace gmshfem::common
