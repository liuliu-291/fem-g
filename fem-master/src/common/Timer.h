// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TIMER
#define H_GMSHFEM_TIMER

#include <string>

namespace gmshfem::common
{
  enum class MessageType;

  template< MessageType T_Type >
  class Message;
} // namespace gmshfem::common

namespace gmshfem::common
{


  class Timer
  {
   private:
    double _t1;
    double _t2;
    double _dt;

   public:
    Timer();
    Timer(const Timer &other);
    Timer(Timer &&other);
    ~Timer();

    Timer &operator=(const Timer &other);
    Timer &operator=(Timer &&other);

    void reset();

    void tick();
    void tock();

    double time() const;

    Timer &operator+=(const double &other);
    Timer &operator+=(const Timer &other);

    Timer &operator-=(const double &other);
    Timer &operator-=(const Timer &other);

    Timer &operator*=(const double &other);
    Timer &operator*=(const Timer &other);

    Timer &operator/=(const double &other);
    Timer &operator/=(const Timer &other);
  };

  Timer operator+(const Timer &a, const Timer &b);
  Timer operator-(const Timer &a, const Timer &b);

  template< MessageType T_Type >
  Message< T_Type > &operator<<(Message< T_Type > &message, const Timer &timer);

  std::string day();
  std::string today();
  std::string hour();


} // namespace gmshfem::common

#endif // H_GMSHFEM_TIMER
