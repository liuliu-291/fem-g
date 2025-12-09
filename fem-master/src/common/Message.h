// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MESSAGE
#define H_GMSHFEM_MESSAGE

#include "MathObject.h"
#include "scalar.h"

#include <complex>
#include <iostream>
#include <sstream>

namespace gmshfem::common
{
  class Timer;
}

namespace gmshfem::common
{


  enum class MessageType {
    None = 0,
    Error = 1,
    Warning = 2,
    Info = 3,
    Debug = 4
  };

  enum class MessageModificatorType {
    Endl = 0,
    Flush = 1,
    Fill = 2,
    Fixed = 3,
    Precision = 4,
    Scientific = 5
  };

  template< MessageType T_Type >
  class Message;

  template< MessageModificatorType T_MType >
  class MessageModificator
  {
   public:
    explicit MessageModificator();
    ~MessageModificator();

    template< MessageType T_Type >
    void run(Message< T_Type > &message) const;
  };

  template<>
  class MessageModificator< MessageModificatorType::Fill >
  {
   private:
    unsigned int _integer;
    char _character;

   public:
    explicit MessageModificator();
    ~MessageModificator();

    MessageModificator &operator()(const unsigned int integer, const char character = ' ');

    template< MessageType T_Type >
    void run(Message< T_Type > &message);
  };

  template<>
  class MessageModificator< MessageModificatorType::Precision >
  {
   private:
    unsigned int _integer;

   public:
    explicit MessageModificator();
    ~MessageModificator();

    MessageModificator &operator()(const unsigned int integer);

    template< MessageType T_Type >
    void run(Message< T_Type > &message);
  };


  template< MessageType T_Type >
  struct MessageInfo {
  };

  template<>
  struct MessageInfo< MessageType::None > {
    constexpr static const char *prefix = "";
    constexpr static const char *color = "";
    constexpr static const char *gmshLabel = "info";
    constexpr static std::ostream &stream = std::cout;
  };

  template<>
  struct MessageInfo< MessageType::Error > {
    constexpr static const char *prefix = "Error\t: ";
    constexpr static const char *color = "\33[1m\33[31m";
    constexpr static const char *gmshLabel = "error";
    constexpr static std::ostream &stream = std::cerr;
  };

  template<>
  struct MessageInfo< MessageType::Warning > {
    constexpr static const char *prefix = "Warning\t: ";
    constexpr static const char *color = "\33[35m";
    constexpr static const char *gmshLabel = "warning";
    constexpr static std::ostream &stream = std::cerr;
  };

  template<>
  struct MessageInfo< MessageType::Info > {
    constexpr static const char *prefix = "Info\t: ";
    constexpr static const char *color = "";
    constexpr static const char *gmshLabel = "info";
    constexpr static std::ostream &stream = std::cout;
  };

  template<>
  struct MessageInfo< MessageType::Debug > {
    constexpr static const char *prefix = "Debug\t: ";
    constexpr static const char *color = "\33[1m\33[33m";
    constexpr static const char *gmshLabel = "info";
    constexpr static std::ostream &stream = std::clog;
  };

  template< MessageType T_Type >
  class Message
  {
   private:
    std::ostringstream _stream;

   public:
    explicit Message();
    virtual ~Message();

    bool isAFile() const;
    virtual std::string getLevel() const;
    virtual std::string getPrefix() const;
    std::ostringstream &stream();

    template< MessageModificatorType T_MType >
    Message &operator<<(MessageModificator< T_MType > &modificator);

    virtual Message &operator<<(const char &data);
    virtual Message &operator<<(const std::string &data);

    virtual Message &operator<<(const short &data);
    virtual Message &operator<<(const unsigned short &data);
    virtual Message &operator<<(const int &data);
    virtual Message &operator<<(const unsigned int &data);
    virtual Message &operator<<(const long &data);
    virtual Message &operator<<(const unsigned long &data);
    virtual Message &operator<<(const long long &data);
    virtual Message &operator<<(const unsigned long long &data);

    virtual Message &operator<<(const float &data);
    virtual Message &operator<<(const double &data);
    virtual Message &operator<<(const long double &data);

    virtual Message &operator<<(const std::complex< float > &data);
    virtual Message &operator<<(const std::complex< double > &data);
    virtual Message &operator<<(const std::complex< long double > &data);

    virtual Message &operator<<(const typename MathObject< std::complex< double >, Degree::Degree1 >::Object &data);
    virtual Message &operator<<(const typename MathObject< double, Degree::Degree1 >::Object &data);
    virtual Message &operator<<(const typename MathObject< std::complex< float >, Degree::Degree1 >::Object &data);
    virtual Message &operator<<(const typename MathObject< float, Degree::Degree1 >::Object &data);

    virtual Message &operator<<(const typename MathObject< std::complex< double >, Degree::Degree2 >::Object &data);
    virtual Message &operator<<(const typename MathObject< double, Degree::Degree2 >::Object &data);
    virtual Message &operator<<(const typename MathObject< std::complex< float >, Degree::Degree2 >::Object &data);
    virtual Message &operator<<(const typename MathObject< float, Degree::Degree2 >::Object &data);

    virtual Message &operator<<(const typename MathObject< std::complex< double >, Degree::Degree4 >::Object &data);
    virtual Message &operator<<(const typename MathObject< double, Degree::Degree4 >::Object &data);
    virtual Message &operator<<(const typename MathObject< std::complex< float >, Degree::Degree4 >::Object &data);
    virtual Message &operator<<(const typename MathObject< float, Degree::Degree4 >::Object &data);
  };


} // namespace gmshfem::common

namespace gmshfem::msg
{


  extern common::MessageModificator< common::MessageModificatorType::Endl > endl;
  extern common::MessageModificator< common::MessageModificatorType::Flush > flush;
  extern common::MessageModificator< common::MessageModificatorType::Fill > fill;
  extern common::MessageModificator< common::MessageModificatorType::Fixed > fixed;
  extern common::MessageModificator< common::MessageModificatorType::Precision > precision;
  extern common::MessageModificator< common::MessageModificatorType::Scientific > scientific;

  extern common::Message< common::MessageType::None > print;
  extern common::Message< common::MessageType::Info > info;
  extern common::Message< common::MessageType::Warning > warning;
  extern common::Message< common::MessageType::Error > error;
  extern common::Message< common::MessageType::Debug > debug;

  unsigned int indent();
  unsigned int unindent();
  void setIndentation(const unsigned int level);


} // namespace gmshfem::msg

#endif // H_GMSHFEM_MESSAGE
