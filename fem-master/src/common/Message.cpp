// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Message.h"

#include "Options.h"
#include "Timer.h"
#include "instantiate.h"

#include <gmsh.h>
#include <iomanip>
#include <sys/stat.h>

namespace gmshfem::common
{


  static unsigned int s_indentation = 0;

  //
  // class MessageModificator
  //

  template< MessageModificatorType T_MType >
  MessageModificator< T_MType >::MessageModificator()
  {
  }

  template< MessageModificatorType T_MType >
  MessageModificator< T_MType >::~MessageModificator()
  {
  }

  template< MessageModificatorType T_MType >
  template< MessageType T_Type >
  void MessageModificator< T_MType >::run(Message< T_Type > &message) const
  {
    switch(T_MType) {
    case MessageModificatorType::Endl: {
      if(common::Options::instance()->interface) {
        if(T_Type != MessageType::Debug || common::Options::instance()->debug) {
          gmsh::logger::write(message.stream().str(), message.getLevel());
        }
      }
      else {
        if(T_Type != MessageType::Debug || common::Options::instance()->debug) {
          MessageInfo< T_Type >::stream << message.getPrefix() << message.stream().str() << (message.isAFile() ? "" : "\33[0m") << std::endl;
        }
      }
      message.stream().unsetf(std::ios::floatfield);
      message.stream().precision(6);
      message.stream() << std::boolalpha;
      message.stream().str("");
      message.stream().clear();
    } break;
    case MessageModificatorType::Flush: {
      if(common::Options::instance()->interface) {
        if(T_Type != MessageType::Debug || common::Options::instance()->debug) {
          gmsh::logger::write(message.stream().str(), message.getLevel());
        }
      }
      else {
        if(T_Type != MessageType::Debug || common::Options::instance()->debug) {
          MessageInfo< T_Type >::stream << message.getPrefix() << message.stream().str() << (message.isAFile() ? "" : "\33[0m") << std::flush;
        }
      }
      message.stream().unsetf(std::ios::floatfield);
      message.stream().precision(6);
      message.stream() << std::boolalpha;
    } break;
    case MessageModificatorType::Fixed: {
      message.stream() << std::fixed;
    } break;
    case MessageModificatorType::Scientific: {
      message.stream() << std::scientific;
    } break;
    }
  }

  // class MessageModificator< MessageModificatorType::Fill >

  MessageModificator< MessageModificatorType::Fill >::MessageModificator() :
    _integer(0), _character(' ')
  {
  }

  MessageModificator< MessageModificatorType::Fill >::~MessageModificator()
  {
  }

  MessageModificator< MessageModificatorType::Fill > &MessageModificator< MessageModificatorType::Fill >::operator()(unsigned int integer, const char character)
  {
    _integer = integer;
    _character = character;
    return *this;
  }

  template< MessageType T_Type >
  void MessageModificator< MessageModificatorType::Fill >::run(Message< T_Type > &message)
  {
    for(auto i = 0; _integer > static_cast< unsigned int >(message.stream().tellp()); ++i) {
      message.stream() << _character;
    }
    _integer = 0;
    _character = ' ';
  }

  // class MessageModificator< MessageModificatorType::Precision >

  MessageModificator< MessageModificatorType::Precision >::MessageModificator() :
    _integer(0)
  {
  }

  MessageModificator< MessageModificatorType::Precision >::~MessageModificator()
  {
  }

  MessageModificator< MessageModificatorType::Precision > &MessageModificator< MessageModificatorType::Precision >::operator()(unsigned int integer)
  {
    _integer = integer;
    return *this;
  }

  template< MessageType T_Type >
  void MessageModificator< MessageModificatorType::Precision >::run(Message< T_Type > &message)
  {
    message.stream().precision(_integer);
    _integer = 0;
  }


  //
  // class Message
  //

  template< MessageType T_Type >
  Message< T_Type >::Message() :
    _stream()
  {
  }

  template< MessageType T_Type >
  Message< T_Type >::~Message()
  {
  }

  template< MessageType T_Type >
  bool Message< T_Type >::isAFile() const
  {
    struct stat streamStat;
    if(fstat(fileno(stderr), &streamStat) == 0) {
      if(streamStat.st_mode & S_IFREG) {
        return true;
      }
    }
    return false;
  }

  template< MessageType T_Type >
  std::string Message< T_Type >::getLevel() const
  {
    return std::string(MessageInfo< T_Type >::gmshLabel);
  }

  template< MessageType T_Type >
  std::string Message< T_Type >::getPrefix() const
  {
    std::string prefix = (isAFile() ? "" : std::string(MessageInfo< T_Type >::color)) + std::string(MessageInfo< T_Type >::prefix);
    if(!common::Options::instance()->interface) {
      for(auto i = 0U; i < s_indentation; ++i) {
        prefix += "\t";
      }
    }
    return prefix;
  }

  template< MessageType T_Type >
  std::ostringstream &Message< T_Type >::stream()
  {
    return _stream;
  }

  template< MessageType T_Type >
  template< MessageModificatorType T_MType >
  Message< T_Type > &Message< T_Type >::operator<<(MessageModificator< T_MType > &modificator)
  {
    if(static_cast< int >(T_Type) > common::Options::instance()->verbose && !common::Options::instance()->debug) {
      return *this;
    }

    modificator.run(*this);
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const char &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const std::string &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const short &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const unsigned short &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const int &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const unsigned int &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const long &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const unsigned long &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const long long &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const unsigned long long &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const float &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const double &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const long double &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      _stream << data;
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const std::complex< float > &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      if(data.imag() == 0.) {
        _stream << data.real();
      }
      else {
        _stream << "(" << data.real() << "," << data.imag() << ")";
      }
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const std::complex< double > &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      if(data.imag() == 0.) {
        _stream << data.real();
      }
      else {
        _stream << "(" << data.real() << "," << data.imag() << ")";
      }
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const std::complex< long double > &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      if(data.imag() == 0.) {
        _stream << data.real();
      }
      else {
        _stream << "(" << data.real() << "," << data.imag() << ")";
      }
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< std::complex< double >, Degree::Degree1 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ " << data(0) << ", " << data(1) << ", " << data(2) << " ]";
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< double, Degree::Degree1 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ " << data(0) << ", " << data(1) << ", " << data(2) << " ]";
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< std::complex< float >, Degree::Degree1 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ " << data(0) << ", " << data(1) << ", " << data(2) << " ]";
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< float, Degree::Degree1 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ " << data(0) << ", " << data(1) << ", " << data(2) << " ]";
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< std::complex< double >, Degree::Degree2 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ " << data(0, 0) << ", " << data(0, 1) << ", " << data(0, 2) << "; "
          << data(1, 0) << ", " << data(1, 1) << ", " << data(1, 2) << "; "
          << data(2, 0) << ", " << data(2, 1) << ", " << data(2, 2) << " ]";
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< double, Degree::Degree2 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ " << data(0, 0) << ", " << data(0, 1) << ", " << data(0, 2) << "; "
          << data(1, 0) << ", " << data(1, 1) << ", " << data(1, 2) << "; "
          << data(2, 0) << ", " << data(2, 1) << ", " << data(2, 2) << " ]";
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< std::complex< float >, Degree::Degree2 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ " << data(0, 0) << ", " << data(0, 1) << ", " << data(0, 2) << "; "
          << data(1, 0) << ", " << data(1, 1) << ", " << data(1, 2) << "; "
          << data(2, 0) << ", " << data(2, 1) << ", " << data(2, 2) << " ]";
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< float, Degree::Degree2 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ " << data(0, 0) << ", " << data(0, 1) << ", " << data(0, 2) << "; "
          << data(1, 0) << ", " << data(1, 1) << ", " << data(1, 2) << "; "
          << data(2, 0) << ", " << data(2, 1) << ", " << data(2, 2) << " ]";
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< std::complex< double >, Degree::Degree4 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ "
          << "A"
          << ", "
          << "B"
          << ", "
          << "C"
          << "; "
          << "D"
          << ", "
          << "E"
          << ", "
          << "F"
          << "; "
          << "G"
          << ", "
          << "H"
          << ", "
          << "I"
          << " ]" << std::endl;
      const std::string subMat[9] = {"A", "B", "C", "D", "E", "F", "G", "H", "I"};
      for(unsigned int i = 0; i < 3; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
          str << "\t " << subMat[i * 3 + j] << " = "
              << "[ " << data(i, j)(0, 0) << ", " << data(i, j)(0, 1) << ", " << data(i, j)(0, 2) << "; " << data(i, j)(1, 0) << ", " << data(i, j)(1, 1) << ", " << data(i, j)(1, 2) << "; " << data(i, j)(2, 0) << ", " << data(i, j)(2, 1) << ", " << data(i, j)(2, 2) << " ]" << std::endl;
          ;
        }
      }
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< double, Degree::Degree4 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ "
          << "A"
          << ", "
          << "B"
          << ", "
          << "C"
          << "; "
          << "D"
          << ", "
          << "E"
          << ", "
          << "F"
          << "; "
          << "G"
          << ", "
          << "H"
          << ", "
          << "I"
          << " ]" << std::endl;
      const std::string subMat[9] = {"A", "B", "C", "D", "E", "F", "G", "H", "I"};
      for(unsigned int i = 0; i < 3; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
          str << "\t " << subMat[i * 3 + j] << " = "
              << "[ " << data(i, j)(0, 0) << ", " << data(i, j)(0, 1) << ", " << data(i, j)(0, 2) << "; " << data(i, j)(1, 0) << ", " << data(i, j)(1, 1) << ", " << data(i, j)(1, 2) << "; " << data(i, j)(2, 0) << ", " << data(i, j)(2, 1) << ", " << data(i, j)(2, 2) << " ]" << std::endl;
          ;
        }
      }
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< std::complex< float >, Degree::Degree4 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ "
          << "A"
          << ", "
          << "B"
          << ", "
          << "C"
          << "; "
          << "D"
          << ", "
          << "E"
          << ", "
          << "F"
          << "; "
          << "G"
          << ", "
          << "H"
          << ", "
          << "I"
          << " ]" << std::endl;
      const std::string subMat[9] = {"A", "B", "C", "D", "E", "F", "G", "H", "I"};
      for(unsigned int i = 0; i < 3; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
          str << "\t " << subMat[i * 3 + j] << " = "
              << "[ " << data(i, j)(0, 0) << ", " << data(i, j)(0, 1) << ", " << data(i, j)(0, 2) << "; " << data(i, j)(1, 0) << ", " << data(i, j)(1, 1) << ", " << data(i, j)(1, 2) << "; " << data(i, j)(2, 0) << ", " << data(i, j)(2, 1) << ", " << data(i, j)(2, 2) << " ]" << std::endl;
          ;
        }
      }
      *this << str.str();
    }
    return *this;
  }

  template< MessageType T_Type >
  Message< T_Type > &Message< T_Type >::operator<<(const typename MathObject< float, Degree::Degree4 >::Object &data)
  {
    if(static_cast< int >(T_Type) <= common::Options::instance()->verbose || common::Options::instance()->debug) {
      std::ostringstream str;
      str << "[ "
          << "A"
          << ", "
          << "B"
          << ", "
          << "C"
          << "; "
          << "D"
          << ", "
          << "E"
          << ", "
          << "F"
          << "; "
          << "G"
          << ", "
          << "H"
          << ", "
          << "I"
          << " ]" << std::endl;
      const std::string subMat[9] = {"A", "B", "C", "D", "E", "F", "G", "H", "I"};
      for(unsigned int i = 0; i < 3; ++i) {
        for(unsigned int j = 0; j < 3; ++j) {
          str << "\t " << subMat[i * 3 + j] << " = "
              << "[ " << data(i, j)(0, 0) << ", " << data(i, j)(0, 1) << ", " << data(i, j)(0, 2) << "; " << data(i, j)(1, 0) << ", " << data(i, j)(1, 1) << ", " << data(i, j)(1, 2) << "; " << data(i, j)(2, 0) << ", " << data(i, j)(2, 1) << ", " << data(i, j)(2, 2) << " ]" << std::endl;
          ;
        }
      }
      *this << str.str();
    }
    return *this;
  }

  INSTANTIATE_CLASS(Message, 5, TEMPLATE_ARGS(common::MessageType::None, common::MessageType::Info, common::MessageType::Warning, common::MessageType::Error, common::MessageType::Debug))

  INSTANTIATE_CLASS_OPLL(Message< TEMPLATE_CLASS_PARAM_1 > &, , Message, 0, 5, common::MessageType, TEMPLATE_ARGS(common::MessageType::None, common::MessageType::Info, common::MessageType::Warning, common::MessageType::Error, common::MessageType::Debug), 5, common::MessageModificatorType, TEMPLATE_ARGS(common::MessageModificatorType::Endl, common::MessageModificatorType::Fill, common::MessageModificatorType::Fixed, common::MessageModificatorType::Precision, common::MessageModificatorType::Scientific), TEMPLATE_PARAMS(MessageModificator< TEMPLATE_PARAM_1 > &modificator))


} // namespace gmshfem::common


namespace gmshfem::msg
{


  common::MessageModificator< common::MessageModificatorType::Endl > endl;
  common::MessageModificator< common::MessageModificatorType::Flush > flush;
  common::MessageModificator< common::MessageModificatorType::Fill > fill;
  common::MessageModificator< common::MessageModificatorType::Fixed > fixed;
  common::MessageModificator< common::MessageModificatorType::Precision > precision;
  common::MessageModificator< common::MessageModificatorType::Scientific > scientific;

  common::Message< common::MessageType::None > print;
  common::Message< common::MessageType::Info > info;
  common::Message< common::MessageType::Warning > warning;
  common::Message< common::MessageType::Error > error;
  common::Message< common::MessageType::Debug > debug;

  unsigned int indent()
  {
    common::s_indentation++;
    return common::s_indentation;
  }

  unsigned int unindent()
  {
    common::s_indentation--;
    return common::s_indentation;
  }

  void setIndentation(const unsigned int level)
  {
    common::s_indentation = level;
  }


} // namespace gmshfem::msg
