// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "CSVio.h"

#include "Exception.h"
#include "Message.h"

#include <iostream>
#include <sstream>

namespace gmshfem::common
{


  //
  // class CSVModificator
  //

  CSVModificator::CSVModificator(const CSVModificator::Type type) :
    _type(type), _integer(0)
  {
  }

  CSVModificator::~CSVModificator()
  {
  }

  CSVModificator::Type CSVModificator::type() const
  {
    return _type;
  }

  int CSVModificator::integer() const
  {
    return _integer;
  }

  void CSVModificator::reset()
  {
    _integer = 0;
  }

  CSVModificator &CSVModificator::operator()(unsigned int integer)
  {
    _integer = integer;
    return *this;
  }

  //
  // class CSVio
  //

  void CSVio::_checks()
  {
    _fluxChecks();
    _firstLineChecks();
  }

  void CSVio::_fluxChecks() const
  {
    if(!_file.is_open()) {
      throw common::Exception("CSVio: no file is opened");
    }
  }

  void CSVio::_firstLineChecks()
  {
    if(_firstColumn) {
      _firstColumn = false;
    }
    else {
      _file << _separator;
    }
  }

  void CSVio::_readNextLine()
  {
    if(isEOF()) {
      msg::warning << "CSVio: end of file is reached" << msg::endl;
    }

    std::string line;

    bool isInsideQuotes = false;
    std::string current;

    do {
      std::getline(_file, line);
      for(auto i = 0ULL; i < line.size(); ++i) {
        if(line[i] == '\"') {
          if(isInsideQuotes) {
            if(line[i + 1] == '\"') {
              current.push_back(line[i]);
              i++;
            }
            else {
              isInsideQuotes = false;
            }
          }
          else {
            isInsideQuotes = true;
          }
        }
        else if(line[i] == _separator) {
          if(isInsideQuotes) {
            current.push_back(line[i]);
          }
          else {
            _line.push(current);
            current.clear();
          }
        }
        else if(line[i] == '\n') {
          if(isInsideQuotes) {
            current.push_back(line[i]);
          }
        }
        else if(line[i] == '#') {
          if(isInsideQuotes) {
            current.push_back(line[i]);
          }
          else {
            break;
          }
        }
        else {
          current.push_back(line[i]);
        }
      }
      if(!isInsideQuotes) {
        // The last line is parse
        if(current.size() != 0) {
          _line.push(current);
          current.clear();
        }
      }
    } while(isInsideQuotes);

    if(_line.size() == 0 && !isEOF()) {
      _readNextLine();
    }
  }

  CSVio::CSVio(const std::string &path, const char separator, const OpeningMode &opMode) :
    _separator(separator), _file(), _firstColumn(true), _line()
  {
    open(path, opMode);
    if(!isOpen()) {
      throw common::Exception("Unable to open CSV file " + path);
    }
  }

  CSVio::~CSVio()
  {
    close();
  }

  void CSVio::separator(const char separator)
  {
    _separator = separator;
  }

  char CSVio::separator() const
  {
    return _separator;
  }

  void CSVio::writeComment(const std::string &comment)
  {
    _checks();
    _file << "# " << comment << std::endl;
  }

  bool CSVio::open(const std::string &path, const OpeningMode &opMode)
  {
    close();
    try {
      if(opMode == OpeningMode::Append) {
        _file.open(path + ".csv", std::fstream::app);
      }
      else if(opMode == OpeningMode::Reading) {
        _file.open(path + ".csv", std::fstream::in);
      }
      else if(opMode == OpeningMode::NewFile) {
        _file.open(path + ".csv", std::fstream::out);
      }
    }
    catch(...) {
      return false;
    }

    return true;
  }

  bool CSVio::isOpen() const
  {
    return _file.is_open();
  }

  bool CSVio::isEOF() const
  {
    return _file.eof();
  }

  bool CSVio::isEOL() const
  {
    return _line.empty();
  }

  void CSVio::close()
  {
    if(_file.is_open()) {
      _file.close();
    }
  }

  CSVio &CSVio::operator<<(CSVModificator &modificator)
  {
    if(!_file.is_open()) {
      throw common::Exception("CSVio: no file is opened");
    }

    if(modificator.type() == CSVModificator::Type::Endl) {
      _firstColumn = true;
      _file << std::endl;
      modificator.reset();
    }
    else if(modificator.type() == CSVModificator::Type::Fixed) {
      _file << std::fixed;
      modificator.reset();
    }
    else if(modificator.type() == CSVModificator::Type::Precision) {
      _file.precision(modificator.integer());
      modificator.reset();
    }
    else if(modificator.type() == CSVModificator::Type::Scientific) {
      _file << std::scientific;
      modificator.reset();
    }

    return *this;
  }

  CSVio &CSVio::operator<<(const char &data)
  {
    _checks();
    _file << "\"" << data << "\"";
    return *this;
  }

  CSVio &CSVio::operator<<(const std::string &data)
  {
    _checks();
    _file << "\"" << data << "\"";
    return *this;
  }

  CSVio &CSVio::operator<<(const short &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const unsigned short &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const int &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const unsigned int &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const long &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const unsigned long &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const long long &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const unsigned long long &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const float &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const double &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const long double &data)
  {
    _checks();
    _file << data;
    return *this;
  }

  CSVio &CSVio::operator<<(const std::complex< float > &data)
  {
    _checks();
    _file << data.real() << _separator << data.imag();
    return *this;
  }

  CSVio &CSVio::operator<<(const std::complex< double > &data)
  {
    _checks();
    _file << data.real() << _separator << data.imag();
    return *this;
  }

  CSVio &CSVio::operator<<(const std::complex< long double > &data)
  {
    _checks();
    _file << data.real() << _separator << data.imag();
    return *this;
  }

  CSVio &CSVio::operator<<(const common::Timer &data)
  {
    _checks();
    _file << data.time();
    return *this;
  }

  CSVio &CSVio::operator<<(const Eigen::Vector3< std::complex< double > > &data)
  {
    return (*this << data(0) << data(1) << data(2));
  }

  CSVio &CSVio::operator<<(const Eigen::Vector3< double > &data)
  {
    return (*this << data(0) << data(1) << data(2));
  }

  CSVio &CSVio::operator<<(const Eigen::Vector3< std::complex< float > > &data)
  {
    return (*this << data(0) << data(1) << data(2));
  }

  CSVio &CSVio::operator<<(const Eigen::Vector3< float > &data)
  {
    return (*this << data(0) << data(1) << data(2));
  }

  CSVio &CSVio::operator<<(const Eigen::Matrix3< std::complex< double > > &data)
  {
    return (*this << data(0, 0) << data(0, 1) << data(0, 2)
                  << data(1, 0) << data(1, 1) << data(1, 2)
                  << data(2, 0) << data(2, 1) << data(2, 2));
  }

  CSVio &CSVio::operator<<(const Eigen::Matrix3< double > &data)
  {
    _checks();
    return (*this << data(0, 0) << data(0, 1) << data(0, 2)
                  << data(1, 0) << data(1, 1) << data(1, 2)
                  << data(2, 0) << data(2, 1) << data(2, 2));
  }

  CSVio &CSVio::operator<<(const Eigen::Matrix3< std::complex< float > > &data)
  {
    return (*this << data(0, 0) << data(0, 1) << data(0, 2)
                  << data(1, 0) << data(1, 1) << data(1, 2)
                  << data(2, 0) << data(2, 1) << data(2, 2));
  }

  CSVio &CSVio::operator<<(const Eigen::Matrix3< float > &data)
  {
    return (*this << data(0, 0) << data(0, 1) << data(0, 2)
                  << data(1, 0) << data(1, 1) << data(1, 2)
                  << data(2, 0) << data(2, 1) << data(2, 2));
  }

  CSVio &CSVio::operator>>(char &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      std::string str = _line.front();
      if(str.size() == 3) {
        if(str[0] == '\"' && str[2] == '\"') {
          data = str[1];
        }
        else {
          throw common::Exception("CSVio: unable to interpret \"" + str + "\" as a 'char'");
        }
      }
      else if(str.size() == 1) {
        data = str[0];
      }
      else {
        throw common::Exception("CSVio: unable to interpret \"" + str + "\" as a 'char'");
      }
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(std::string &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = _line.front();
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(short &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      int integer = stoi(_line.front());
      data = static_cast< short >(integer);
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(unsigned short &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      unsigned long integer = stoul(_line.front());
      data = static_cast< unsigned short >(integer);
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(int &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = stoi(_line.front());
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(unsigned int &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = stoi(_line.front());
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(long &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = stol(_line.front());
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(unsigned long &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = stoul(_line.front());
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(long long &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = stoll(_line.front());
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(unsigned long long &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = stoull(_line.front());
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(float &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = stof(_line.front());
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(double &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = stod(_line.front());
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(long double &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() != 0) {
      data = stold(_line.front());
      _line.pop();
    }
    return *this;
  }

  CSVio &CSVio::operator>>(std::complex< float > &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 2) {
      data = stof(_line.front());
      _line.pop();
      data += std::complex< float >(0.f, 1.f) * stof(_line.front());
      _line.pop();
    }
    else {
      msg::warning << "Not enouth data to fill a 'std::complex< float >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(std::complex< double > &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 2) {
      data = stod(_line.front());
      _line.pop();
      data += std::complex< double >(0., 1.) * stod(_line.front());
      _line.pop();
    }
    else {
      msg::warning << "Not enouth data to fill a 'std::complex< double >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(std::complex< long double > &data)
  {
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 2) {
      data = stold(_line.front());
      _line.pop();
      data += std::complex< long double >(0.L, 1.L) * stold(_line.front());
      _line.pop();
    }
    else {
      msg::warning << "Not enouth data to fill a 'std::complex< long double >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(Eigen::Vector3< std::complex< double > > &data)
  {
    std::complex< double > comp;
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 6) {
      for(auto i = 0; i < 3; ++i) {
        comp = stod(_line.front());
        _line.pop();
        comp += std::complex< double >(0., 1.) * stod(_line.front());
        _line.pop();
        data(i) = comp;
      }
    }
    else {
      msg::warning << "Not enouth data to fill a 'Vector3< std::complex< double > >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(Eigen::Vector3< double > &data)
  {
    double dbl;
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 3) {
      for(auto i = 0; i < 3; ++i) {
        dbl = stod(_line.front());
        _line.pop();
        data(i) = dbl;
      }
    }
    else {
      msg::warning << "Not enouth data to fill a 'Vector3< double >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(Eigen::Vector3< std::complex< float > > &data)
  {
    std::complex< float > comp;
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 6) {
      for(auto i = 0; i < 3; ++i) {
        comp = stof(_line.front());
        _line.pop();
        comp += std::complex< float >(0.f, 1.f) * stof(_line.front());
        _line.pop();
        data(i) = comp;
      }
    }
    else {
      msg::warning << "Not enouth data to fill a 'Vector3< std::complex< float > >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(Eigen::Vector3< float > &data)
  {
    float dbl;
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 3) {
      for(auto i = 0; i < 3; ++i) {
        dbl = stof(_line.front());
        _line.pop();
        data(i) = dbl;
      }
    }
    else {
      msg::warning << "Not enouth data to fill a 'Vector3< float >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(Eigen::Matrix3< std::complex< double > > &data)
  {
    std::complex< double > comp;
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 18) {
      for(auto i = 0; i < 3; ++i) {
        for(auto j = 0; j < 3; ++j) {
          comp = stod(_line.front());
          _line.pop();
          comp += std::complex< double >(0., 1.) * stod(_line.front());
          _line.pop();
          data(i, j) = comp;
        }
      }
    }
    else {
      msg::warning << "Not enouth data to fill a 'Matrix3< std::complex< double > >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(Eigen::Matrix3< double > &data)
  {
    double dbl;
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 9) {
      for(auto i = 0; i < 3; ++i) {
        for(auto j = 0; j < 3; ++j) {
          dbl = stod(_line.front());
          _line.pop();
          data(i, j) = dbl;
        }
      }
    }
    else {
      msg::warning << "Not enouth data to fill a 'Matrix3< double >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(Eigen::Matrix3< std::complex< float > > &data)
  {
    std::complex< float > comp;
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 18) {
      for(auto i = 0; i < 3; ++i) {
        for(auto j = 0; j < 3; ++j) {
          comp = stof(_line.front());
          _line.pop();
          comp += std::complex< float >(0.f, 1.f) * stof(_line.front());
          _line.pop();
          data(i, j) = comp;
        }
      }
    }
    else {
      msg::warning << "Not enouth data to fill a 'Matrix3< std::complex< float > >'" << msg::endl;
    }
    return *this;
  }

  CSVio &CSVio::operator>>(Eigen::Matrix3< float > &data)
  {
    float dbl;
    if(_line.size() == 0) _readNextLine();
    if(_line.size() >= 9) {
      for(auto i = 0; i < 3; ++i) {
        for(auto j = 0; j < 3; ++j) {
          dbl = stof(_line.front());
          _line.pop();
          data(i, j) = dbl;
        }
      }
    }
    else {
      msg::warning << "Not enouth data to fill a 'Matrix3< float >'" << msg::endl;
    }
    return *this;
  }


} // namespace gmshfem::common


namespace gmshfem::csv
{


  common::CSVModificator endl(common::CSVModificator::Type::Endl);
  common::CSVModificator fixed(common::CSVModificator::Type::Fixed);
  common::CSVModificator precision(common::CSVModificator::Type::Precision);
  common::CSVModificator scientific(common::CSVModificator::Type::Scientific);


} // namespace gmshfem::csv
