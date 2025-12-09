// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_CSVIO
#define H_GMSHFEM_CSVIO

#include "MathObject.h"
#include "Timer.h"
#include "io.h"

#include <complex>
#include <fstream>
#include <initializer_list>
#include <queue>
#include <string>

namespace gmshfem::common
{


  class CSVModificator
  {
   public:
    enum class Type {
      Endl,
      Fixed,
      Precision,
      Scientific
    };

   private:
    const Type _type;
    int _integer;

   public:
    explicit CSVModificator(const CSVModificator::Type type);
    ~CSVModificator();

    CSVModificator::Type type() const;
    int integer() const;

    void reset();

    CSVModificator &operator()(unsigned int integer);
  };

  class CSVio
  {
   private:
    char _separator;
    std::fstream _file;
    bool _firstColumn;
    std::queue< std::string > _line;

    void _checks();
    void _fluxChecks() const;
    void _firstLineChecks();
    void _readNextLine();

   public:
    CSVio();
    CSVio(const std::string &path, const char separator = ';', const OpeningMode &opMode = OpeningMode::NewFile);
    ~CSVio();

    void separator(const char separator);
    char separator() const;

    void writeComment(const std::string &comment);

    bool open(const std::string &path, const OpeningMode &opMode = OpeningMode::NewFile);
    bool isOpen() const;
    bool isEOF() const;
    bool isEOL() const;
    void close();

    CSVio &operator<<(CSVModificator &modificator);

    // writing

    CSVio &operator<<(const char &data);
    CSVio &operator<<(const std::string &data);

    CSVio &operator<<(const short &data);
    CSVio &operator<<(const unsigned short &data);
    CSVio &operator<<(const int &data);
    CSVio &operator<<(const unsigned int &data);
    CSVio &operator<<(const long &data);
    CSVio &operator<<(const unsigned long &data);
    CSVio &operator<<(const long long &data);
    CSVio &operator<<(const unsigned long long &data);

    CSVio &operator<<(const float &data);
    CSVio &operator<<(const double &data);
    CSVio &operator<<(const long double &data);

    CSVio &operator<<(const std::complex< float > &data);
    CSVio &operator<<(const std::complex< double > &data);
    CSVio &operator<<(const std::complex< long double > &data);

    CSVio &operator<<(const common::Timer &data);

    CSVio &operator<<(const Eigen::Vector3< std::complex< double > > &data);
    CSVio &operator<<(const Eigen::Vector3< double > &data);
    CSVio &operator<<(const Eigen::Vector3< std::complex< float > > &data);
    CSVio &operator<<(const Eigen::Vector3< float > &data);

    CSVio &operator<<(const Eigen::Matrix3< std::complex< double > > &data);
    CSVio &operator<<(const Eigen::Matrix3< double > &data);
    CSVio &operator<<(const Eigen::Matrix3< std::complex< float > > &data);
    CSVio &operator<<(const Eigen::Matrix3< float > &data);

    // reading

    CSVio &operator>>(char &data);
    CSVio &operator>>(std::string &data);

    CSVio &operator>>(short &data);
    CSVio &operator>>(unsigned short &data);
    CSVio &operator>>(int &data);
    CSVio &operator>>(unsigned int &data);
    CSVio &operator>>(long &data);
    CSVio &operator>>(unsigned long &data);
    CSVio &operator>>(long long &data);
    CSVio &operator>>(unsigned long long &data);

    CSVio &operator>>(float &data);
    CSVio &operator>>(double &data);
    CSVio &operator>>(long double &data);

    CSVio &operator>>(std::complex< float > &data);
    CSVio &operator>>(std::complex< double > &data);
    CSVio &operator>>(std::complex< long double > &data);

    CSVio &operator>>(Eigen::Vector3< std::complex< double > > &data);
    CSVio &operator>>(Eigen::Vector3< double > &data);
    CSVio &operator>>(Eigen::Vector3< std::complex< float > > &data);
    CSVio &operator>>(Eigen::Vector3< float > &data);

    CSVio &operator>>(Eigen::Matrix3< std::complex< double > > &data);
    CSVio &operator>>(Eigen::Matrix3< double > &data);
    CSVio &operator>>(Eigen::Matrix3< std::complex< float > > &data);
    CSVio &operator>>(Eigen::Matrix3< float > &data);
  };


} // namespace gmshfem::common

namespace gmshfem::csv
{


  extern common::CSVModificator endl;
  extern common::CSVModificator fixed;
  extern common::CSVModificator precision;
  extern common::CSVModificator scientific;


} // namespace gmshfem::csv


#endif // H_GMSHFEM_CSVIO
