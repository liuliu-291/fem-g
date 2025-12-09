// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "PPMio.h"

#include "Exception.h"
#include "Message.h"
#include "OmpInterface.h"

namespace gmshfem::common
{


  //
  // class PPMio
  //

  PPMio::PPMio() :
    _file()
  {
  }

  PPMio::PPMio(const std::string &path, const PPMFormat &format, const OpeningMode &opMode) :
    _file(), _format(format), _tmp()
  {
    open(path, opMode);
    if(!isOpen()) {
      throw common::Exception("Unable to open PPM file " + path);
    }
  }

  PPMio::~PPMio()
  {
    close();
  }

  bool PPMio::open(const std::string &path, const OpeningMode &opMode)
  {
    close();

    std::string extension = "";
    if(_format == PPMFormat::PPM) {
      extension = ".ppm";
    }
    else if(_format == PPMFormat::PGM) {
      extension = ".pgm";
    }
    else if(_format == PPMFormat::PBM) {
      extension = ".pbm";
    }

    try {
      if(opMode == OpeningMode::Append) {
        _file.open(path + extension, std::fstream::app | std::fstream::binary);
      }
      else if(opMode == OpeningMode::Reading) {
        _file.open(path + extension, std::fstream::in | std::fstream::binary);
      }
      else if(opMode == OpeningMode::NewFile) {
        _file.open(path + extension, std::fstream::out | std::fstream::binary);
      }
    }
    catch(...) {
      return false;
    }

    return true;
  }

  bool PPMio::isOpen() const
  {
    return _file.is_open();
  }

  void PPMio::close()
  {
    if(_file.is_open()) {
      _file.close();
    }
  }

  void PPMio::writeHeader(const unsigned int sizeX, const unsigned int sizeY)
  {
    if(!isOpen()) {
      msg::error << "Unable to write PPM file: no file is opened" << msg::endl;
      return;
    }

    if(_format == PPMFormat::PPM) {
      _file << "P6" << std::endl;
    }
    else if(_format == PPMFormat::PGM) {
      _file << "P5" << std::endl;
    }
    else if(_format == PPMFormat::PBM) {
      _file << "P4" << std::endl;
    }

    _file << sizeX << " " << sizeY << std::endl;
    if(_format != PPMFormat::PBM) {
      _file << "255" << std::endl;
    }
  }

  void PPMio::writeLine(const std::vector< Color > &map)
  {
    if(!isOpen()) {
      msg::error << "Unable to write PPM file: no file is opened" << msg::endl;
      return;
    }

    if(_format == PPMFormat::PPM) {
      if(_tmp.size() != 3 * map.size()) {
        _tmp.resize(3 * map.size());
      }
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < map.size(); ++i) {
        map[i].extractRGB(&_tmp[3 * i]);
      }
    }
    else if(_format == PPMFormat::PGM) {
      if(_tmp.size() != map.size()) {
        _tmp.resize(map.size());
      }
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < map.size(); ++i) {
        map[i].extractGrey(&_tmp[i]);
      }
    }
    else if(_format == PPMFormat::PBM) {
      if(_tmp.size() != (map.size() + 7) / 8) {
        _tmp.resize((map.size() + 7) / 8);
      }
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < _tmp.size(); ++i) {
        _tmp[i] = 0;
      }
#pragma omp parallel for num_threads(omp::getMaxThreads()) schedule(static, 512)
      for(auto i = 0ULL; i < map.size(); ++i) {
        map[i].extractBW(&_tmp[i / 8], i % 8);
      }
    }

    _file.write(reinterpret_cast< const char * >(&_tmp[0]), _tmp.size());
  }


} // namespace gmshfem::common
