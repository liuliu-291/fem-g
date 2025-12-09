// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_PPMIO
#define H_GMSHFEM_PPMIO

#include "Color.h"
#include "io.h"

#include <fstream>
#include <vector>

namespace gmshfem::common
{


  enum class PPMFormat {
    PPM,
    PGM,
    PBM
  };

  class PPMio
  {
   private:
    std::fstream _file;
    PPMFormat _format;
    std::vector< unsigned char > _tmp;

   public:
    PPMio();
    PPMio(const std::string &path, const PPMFormat &format = PPMFormat::PPM, const OpeningMode &opMode = OpeningMode::NewFile);
    ~PPMio();

    bool open(const std::string &path, const OpeningMode &opMode = OpeningMode::NewFile);
    bool isOpen() const;
    void close();

    void writeHeader(const unsigned int sizeX, const unsigned int sizeY);
    void writeLine(const std::vector< Color > &map);
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_PPMIO
