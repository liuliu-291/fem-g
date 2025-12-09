// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MATRIXOPTIONS
#define H_GMSHFEM_MATRIXOPTIONS

#include <string>

namespace gmshfem::system
{

  class MatrixOptions
  {
   private:
    bool _symmetric;
    bool _hermitian;

    void _parseOptions(const std::string &options);
    void _activateOptions(const std::string &options);

   public:
    explicit MatrixOptions(const std::string &options);
    MatrixOptions(const MatrixOptions &other);
    ~MatrixOptions();

    void init(const std::string &options);
    // Return to its state after initialization
    void clean();

    bool symmetric() const;
    bool hermitian() const;
  };


} // namespace gmshfem::system

#endif // H_GMSHFEM_MATRIXOPTIONS
