// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "MatrixOptions.h"

#include "Message.h"
#include "Options.h"

namespace gmshfem::system
{


  void MatrixOptions::_activateOptions(const std::string &options)
  {
    if(options == "Symmetric") {
      _symmetric = true;
    }
    else if(options == "Hermitian") {
      _hermitian = true;
    }
  }

  void MatrixOptions::_parseOptions(const std::string &options)
  {
    std::string currentOption;
    for(auto i = 0ULL; i < options.size(); ++i) {
      if(options[i] == ';') {
        _activateOptions(options);
        currentOption.clear();
      }
      else {
        currentOption += options[i];
      }
    }
    _activateOptions(options);
  }

  MatrixOptions::MatrixOptions(const std::string &options)
  {
    init(options);
  }

  MatrixOptions::MatrixOptions(const MatrixOptions &other) :
    _symmetric(other._symmetric), _hermitian(other._hermitian)
  {
  }

  MatrixOptions::~MatrixOptions()
  {
  }

  void MatrixOptions::init(const std::string &options)
  {
    _symmetric = common::Options::instance()->forceSymmetric;
    _hermitian = common::Options::instance()->forceHermitian;
    _parseOptions(options);
  }

  void MatrixOptions::clean()
  {
    _symmetric = common::Options::instance()->forceSymmetric;
    _hermitian = common::Options::instance()->forceHermitian;
  }

  bool MatrixOptions::symmetric() const
  {
    return _symmetric;
  }

  bool MatrixOptions::hermitian() const
  {
    return _hermitian;
  }


} // namespace gmshfem::system
