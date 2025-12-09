// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "SaveField.h"

#include "FieldInterface.h"
#include "OmpInterface.h"
#include "Options.h"
#include "instantiate.h"

#include <gmsh.h>

namespace gmshfem::post
{


  template< class T_Scalar, field::Form T_Form >
  SaveField< T_Scalar, T_Form >::SaveField(const field::Field< T_Scalar, T_Form > &field, const std::string &format, const std::string &path, const double memorize, const int step, const double time, const int partition) :
    PostInterface(), _field(&field), _format(format), _path(path), _memorize(memorize), _step(step), _time(time), _partition(partition)
  {
  }

  template< class T_Scalar, field::Form T_Form >
  SaveField< T_Scalar, T_Form >::~SaveField()
  {
  }

  template< class T_Scalar >
  static bool s_checkField(const field::Field< T_Scalar, field::Form::Form0 > *const field)
  {
    if(field->getFunctionSpace()->type() == field::FunctionSpaceTypeForm0::Lagrange) {
      return true;
    }
    else if(field->getFunctionSpace()->type() == field::FunctionSpaceTypeForm0::HierarchicalH1 && field->getFunctionSpace()->order() == 1) {
      return true;
    }
    return false;
  }

  template< class T_Scalar, field::Form T_Form >
  int SaveField< T_Scalar, T_Form >::run()
  {
    if(!s_checkField(_field)) {
      msg::warning << "The fast save function cannot be used on " << field::NameOfForm< T_Form >::value << " field using " << _field->getFunctionSpace()->getGmshFemName(false) << " function space, use the general one instead." << msg::endl;
      return -1;
    }
    const int tag = gmsh::view::add(_field->name());
    std::vector< std::size_t > tags;
    std::vector< T_Scalar > data;
    tags.reserve(_field->numberOfDofs());
    data.reserve(_field->numberOfDofs());

    for(auto it = _field->begin(); it != _field->end(); ++it) {
      tags.push_back(it->first->entity());
      if(it->first->type() == dofs::Type::Linked) {
        data.push_back(it->second * _field->getValue(static_cast< const dofs::LinkedDof * >(it->first)->master()));
      }
      else {
        data.push_back(it->second);
      }
    }

    std::vector< std::vector< double > > gmshData(data.size());
    std::string dataType = "NodeData";

    if(scalar::IsComplex< T_Scalar >::value) {
#pragma omp parallel num_threads(omp::getMaxThreads())
      {
        // real part
#pragma omp for
        for(auto i = 0ULL; i < data.size(); ++i) {
          gmshData[i] = {std::real(data[i])};
        }
#pragma omp single
        gmsh::view::addModelData(tag, 2 * _step, "", dataType, tags, gmshData, _time, 1, _partition);

        // imaginary part
#pragma omp for
        for(auto i = 0ULL; i < data.size(); ++i) {
          gmshData[i] = {std::imag(data[i])};
        }
#pragma omp single
        gmsh::view::addModelData(tag, 2 * _step + 1, "", dataType, tags, gmshData, _time, 1, _partition);
      }
    }
    else {
#pragma omp parallel for num_threads(omp::getMaxThreads())
      for(auto i = 0ULL; i < data.size(); ++i) {
        gmshData[i] = {std::real(data[i])};
      }
      if(gmshData.size() != 0) {
        gmsh::view::addModelData(tag, _step, "", dataType, tags, gmshData, _time, 1, _partition);
      }
    }

    if(_format.size()) {
      gmsh::view::write(tag, _path + _field->name() + "." + _format);
    }

    if(!_memorize && !common::Options::instance()->interface) {
      gmsh::view::remove(tag);
      return 0;
    }
    else {
      return tag;
    }
  }

  INSTANTIATE_CLASS_2(SaveField, 4, 1, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(field::Form::Form0))


} // namespace gmshfem::post
