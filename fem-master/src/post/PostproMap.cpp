// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "PostproMap.h"

#include "FieldInterface.h"
#include "OmpInterface.h"
#include "Options.h"
#include "instantiate.h"

#include <gmsh.h>

namespace gmshfem::post
{


  template< class T_Scalar, Degree T_Degree >
  PostproMap< T_Scalar, T_Degree >::PostproMap(const std::string &name, const int tag) :
    _name(name), _tag(tag)
  {
    _tag = gmsh::view::add(_name, _tag);
  }

  template< class T_Scalar, Degree T_Degree >
  PostproMap< T_Scalar, T_Degree >::~PostproMap()
  {
    gmsh::view::remove(_tag);
  }

  template< class T_Scalar, Degree T_Degree >
  static void s_addModelData(const int tag, const int step, const double time, const int partition, const unsigned int nbrNodesByElements, std::vector< std::size_t > &gmshElementsTags, std::vector< std::vector< double > > &gmshData, std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > &values)
  {
    if(scalar::IsComplex< T_Scalar >::value) {
      // real part
#pragma omp for
      for(auto i = 0ULL; i < gmshElementsTags.size(); ++i) {
        for(auto j = 0U; j < nbrNodesByElements; ++j) {
          gmshData[i].push_back(std::real(values[i * nbrNodesByElements + j]));
        }
      }
#pragma omp single
      gmsh::view::addModelData(tag, 2 * step, "", "ElementNodeData", gmshElementsTags, gmshData, time, 1, partition);

      // imaginary part
#pragma omp for
      for(auto i = 0ULL; i < gmshElementsTags.size(); ++i) {
        gmshData[i].clear();
        for(auto j = 0U; j < nbrNodesByElements; ++j) {
          gmshData[i].push_back(std::imag(values[i * nbrNodesByElements + j]));
        }
      }
#pragma omp single
      gmsh::view::addModelData(tag, 2 * step + 1, "", "ElementNodeData", gmshElementsTags, gmshData, time, 1, partition);
    }
    else {
#pragma omp for
      for(auto i = 0ULL; i < gmshElementsTags.size(); ++i) {
        for(auto j = 0U; j < nbrNodesByElements; ++j) {
          gmshData[i].push_back(std::real(values[i * nbrNodesByElements + j]));
        }
      }
#pragma omp single
      if(gmshData.size() != 0) {
        gmsh::view::addModelData(tag, step, "", "ElementNodeData", gmshElementsTags, gmshData, time, 1, partition);
      }
    }
  }

  template< class T_Scalar, Degree T_Degree >
  static void s_addModelData(const int tag, const int step, const double time, const int partition, const unsigned int nbrNodesByElements, std::vector< std::size_t > &gmshElementsTags, std::vector< std::vector< double > > &gmshData, std::vector< typename MathObject< T_Scalar, Degree::Degree1 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree1 >::Object > > &values)
  {
    if(scalar::IsComplex< T_Scalar >::value) {
      // real part
#pragma omp for
      for(auto i = 0ULL; i < gmshElementsTags.size(); ++i) {
        gmshData[i].reserve(3);
        for(auto j = 0U; j < nbrNodesByElements; ++j) {
          for(auto k = 0; k < 3; ++k) {
            gmshData[i].push_back(std::real(values[i * nbrNodesByElements + j](k)));
          }
        }
      }
#pragma omp single
      gmsh::view::addModelData(tag, 2 * step, "", "ElementNodeData", gmshElementsTags, gmshData, time, 3, partition);

      // imaginary part
#pragma omp for
      for(auto i = 0ULL; i < gmshElementsTags.size(); ++i) {
        gmshData[i].clear();
        gmshData[i].reserve(3);
        for(auto j = 0U; j < nbrNodesByElements; ++j) {
          for(auto k = 0; k < 3; ++k) {
            gmshData[i].push_back(std::imag(values[i * nbrNodesByElements + j](k)));
          }
        }
      }
#pragma omp single
      gmsh::view::addModelData(tag, 2 * step + 1, "", "ElementNodeData", gmshElementsTags, gmshData, time, 3, partition);
    }
    else {
#pragma omp for
      for(auto i = 0ULL; i < gmshElementsTags.size(); ++i) {
        gmshData[i].reserve(3);
        for(auto j = 0U; j < nbrNodesByElements; ++j) {
          for(auto k = 0; k < 3; ++k) {
            gmshData[i].push_back(std::real(values[i * nbrNodesByElements + j](k)));
          }
        }
      }
#pragma omp single
      if(gmshData.size() != 0) {
        gmsh::view::addModelData(tag, step, "", "ElementNodeData", gmshElementsTags, gmshData, time, 3, partition);
      }
    }
  }

  template< class T_Scalar, Degree T_Degree >
  static void s_addModelData(const int tag, const int step, const double time, const int partition, const unsigned int nbrNodesByElements, std::vector< std::size_t > &gmshElementsTags, std::vector< std::vector< double > > &gmshData, std::vector< typename MathObject< T_Scalar, Degree::Degree2 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree2 >::Object > > &values)
  {
    if(scalar::IsComplex< T_Scalar >::value) {
      // real part
#pragma omp for
      for(auto i = 0ULL; i < gmshElementsTags.size(); ++i) {
        gmshData[i].reserve(9);
        for(auto j = 0U; j < nbrNodesByElements; ++j) {
          for(auto k = 0; k < 3; ++k) {
            for(auto l = 0; l < 3; ++l) {
              gmshData[i].push_back(std::real(values[i * nbrNodesByElements + j](k, l)));
            }
          }
        }
      }
#pragma omp single
      gmsh::view::addModelData(tag, 2 * step, "", "ElementNodeData", gmshElementsTags, gmshData, time, 9, partition);

      // imaginary part
#pragma omp for
      for(auto i = 0ULL; i < gmshElementsTags.size(); ++i) {
        gmshData[i].clear();
        gmshData[i].reserve(9);
        for(auto j = 0U; j < nbrNodesByElements; ++j) {
          for(auto k = 0; k < 3; ++k) {
            for(auto l = 0; l < 3; ++l) {
              gmshData[i].push_back(std::imag(values[i * nbrNodesByElements + j](k, l)));
            }
          }
        }
      }
#pragma omp single
      gmsh::view::addModelData(tag, 2 * step + 1, "", "ElementNodeData", gmshElementsTags, gmshData, time, 9, partition);
    }
    else {
#pragma omp for
      for(auto i = 0ULL; i < gmshElementsTags.size(); ++i) {
        gmshData[i].reserve(9);
        for(auto j = 0U; j < nbrNodesByElements; ++j) {
          for(auto k = 0; k < 3; ++k) {
            for(auto l = 0; l < 3; ++l) {
              gmshData[i].push_back(std::real(values[i * nbrNodesByElements + j](k, l)));
            }
          }
        }
      }
#pragma omp single
      if(gmshData.size() != 0) {
        gmsh::view::addModelData(tag, step, "", "ElementNodeData", gmshElementsTags, gmshData, time, 9, partition);
      }
    }
  }

  template< class T_Scalar, Degree T_Degree >
  void PostproMap< T_Scalar, T_Degree >::append(const function::Function< T_Scalar, T_Degree > &function, const domain::GeometricObject &domain, const int step, const double time, const int partition)
  {
    for(auto itEntity = domain.cbegin(); itEntity != domain.cend(); ++itEntity) {
      std::vector< int > elementTypes;
      gmsh::model::mesh::getElementTypes(elementTypes, itEntity->first, itEntity->second);

      for(auto typeIndex = 0ULL; typeIndex < elementTypes.size(); ++typeIndex) {
        std::vector< double > gmshNodesCoord;
        std::string name;
        int dim;
        int degree;
        int gmshNbrNodesByElements;
        int numPrimaryNodes;
        gmsh::model::mesh::getElementProperties(elementTypes[typeIndex], name, dim, degree, gmshNbrNodesByElements, gmshNodesCoord, numPrimaryNodes);
        unsigned int nbrNodesByElements = gmshNbrNodesByElements;
        std::vector< scalar::Precision< T_Scalar > > nodesCoord(nbrNodesByElements * 3, 0.);
        for(auto i = 0U; i < nbrNodesByElements; ++i) {
          for(auto j = 0; j < dim; ++j) {
            nodesCoord[3 * i + j] = gmshNodesCoord[i * dim + j];
          }
        }

        std::vector< typename MathObject< T_Scalar, T_Degree >::Object, numa::allocator< typename MathObject< T_Scalar, T_Degree >::Object > > values;

        std::vector< std::size_t > gmshNodesTags;
        std::vector< double > gmshNodesParametricCoord;
        std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > coord;
        gmsh::model::mesh::getNodesByElementType(elementTypes[typeIndex], gmshNodesTags, gmshNodesCoord, gmshNodesParametricCoord, itEntity->second, false);
        numa::copy(coord, gmshNodesCoord);

        std::vector< std::size_t > gmshElementsTags;
        gmsh::model::mesh::preallocateElementsByType(elementTypes[typeIndex], true, false, gmshElementsTags, gmshNodesTags, itEntity->second);

        std::vector< std::vector< double > > gmshData;
#pragma omp parallel num_threads(omp::getMaxThreads())
        {
          const unsigned int numThreads = omp::getNumThreads();
          const unsigned int myThreadID = omp::getThreadNum();

          gmsh::model::mesh::getElementsByType(elementTypes[typeIndex], gmshElementsTags, gmshNodesTags, itEntity->second, myThreadID, numThreads);
          function.evaluate(values, coord, nodesCoord, elementTypes[typeIndex], *itEntity);
#pragma omp single
          gmshData.resize(gmshElementsTags.size());
          s_addModelData< T_Scalar, T_Degree >(_tag, step, time, partition, nbrNodesByElements, gmshElementsTags, gmshData, values);
        }
      }
    }
  }

  template< class T_Scalar, Degree T_Degree >
  template< field::Form T_Form, class >
  void PostproMap< T_Scalar, T_Degree >::append(const field::Field< T_Scalar, T_Form > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition)
  {
    append(field.getEvaluableFunction(), domain, step, time, partition);
  }


  template< class T_Scalar, Degree T_Degree >
  template< field::Form T_Form, unsigned int T_NumFields, class >
  void PostproMap< T_Scalar, T_Degree >::append(const field::CompoundField< T_Scalar, T_Form, T_NumFields > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition)
  {
    append(field.getEvaluableFunction(), domain, step, time, partition);
  }

  template< class T_Scalar, Degree T_Degree >
  void PostproMap< T_Scalar, T_Degree >::write(const std::string &format, const std::string &path) const
  {
    gmsh::view::write(_tag, path + _name + "." + format);
  }


  INSTANTIATE_CLASS_2(PostproMap, 4, 3, TEMPLATE_ARGS(std::complex< double >, std::complex< float >, double, float), TEMPLATE_ARGS(Degree::Degree0, Degree::Degree1, Degree::Degree2))


  template void PostproMap< std::complex< double >, Degree::Degree0 >::append(const field::Field< std::complex< double >, field::Form::Form0 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree0 >::append(const field::Field< std::complex< float >, field::Form::Form0 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree0 >::append(const field::Field< double, field::Form::Form0 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree0 >::append(const field::Field< float, field::Form::Form0 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree1 >::append(const field::Field< std::complex< double >, field::Form::Form1 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree1 >::append(const field::Field< std::complex< float >, field::Form::Form1 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree1 >::append(const field::Field< double, field::Form::Form1 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree1 >::append(const field::Field< float, field::Form::Form1 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree1 >::append(const field::Field< std::complex< double >, field::Form::Form2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree1 >::append(const field::Field< std::complex< float >, field::Form::Form2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree1 >::append(const field::Field< double, field::Form::Form2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree1 >::append(const field::Field< float, field::Form::Form2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree0 >::append(const field::Field< std::complex< double >, field::Form::Form3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree0 >::append(const field::Field< std::complex< float >, field::Form::Form3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree0 >::append(const field::Field< double, field::Form::Form3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree0 >::append(const field::Field< float, field::Form::Form3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree1 >::append(const field::CompoundField< std::complex< double >, field::Form::Form0, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree1 >::append(const field::CompoundField< std::complex< float >, field::Form::Form0, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree1 >::append(const field::CompoundField< double, field::Form::Form0, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree1 >::append(const field::CompoundField< float, field::Form::Form0, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree2 >::append(const field::CompoundField< std::complex< double >, field::Form::Form1, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree2 >::append(const field::CompoundField< std::complex< float >, field::Form::Form1, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree2 >::append(const field::CompoundField< double, field::Form::Form1, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree2 >::append(const field::CompoundField< float, field::Form::Form1, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree2 >::append(const field::CompoundField< std::complex< double >, field::Form::Form2, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree2 >::append(const field::CompoundField< std::complex< float >, field::Form::Form2, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree2 >::append(const field::CompoundField< double, field::Form::Form2, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree2 >::append(const field::CompoundField< float, field::Form::Form2, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree1 >::append(const field::CompoundField< std::complex< double >, field::Form::Form3, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree1 >::append(const field::CompoundField< std::complex< float >, field::Form::Form3, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree1 >::append(const field::CompoundField< double, field::Form::Form3, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree1 >::append(const field::CompoundField< float, field::Form::Form3, 2 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree1 >::append(const field::CompoundField< std::complex< double >, field::Form::Form0, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree1 >::append(const field::CompoundField< std::complex< float >, field::Form::Form0, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree1 >::append(const field::CompoundField< double, field::Form::Form0, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree1 >::append(const field::CompoundField< float, field::Form::Form0, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree2 >::append(const field::CompoundField< std::complex< double >, field::Form::Form1, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree2 >::append(const field::CompoundField< std::complex< float >, field::Form::Form1, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree2 >::append(const field::CompoundField< double, field::Form::Form1, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree2 >::append(const field::CompoundField< float, field::Form::Form1, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree2 >::append(const field::CompoundField< std::complex< double >, field::Form::Form2, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree2 >::append(const field::CompoundField< std::complex< float >, field::Form::Form2, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree2 >::append(const field::CompoundField< double, field::Form::Form2, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree2 >::append(const field::CompoundField< float, field::Form::Form2, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


  template void PostproMap< std::complex< double >, Degree::Degree1 >::append(const field::CompoundField< std::complex< double >, field::Form::Form3, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< std::complex< float >, Degree::Degree1 >::append(const field::CompoundField< std::complex< float >, field::Form::Form3, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< double, Degree::Degree1 >::append(const field::CompoundField< double, field::Form::Form3, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);
  template void PostproMap< float, Degree::Degree1 >::append(const field::CompoundField< float, field::Form::Form3, 3 > &field, const domain::GeometricObject &domain, const int step, const double time, const int partition);


} // namespace gmshfem::post
