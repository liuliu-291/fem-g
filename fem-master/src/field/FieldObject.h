// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FIELDOBJECT
#define H_GMSHFEM_FIELDOBJECT

#include "MathObject.h"

namespace gmshfem::field
{
  enum class FunctionSpaceTypeForm0;
  enum class FunctionSpaceTypeForm1;
  enum class FunctionSpaceTypeForm2;
  enum class FunctionSpaceTypeForm3;
} // namespace gmshfem::field

namespace gmshfem::field
{


  enum class Form { // Mathematical differential form
    Form0 = 0, // H^1
    Form1 = 1, // H(curl)
    Form2 = 2, // H(div)
    Form3 = 3 // H^0
  };


  // Degree associate to forms
  template< Form T_Form >
  struct fieldDegreeOfForm {
    constexpr static gmshfem::Degree value = gmshfem::Degree::Degree0;
  };

  template<>
  struct fieldDegreeOfForm< Form::Form1 > {
    constexpr static gmshfem::Degree value = gmshfem::Degree::Degree1;
  };

  template<>
  struct fieldDegreeOfForm< Form::Form2 > {
    constexpr static gmshfem::Degree value = gmshfem::Degree::Degree1;
  };

  template< Form T_Form >
  using DegreeOfForm = fieldDegreeOfForm< T_Form >;

  // Degree associate to compound forms
  template< Form T_Form >
  struct fieldDegreeOfCompoundForm {
    constexpr static gmshfem::Degree value = gmshfem::Degree::Degree1;
  };

  template<>
  struct fieldDegreeOfCompoundForm< Form::Form1 > {
    constexpr static gmshfem::Degree value = gmshfem::Degree::Degree2;
  };

  template<>
  struct fieldDegreeOfCompoundForm< Form::Form2 > {
    constexpr static gmshfem::Degree value = gmshfem::Degree::Degree2;
  };

  template< Form T_Form >
  using DegreeOfCompoundForm = fieldDegreeOfCompoundForm< T_Form >;


  // Resulting Form of the exterior derivative operator
  template< Form T_Form >
  struct fieldNextForm {
    constexpr static gmshfem::field::Form value = gmshfem::field::Form::Form0;
  };

  template<>
  struct fieldNextForm< Form::Form0 > {
    constexpr static gmshfem::field::Form value = gmshfem::field::Form::Form1;
  };

  template<>
  struct fieldNextForm< Form::Form1 > {
    constexpr static gmshfem::field::Form value = gmshfem::field::Form::Form2;
  };

  template<>
  struct fieldNextForm< Form::Form2 > {
    constexpr static gmshfem::field::Form value = gmshfem::field::Form::Form3;
  };

  template< Form T_Form >
  using NextForm = fieldNextForm< T_Form >;


  // Resulting Form of the exterior derivative operator
  template< Form T_Form >
  struct fieldPastForm {
  };

  template<>
  struct fieldPastForm< Form::Form0 > {
    constexpr static gmshfem::field::Form value = gmshfem::field::Form::Form2; // set such that d(Form2) has the same degree as Form0
  };

  template<>
  struct fieldPastForm< Form::Form1 > {
    constexpr static gmshfem::field::Form value = gmshfem::field::Form::Form0;
  };

  template<>
  struct fieldPastForm< Form::Form2 > {
    constexpr static gmshfem::field::Form value = gmshfem::field::Form::Form1;
  };

  template<>
  struct fieldPastForm< Form::Form3 > {
    constexpr static gmshfem::field::Form value = gmshfem::field::Form::Form2;
  };

  template< Form T_Form >
  using PastForm = fieldPastForm< T_Form >;


  // Name of the form
  template< Form T_Form >
  struct nameOfFrom {
  };

  template<>
  struct nameOfFrom< Form::Form0 > {
    constexpr static const char *value = "0-form";
  };

  template<>
  struct nameOfFrom< Form::Form1 > {
    constexpr static const char *value = "1-form";
  };

  template<>
  struct nameOfFrom< Form::Form2 > {
    constexpr static const char *value = "2-form";
  };

  template<>
  struct nameOfFrom< Form::Form3 > {
    constexpr static const char *value = "3-form";
  };

  template< Form T_Form >
  using NameOfForm = nameOfFrom< T_Form >;


  // Number of components associate to forms
  template< Form T_Form >
  struct fieldNumberOfComponentsOfForm {
    constexpr static unsigned int value = 1;
  };

  template<>
  struct fieldNumberOfComponentsOfForm< Form::Form1 > {
    constexpr static unsigned int value = 3;
  };

  template<>
  struct fieldNumberOfComponentsOfForm< Form::Form2 > {
    constexpr static unsigned int value = 3;
  };

  template< Form T_Form >
  using NumberOfComponentsOfForm = fieldNumberOfComponentsOfForm< T_Form >;


  // Function space of form
  template< Form T_Form >
  struct fieldFunctionSpaceOfForm {
  };

  template<>
  struct fieldFunctionSpaceOfForm< Form::Form0 > {
    typedef gmshfem::field::FunctionSpaceTypeForm0 T_FunctionSpace;
  };

  template<>
  struct fieldFunctionSpaceOfForm< Form::Form1 > {
    typedef gmshfem::field::FunctionSpaceTypeForm1 T_FunctionSpace;
  };

  template<>
  struct fieldFunctionSpaceOfForm< Form::Form2 > {
    typedef gmshfem::field::FunctionSpaceTypeForm2 T_FunctionSpace;
  };

  template<>
  struct fieldFunctionSpaceOfForm< Form::Form3 > {
    typedef gmshfem::field::FunctionSpaceTypeForm3 T_FunctionSpace;
  };

  template< Form T_Form >
  using FunctionSpaceOfForm = typename fieldFunctionSpaceOfForm< T_Form >::T_FunctionSpace;


} // namespace gmshfem::field

namespace gmshfem
{
  typedef field::Form form;
}


#endif // H_GMSHFEM_FIELDOBJECT
