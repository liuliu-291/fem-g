// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_SCALARFUNCTION
#define H_GMSHFEM_SCALARFUNCTION

#include "Function.h"
#include "GeneralEvaluableObject.h"
#include "numa.h"

namespace gmshfem::field
{
  template< class T_Scalar, field::Form >
  class Field;
}

namespace gmshfem::function
{
  template< class T_Scalar, Degree T_Degree >
  class ExecutionTreeWithDegree;
}

namespace gmshfem::function
{


  // ************************************
  // Function< T_Scalar, Degree Degree0 >
  // ************************************

  template< class T_Scalar >
  class Function< T_Scalar, Degree::Degree0 > : public GeneralEvaluableObject< T_Scalar, Degree::Degree0 >
  {
   protected:
    const ExecutionTreeWithDegree< T_Scalar, Degree::Degree0 > *_tree;

   public:
    Function();
    explicit Function(const ExecutionTreeWithDegree< T_Scalar, Degree::Degree0 > *tree);
    Function(const Function &other);
    Function(Function &&other);
    template< class T_SFINAE = T_Scalar, class = typename std::enable_if< scalar::IsComplex< T_SFINAE >::value >::type >
    Function(const Function< scalar::Precision< T_SFINAE >, Degree::Degree0 > &other);
    Function(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value);
    template< class T_SFINAE = T_Scalar, class = typename std::enable_if< scalar::IsComplex< T_SFINAE >::value >::type >
    Function(const typename MathObject< scalar::Precision< T_SFINAE >, Degree::Degree0 >::Object &value);
    Function(const field::Field< T_Scalar, field::Form::Form0 > &field);
    Function(const field::Field< T_Scalar, field::Form::Form3 > &field);

    Function &operator=(const Function &other);
    Function &operator=(Function &&other);
    template< class T_SFINAE = T_Scalar, class = typename std::enable_if< scalar::IsComplex< T_SFINAE >::value >::type >
    Function &operator=(const Function< scalar::Precision< T_SFINAE >, Degree::Degree0 > &other);
    Function &operator=(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &value);
    template< class T_SFINAE = T_Scalar, class = typename std::enable_if< scalar::IsComplex< T_SFINAE >::value >::type >
    Function &operator=(const typename MathObject< scalar::Precision< T_SFINAE >, Degree::Degree0 >::Object &value);
    Function &operator=(const field::Field< T_Scalar, field::Form::Form0 > &field);
    Function &operator=(const field::Field< T_Scalar, field::Form::Form3 > &field);

    virtual ~Function();

    virtual ExecutionTreeIterator begin() const;
    virtual ExecutionTreeIterator end() const;

    virtual const ExecutionTreeWithDegree< T_Scalar, Degree::Degree0 > *getTree() const;
    virtual const ExecutionTreeWithDegree< T_Scalar, Degree::Degree0 > *copy() const;

    virtual void evaluate(std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const;
    virtual void evaluate(typename MathObject< T_Scalar, Degree::Degree0 >::Object &value, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const;
    virtual common::Timer evaluate(std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree0 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree0 >::Object > > > &values, const domain::Domain &domain, const std::string &gauss) const;
    virtual bool isConstant(const std::pair< int, int > &entity) const;
    virtual void setPointEvaluation() const;

    void exportExecutionTree(const domain::Domain &domain, const std::string &name, const std::string &path = "");

    // GeneralEvaluableObject
    virtual const Function< T_Scalar, Degree::Degree0 > &getEvaluableFunction() const override;
  };

  namespace literals
  {


    Function< std::complex< double >, Degree::Degree0 > operator"" _cd_sf(const long double value);
    Function< double, Degree::Degree0 > operator"" _d_sf(const long double value);
    Function< std::complex< float >, Degree::Degree0 > operator"" _cf_sf(const long double value);
    Function< float, Degree::Degree0 > operator"" _f_sf(const long double value);


  } // namespace literals

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b);


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b);


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b);


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const Function< T_Scalar, Degree::Degree1 > &a, const Function< T_Scalar, Degree::Degree1 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const typename MathObject< T_Scalar, Degree::Degree1 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree1 > &a, const typename MathObject< T_Scalar, Degree::Degree1 >::Object &b);


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator/(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator/(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator/(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree0 > operator/(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b);


} // namespace gmshfem::function

#endif // H_GMSHFEM_SCALARFUNCTION
