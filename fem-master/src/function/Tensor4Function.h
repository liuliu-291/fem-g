// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_TENSOR4FUNCTION
#define H_GMSHFEM_TENSOR4FUNCTION

#include "Function.h"
#include "GeneralEvaluableObject.h"
#include "numa.h"

namespace gmshfem::function
{
  template< class T_Scalar, Degree T_Degree >
  class ExecutionTreeWithDegree;
}

namespace gmshfem::function
{


  // **********************************
  // Function< T_Scalar, Degree Degree4 >
  // **********************************

  template< class T_Scalar >
  class Function< T_Scalar, Degree::Degree4 > : public GeneralEvaluableObject< T_Scalar, Degree::Degree4 >
  {
   protected:
    const ExecutionTreeWithDegree< T_Scalar, Degree::Degree4 > *_tree;

   public:
    Function();
    explicit Function(const ExecutionTreeWithDegree< T_Scalar, Degree::Degree4 > *tree);
    Function(const Function &other);
    Function(Function &&other);
    template< class T_SFINAE = T_Scalar, class = typename std::enable_if< scalar::IsComplex< T_SFINAE >::value >::type >
    Function(const Function< scalar::Precision< T_SFINAE >, Degree::Degree4 > &other);
    Function(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &value);

    Function &operator=(const Function &other);
    Function &operator=(Function &&other);
    template< class T_SFINAE = T_Scalar, class = typename std::enable_if< scalar::IsComplex< T_SFINAE >::value >::type >
    Function &operator=(const Function< scalar::Precision< T_SFINAE >, Degree::Degree4 > &other);
    Function &operator=(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &value);

    virtual ~Function();

    virtual ExecutionTreeIterator begin() const;
    virtual ExecutionTreeIterator end() const;

    virtual const ExecutionTreeWithDegree< T_Scalar, Degree::Degree4 > *getTree() const;
    virtual const ExecutionTreeWithDegree< T_Scalar, Degree::Degree4 > *copy() const;

    virtual void evaluate(std::vector< typename MathObject< T_Scalar, Degree::Degree4 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree4 >::Object > > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points, const std::vector< scalar::Precision< T_Scalar > > &gaussPoints, const int elementType, const std::pair< int, int > &entity) const;
    virtual void evaluate(typename MathObject< T_Scalar, Degree::Degree4 >::Object &value, const scalar::Precision< T_Scalar > x, const scalar::Precision< T_Scalar > y, const scalar::Precision< T_Scalar > z, const std::pair< int, int > &entity) const;
    virtual common::Timer evaluate(std::vector< std::vector< typename MathObject< T_Scalar, Degree::Degree4 >::Object, numa::allocator< typename MathObject< T_Scalar, Degree::Degree4 >::Object > > > &values, const domain::Domain &domain, const std::string &gauss) const;
    virtual bool isConstant(const std::pair< int, int > &entity) const;
    virtual void setPointEvaluation() const;

    void exportExecutionTree(const domain::Domain &domain, const std::string &name, const std::string &path = "");

    // GeneralEvaluableObject
    virtual const Function< T_Scalar, Degree::Degree4 > &getEvaluableFunction() const override;
  };


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator+(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree4 >::Object &b);


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator-(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree4 >::Object &b);


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const Function< T_Scalar, Degree::Degree0 > &a, const Function< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const typename MathObject< T_Scalar, Degree::Degree0 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &a, const typename MathObject< T_Scalar, Degree::Degree4 >::Object &b);


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b);
  
  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator*(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree4 >::Object &b);


  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator/(const Function< T_Scalar, Degree::Degree4 > &a, const Function< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator/(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator/(const typename MathObject< T_Scalar, Degree::Degree4 >::Object &a, const GeneralEvaluableObject< T_Scalar, Degree::Degree0 > &b);

  template< class T_Scalar >
  Function< T_Scalar, Degree::Degree4 > operator/(const GeneralEvaluableObject< T_Scalar, Degree::Degree4 > &a, const typename MathObject< T_Scalar, Degree::Degree0 >::Object &b);


} // namespace gmshfem::function

#endif // H_GMSHFEM_TENSOR4FUNCTION
