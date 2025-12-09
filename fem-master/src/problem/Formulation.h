// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FORMULATION
#define H_GMSHFEM_FORMULATION

#include "AlgebraicFunctions.h"
#include "DofsManager.h"
#include "Domain.h"
#include "Equation.h"
#include "FieldInterface.h"
#include "Matrix.h"
#include "Options.h"
#include "Pattern.h"
#include "Post.h"
#include "Timer.h"
#include "UnknownField.h"
#include "Vector.h"
#include "optionsEnums.h"

#include <complex>
#include <initializer_list>
#include <string>
#include <vector>

namespace gmshfem::system
{
  template< class T_Scalar >
  class MatrixFactory;

  template< class T_Scalar >
  class VectorFactory;

  template< class T_Scalar >
  class Solver;
} // namespace gmshfem::system

namespace gmshfem::term
{
  template< class T_Scalar >
  class Term;
}

namespace gmshfem::problem
{


  template< class T_Scalar >
  class Formulation
  {
   protected:
    const std::string _name;
    std::vector< term::Term< T_Scalar > * > _terms;
    dofs::DofsManager< T_Scalar > _dofs;
    system::MatrixFactory< T_Scalar > *_A;
    std::vector<system::VectorFactory<T_Scalar>> _multiB;
    system::Solver< T_Scalar > *_solver;
    std::vector< std::pair< unsigned int, field::FieldInterface< T_Scalar > * > > _unknownFields;
    gmsh::vectorpair _entities[4];
    std::vector< int > _elementTypes[4];
    scalar::Precision< T_Scalar > _frequency;
    std::unordered_map< std::string, void * > _attributes;

    void _addField(field::FieldInterface< T_Scalar > *field);
    common::Timer _assembleDim(const int dimToAssembly, const ElementsSort::Algorithm algo);
    common::Timer _buildPattern(const ElementsSort::Algorithm algo);
    std::vector< term::Term< T_Scalar > * > _getTermOnEntity(const std::pair< int, int > &entity) const;

    unsigned _currentRHS = 0;
    unsigned _numRHS = 1;
    std::vector<std::vector<T_Scalar>> _solutions;

   public:
    Formulation(const std::string &name, const std::string &matrixOptions = "");
    virtual ~Formulation();

    Formulation(Formulation< T_Scalar > &&other);

    virtual void initSystem(const std::string &matrixOptions = "");

    virtual void removeTerms();
    virtual void removeSystem();

    virtual void setSystemToZero();
    virtual void setRHSToZero();
    std::string name() const;

    typename std::vector< term::Term< T_Scalar > * >::iterator begin();
    typename std::vector< term::Term< T_Scalar > * >::iterator end();

    template< class T_Object >
    void setAttribute(const std::string &name, const T_Object &attribute);
    template< class T_Object >
    void getAttribute(const std::string &name, T_Object &attribute) const;

    // Multi RHS interface
    unsigned numRHS() const;
    void setCurrentRHS(unsigned idx);

    // Bilinear term
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);


    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);


    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);


    // Linear terms
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree0 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree0 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree0 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree0 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree0, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const function::Function< T_Scalar, Degree::Degree1 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree1 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree1 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree1 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree1, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);

    unsigned int integral(const function::Function< T_Scalar, Degree::Degree2 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form0 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree2 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form1 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree2 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form2 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);
    unsigned int integral(const function::Function< T_Scalar, Degree::Degree2 > &functionLhs, const equation::EquationInterface< T_Scalar, Degree::Degree2, field::Form::Form3 > &equationRhs, const domain::GeometricObject &domain, const std::string &integrationType, const term::ProductType productType = term::ProductType::Hermitian);


    // Global quantity terms
    void globalTerm(field::GlobalQuantity< T_Scalar > &globalQuantity, const field::FixedComponent &fixComponent, const T_Scalar &fixedValue);

    virtual unsigned long long getTotalNumberOfDof() const;
    virtual unsigned long long getNumberOfUnknownDof() const;
    virtual unsigned long long getNumberOfFixedDof() const;
    virtual void setAngularFrequency(const scalar::Precision< T_Scalar > &frequency);
    virtual scalar::Precision< T_Scalar > getAngularFrequency() const;
    virtual field::FieldInterface< T_Scalar > *getField(const std::string &name) const;
    common::Memory getEstimatedFactorizationMemoryUsage() const;
    virtual unsigned long long getNumberOfNonZeros() const;

    // Pre-processing
    virtual common::Timer pre(const DofsSort::Algorithm algo = common::Options::instance()->dofsSortAlgorithm);

    // Processing
    virtual common::Timer assemble(const bool separate = false, const ElementsSort::Algorithm algo = common::Options::instance()->elementsSortAlgorithm); // assemble A and b
    virtual common::Timer solveAll(const bool reusePreconditioner = false); // solve A X = B and store all solutions
    virtual void loadSolution(unsigned idx);
    virtual common::Timer solve(const bool reusePreconditioner = false); // solve A x = b
    virtual common::Timer solve(const bool reusePreconditioner,
                                algebra::Vector< T_Scalar > &rawSolution); // same as above, but exposes raw solution vector
    virtual common::Timer eigensolve(algebra::Vector< scalar::ComplexPrecision< T_Scalar > > &eigenvalues, const bool computeEigenvectors = false, const unsigned long long numberOfEigenvalues = 0, const scalar::ComplexPrecision< T_Scalar > target = 0); // solve A x = \lambda x, K x = \lambda M x, K x = \lambda C x or (\lambda^2 M + \lambda C + K) x = 0

    void setSolutionIntoFields(const std::vector< T_Scalar > &solution);

    virtual scalar::Precision< T_Scalar > getResidual() const; // return b - Ax

    void getLHS(algebra::Matrix< T_Scalar > &matrix) const;
    void getMass(algebra::Matrix< T_Scalar > &matrix) const;
    void getDamping(algebra::Matrix< T_Scalar > &matrix) const;
    void getStiffness(algebra::Matrix< T_Scalar > &matrix) const;

    void getLHSBlock(algebra::Matrix< T_Scalar > &matrix, const field::FieldInterface< T_Scalar > &dof, const field::FieldInterface< T_Scalar > &tf) const;
    void getLHSBlock(algebra::Matrix< T_Scalar > &matrix, const std::vector< const field::FieldInterface< T_Scalar > * > &dofs, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const;
    void getMassBlock(algebra::Matrix< T_Scalar > &matrix, const field::FieldInterface< T_Scalar > &dof, const field::FieldInterface< T_Scalar > &tf) const;
    void getMassBlock(algebra::Matrix< T_Scalar > &matrix, const std::vector< const field::FieldInterface< T_Scalar > * > &dofs, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const;
    void getDampingBlock(algebra::Matrix< T_Scalar > &matrix, const field::FieldInterface< T_Scalar > &dof, const field::FieldInterface< T_Scalar > &tf) const;
    void getDampingBlock(algebra::Matrix< T_Scalar > &matrix, const std::vector< const field::FieldInterface< T_Scalar > * > &dofs, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const;
    void getStiffnessBlock(algebra::Matrix< T_Scalar > &matrix, const field::FieldInterface< T_Scalar > &dof, const field::FieldInterface< T_Scalar > &tf) const;
    void getStiffnessBlock(algebra::Matrix< T_Scalar > &matrix, const std::vector< const field::FieldInterface< T_Scalar > * > &dofs, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const;

    void getRHS(algebra::Vector< T_Scalar > &vector) const;
    void getRHSBlock(algebra::Vector< T_Scalar > &vector, const field::FieldInterface< T_Scalar > &tf) const;
    void getRHSBlock(algebra::Vector< T_Scalar > &vector, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs) const;

    void setRHS(const algebra::Vector< T_Scalar > &vector);
    void setRHSBlock(const algebra::Vector< T_Scalar > &vector, const field::FieldInterface< T_Scalar > &tf);
    void setRHSBlock(const algebra::Vector< T_Scalar > &vector, const std::vector< const field::FieldInterface< T_Scalar > * > &tfs);
  };


} // namespace gmshfem::problem

#endif // H_GMSHFEM_FORMULATION
