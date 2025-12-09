// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MATRIXMODULE
#define H_GMSHFEM_MATRIXMODULE

#include "Memory.h"
#include "scalar.h"

#include <string>
#include <vector>

namespace gmshfem::system
{
  template< class T_Scalar >
  class MatrixFactory;
}

namespace gmshfem::system
{


  template< class T_Scalar >
  class MatrixModule
  {
   protected:
    MatrixFactory< T_Scalar > *const _parentMatrix;

   public:
    MatrixModule(MatrixFactory< T_Scalar > *parentMatrix);
    virtual ~MatrixModule();

    virtual void setToZero() = 0;
    virtual std::string name() const = 0;
    virtual common::Memory memory() const = 0;

    virtual const T_Scalar &operator[](const unsigned int index) const = 0;
    virtual T_Scalar &operator[](const unsigned int index) = 0;
    virtual void activate(const char &param) const = 0;
    virtual const T_Scalar *getMatrix() const = 0;
  };


  template< class T_Scalar >
  class AModule : public MatrixModule< T_Scalar >
  {
   protected:
    std::vector< T_Scalar > _a;

   public:
    AModule(MatrixFactory< T_Scalar > *parentMatrix);
    virtual ~AModule();

    virtual void setToZero() override;
    virtual std::string name() const override;
    virtual common::Memory memory() const override;

    virtual const T_Scalar &operator[](const unsigned int index) const override;
    virtual T_Scalar &operator[](const unsigned int index) override;
    virtual void activate(const char &param) const override;
    virtual const T_Scalar *getMatrix() const override;
  };


  template< class T_CPScalar >
  class AFrequencyModule : public MatrixModule< T_CPScalar >
  {
   protected:
    std::vector< T_CPScalar > _a;
    scalar::Precision< T_CPScalar > _frequency;

   public:
    AFrequencyModule(MatrixFactory< T_CPScalar > *parentMatrix);
    AFrequencyModule(MatrixFactory< T_CPScalar > *parentMatrix, std::vector< T_CPScalar > &matrix);
    virtual ~AFrequencyModule();

    virtual void setToZero() override;
    virtual std::string name() const override;
    virtual common::Memory memory() const override;

    virtual const T_CPScalar &operator[](const unsigned int index) const override;
    virtual T_CPScalar &operator[](const unsigned int index) override;
    virtual void activate(const char &param) const override;
    virtual const T_CPScalar *getMatrix() const override;

    void setFrequency(const scalar::Precision< T_CPScalar > &frequency);
    scalar::Precision< T_CPScalar > getFrequency() const;
  };


  template< class T_Scalar >
  class MKModule : public MatrixModule< T_Scalar >
  {
   protected:
    std::vector< std::vector< T_Scalar > > _a; // M, K
    mutable int _activeComponents;

   public:
    MKModule(MatrixFactory< T_Scalar > *parentMatrix);
    virtual ~MKModule();

    virtual void setToZero() override;
    virtual std::string name() const override;
    virtual common::Memory memory() const override;

    virtual const T_Scalar &operator[](const unsigned int index) const override;
    virtual T_Scalar &operator[](const unsigned int index) override;
    virtual void activate(const char &param) const override;
    virtual const T_Scalar *getMatrix() const override;
  };


  template< class T_Scalar >
  class CKModule : public MatrixModule< T_Scalar >
  {
   protected:
    std::vector< std::vector< T_Scalar > > _a; // C, K
    mutable int _activeComponents;

   public:
    CKModule(MatrixFactory< T_Scalar > *parentMatrix);
    virtual ~CKModule();

    virtual void setToZero() override;
    virtual std::string name() const override;
    virtual common::Memory memory() const override;

    virtual const T_Scalar &operator[](const unsigned int index) const override;
    virtual T_Scalar &operator[](const unsigned int index) override;
    virtual void activate(const char &param) const override;
    virtual const T_Scalar *getMatrix() const override;
  };


  template< class T_Scalar >
  class MCKModule : public MatrixModule< T_Scalar >
  {
   protected:
    std::vector< std::vector< T_Scalar > > _a; // M, C, K
    mutable int _activeComponents;

   public:
    MCKModule(MatrixFactory< T_Scalar > *parentMatrix);
    virtual ~MCKModule();

    virtual void setToZero() override;
    virtual std::string name() const override;
    virtual common::Memory memory() const override;

    virtual const T_Scalar &operator[](const unsigned int index) const override;
    virtual T_Scalar &operator[](const unsigned int index) override;
    virtual void activate(const char &param) const override;
    virtual const T_Scalar *getMatrix() const override;
  };


} // namespace gmshfem::system


#endif // H_GMSHFEM_MATRIXMODULE
