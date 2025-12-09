// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FIELDMODULE
#define H_GMSHFEM_FIELDMODULE

#include "Memory.h"
#include "scalar.h"

#include <string>
#include <unordered_map>

namespace gmshfem::field
{
  template< class T_Scalar >
  class FieldInterface;
}

namespace gmshfem::field
{


  template< class T_Scalar >
  class FieldModule
  {
   protected:
    FieldInterface< T_Scalar > *const _parentField;

   public:
    FieldModule(FieldInterface< T_Scalar > *parentField);
    FieldModule(const FieldModule< T_Scalar > &other);
    virtual ~FieldModule();

    virtual void clear() = 0;
    virtual std::string name() const = 0;
    virtual common::Memory memory() const = 0;
    virtual FieldModule< T_Scalar > *copy() const = 0;
  };


  template< class T_Scalar >
  class EigenpairModule final : public FieldModule< T_Scalar >
  {
   protected:
    std::unordered_map< unsigned int, unsigned int > _mapping;
    std::vector< std::pair< scalar::ComplexPrecision< T_Scalar >, std::vector< scalar::ComplexPrecision< T_Scalar > > > > _eigenpairs;

   public:
    EigenpairModule(FieldInterface< T_Scalar > *parentField);
    EigenpairModule(const EigenpairModule< T_Scalar > &other);
    virtual ~EigenpairModule();

    virtual void clear() override;
    virtual std::string name() const override;
    virtual common::Memory memory() const override;
    virtual EigenpairModule< T_Scalar > *copy() const override;

    void assignEigenpair(const scalar::ComplexPrecision< T_Scalar > &eigenvalue, const std::vector< scalar::ComplexPrecision< T_Scalar > > &eigenvector);
    typename std::vector< std::pair< scalar::ComplexPrecision< T_Scalar >, std::vector< scalar::ComplexPrecision< T_Scalar > > > >::const_iterator begin() const;
    typename std::vector< std::pair< scalar::ComplexPrecision< T_Scalar >, std::vector< scalar::ComplexPrecision< T_Scalar > > > >::const_iterator end() const;

    void getValues(const unsigned int eigenTag, const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, std::vector< scalar::ComplexPrecision< T_Scalar > > &values, const unsigned int begin, const unsigned int end) const;
  };


} // namespace gmshfem::field


#endif // H_GMSHFEM_FIELDMODULE
