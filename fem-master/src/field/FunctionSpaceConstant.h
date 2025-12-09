// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FUNCTIONSPACECONSTANT
#define H_GMSHFEM_FUNCTIONSPACECONSTANT

#include "FunctionSpace.h"

#include <vector>

namespace gmshfem::field
{


  template< class T_PScalar >
  class FunctionSpaceConstant final : public FunctionSpace< T_PScalar, Form::Form3 >
  {
   private:
    const bool _pForm;

   public:
    FunctionSpaceConstant(const bool pForm);
    ~FunctionSpaceConstant();

    virtual std::string getGmshFemName(bool derivative) const override;
    virtual std::string getGmshFemOrientationName() const override;
    virtual FunctionSpaceTypeForm3 type() const override;
    virtual unsigned int order() const override;
    virtual bool isConstantByElements() const override;
    virtual bool isFormP() const override;

    virtual FunctionSpaceConstant *copy() const override;

    virtual bool operator==(const FunctionSpaceInterface< T_PScalar > &other) const override;

    virtual unsigned int getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientation, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const override;
    virtual unsigned int getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation = 0) const noexcept override;

    virtual void getKeys(const bool withCoordinates, std::vector< int > &typeKeys, std::vector< unsigned long long > &entityKeys, std::vector< T_PScalar > &coordinates, const int elementType, const std::pair< int, int > &entity) const override;
    virtual void getKeysOnElement(const bool withCoordinates, std::vector< int > &typeKeys, std::vector< unsigned long long > &entityKeys, std::vector< T_PScalar > &coordinates, const unsigned int elementTag) const override;
    virtual bool getPeriodicKeys(const bool withCoordinates, std::vector< int > &slaveTypeKeys, std::vector< unsigned long long > &slaveEntityKeys, std::vector< T_PScalar > &slaveCoordinates, std::vector< int > &masterTypeKeys, std::vector< unsigned long long > &masterEntityKeys, std::vector< T_PScalar > &masterCoordinates, const int elementType, const std::pair< int, int > &entity) const override;
    virtual unsigned int getNumberOfKeysByElement(const int elementType) const override;
    virtual unsigned int getNumberOfKeysByNode() const override;
    virtual unsigned int getNumberOfKeysByEdge() const override;
    virtual unsigned int getNumberOfKeysByTriangularFace() const override;
    virtual unsigned int getNumberOfKeysByQuadrangularFace() const override;
    virtual unsigned int getNumberOfOrientations(const int elementType) const override;
    virtual void getKeyInformation(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, const int elementType, std::vector< std::pair< int, int > > &infoKey) const override;
    virtual unsigned int getOrientationOfElement(const unsigned int elementTag) const override;
  };


} // namespace gmshfem::field

#endif // H_GMSHFEM_FUNCTIONSPACECONSTANT
