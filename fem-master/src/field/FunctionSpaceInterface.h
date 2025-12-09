// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_FUNCTIONSPACEINTERFACE
#define H_GMSHFEM_FUNCTIONSPACEINTERFACE

#include "FieldObject.h"
#include "scalar.h"

namespace gmshfem::field
{


  template< class T_PScalar >
  class FunctionSpaceInterface
  {
   protected:
    void _errorIfToManyDofsByElements(const unsigned int nbrDofsByElement) const;

   public:
    FunctionSpaceInterface();
    virtual ~FunctionSpaceInterface();

    void draw(const bool derivative, const std::string &elementName, const unsigned int orientation = 0, const double meshSize = 0.1) const;

    virtual std::string getGmshFemName(bool derivative) const = 0;
    virtual std::string _getGmshName(bool derivative) const;
    virtual std::string getGmshFemOrientationName() const = 0;
    virtual bool isConstantByElements() const = 0;
    virtual bool isFormP() const = 0;
    virtual field::Form form() const = 0;

    virtual FunctionSpaceInterface *copy() const = 0;

    virtual unsigned int numberOfComponents(const bool derivative) const = 0;
    virtual unsigned int getBasisFunctions(const bool derivative, std::vector< T_PScalar > &functions, std::vector< int > &orientations, const std::vector< T_PScalar > &integrationPoints, const int elementType, const std::pair< int, int > &entity) const = 0;
    virtual unsigned int getBasisFunction(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation = 0) const;
    virtual unsigned int getBasisFunction_noexcept(const bool derivative, std::vector< T_PScalar > &functions, const std::vector< T_PScalar > &integrationPoints, const int elementType, const unsigned int orientation = 0) const noexcept = 0;
    virtual void getKeys(const bool withCoordinates, std::vector< int > &typeKeys, std::vector< unsigned long long > &entityKeys, std::vector< T_PScalar > &coordinates, const int elementType, const std::pair< int, int > &entity) const;
    virtual void getKeysOnElement(const bool withCoordinates, std::vector< int > &typeKeys, std::vector< unsigned long long > &entityKeys, std::vector< T_PScalar > &coordinates, const unsigned int elementTag) const;
    virtual bool getPeriodicKeys(const bool withCoordinates, std::vector< int > &slaveTypeKeys, std::vector< unsigned long long > &slaveEntityKeys, std::vector< T_PScalar > &slaveCoordinates, std::vector< int > &masterTypeKeys, std::vector< unsigned long long > &masterEntityKeys, std::vector< T_PScalar > &masterCoordinates, const int elementType, const std::pair< int, int > &entity) const;
    virtual unsigned int getNumberOfKeysByElement(const int elementType) const;
    virtual unsigned int getNumberOfKeysByNode() const;
    virtual unsigned int getNumberOfKeysByEdge() const;
    virtual unsigned int getNumberOfKeysByTriangularFace() const;
    virtual unsigned int getNumberOfKeysByQuadrangularFace() const;
    virtual unsigned int getNumberOfOrientations(const int elementType) const;
    virtual void getKeyInformation(const std::vector< int > &typeKeys, const std::vector< unsigned long long > &entityKeys, const int elementType, std::vector< std::pair< int, int > > &infoKey) const;
    virtual unsigned int getOrientationOfElement(const unsigned int elementTag) const;

    virtual bool operator==(const FunctionSpaceInterface< T_PScalar > &other) const;
    virtual bool operator!=(const FunctionSpaceInterface< T_PScalar > &other) const;
    virtual bool operator<(const FunctionSpaceInterface< T_PScalar > &other) const;

    virtual unsigned int getGaussInfo(const std::string &integrationType, const int elementType, std::vector< T_PScalar > &weights, std::vector< T_PScalar > &points) const;

    static unsigned int GetGaussInfo(const std::string &integrationType, const int elementType, std::vector< T_PScalar > &weights, std::vector< T_PScalar > &points);
  };

  unsigned int getDimOfElementType(const int elementType);
  std::vector< int > getDenseOrientations(std::vector< int > &orientations, const unsigned int numberOfOrientations);


} // namespace gmshfem::field


#include "FunctionSpace.h"

#endif // H_GMSHFEM_FUNCTIONSPACEINTERFACE
