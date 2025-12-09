// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MAXWELL2D
#define H_GMSHFEM_MAXWELL2D

#include "AnalyticalNode.h"
#include "Helmholtz2D.h"

namespace gmshfem::analytics::maxwell2D
{


  //
  //  Scattering by a PEC (perfect electric conductor) or PMC (perfect magnetic conductor) cylinder:
  //   * Wave sign convention : e^-iwt
  //   * Plane wave : e^ikx
  //

  //***************************
  // ScatteringByAPECCylinderTE
  //***************************

  template< class T_Scalar >
  class ScatteringByAPECCylinderTE final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    const helmholtz2D::ScatteringByASoftCylinder< T_Scalar > _helmholtz;
    mutable function::OutputVector< T_Scalar, Degree::Degree0 > _values;

   public:
    ScatteringByAPECCylinderTE(const ScatteringByAPECCylinderTE &other);
    ScatteringByAPECCylinderTE(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    ~ScatteringByAPECCylinderTE();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //***************************
  // ScatteringByAPECCylinderTM
  //***************************

  template< class T_Scalar >
  class ScatteringByAPECCylinderTM final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    const scalar::Precision< T_Scalar > _k;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    ScatteringByAPECCylinderTM(const ScatteringByAPECCylinderTM &other);
    ScatteringByAPECCylinderTM(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    ~ScatteringByAPECCylinderTM();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //***************************
  // ScatteringByAPMCCylinderTE
  //***************************

  template< class T_Scalar >
  class ScatteringByAPMCCylinderTE final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    const helmholtz2D::ScatteringByAHardCylinder< T_Scalar > _helmholtz;
    mutable function::OutputVector< T_Scalar, Degree::Degree0 > _values;

   public:
    ScatteringByAPMCCylinderTE(const ScatteringByAPMCCylinderTE &other);
    ScatteringByAPMCCylinderTE(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    ~ScatteringByAPMCCylinderTE();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //********************************
  // dr_z_ScatteringByAPMCCylinderTE
  //********************************

  template< class T_Scalar >
  class dr_ScatteringByAPMCCylinderTE final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    const helmholtz2D::dr_ScatteringByAHardCylinder< T_Scalar > _helmholtz;
    mutable function::OutputVector< T_Scalar, Degree::Degree0 > _values;

   public:
    dr_ScatteringByAPMCCylinderTE(const dr_ScatteringByAPMCCylinderTE &other);
    dr_ScatteringByAPMCCylinderTE(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    ~dr_ScatteringByAPMCCylinderTE();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };


} // namespace gmshfem::analytics::maxwell2D

#endif // H_GMSHFEM_MAXWELL2D
