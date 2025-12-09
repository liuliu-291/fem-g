// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_HELMHOLTZ3D
#define H_GMSHFEM_HELMHOLTZ3D

#include "AnalyticalNode.h"

namespace gmshfem::analytics::helmholtz3D
{


  //************************
  // ScatteringByASoftSphere
  //************************

  template< class T_Scalar >
  class ScatteringByASoftSphere final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const scalar::Precision< T_Scalar > _z;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc; // azimuthal angle [0, 2*pi]
    const scalar::Precision< T_Scalar > _phiInc; // polar angle [0, pi]

   public:
    ScatteringByASoftSphere(const ScatteringByASoftSphere &other);
    ScatteringByASoftSphere(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const scalar::Precision< T_Scalar > z = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0., const scalar::Precision< T_Scalar > phiInc = 3.14159265358979323846 / 2.);
    ~ScatteringByASoftSphere();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //************************
  // ScatteringByAHardSphere
  //************************

  template< class T_Scalar >
  class ScatteringByAHardSphere final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const scalar::Precision< T_Scalar > _z;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc; // azimuthal angle [0, 2*pi]
    const scalar::Precision< T_Scalar > _phiInc; // polar angle [0, pi]

   public:
    ScatteringByAHardSphere(const ScatteringByAHardSphere &other);
    ScatteringByAHardSphere(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const scalar::Precision< T_Scalar > z = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0., const scalar::Precision< T_Scalar > phiInc = 3.14159265358979323846 / 2.);
    ~ScatteringByAHardSphere();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //***************************
  // dr_ScatteringByAHardSphere
  //***************************

  template< class T_Scalar >
  class dr_ScatteringByAHardSphere final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const scalar::Precision< T_Scalar > _z;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc; // azimuthal angle [0, 2*pi]
    const scalar::Precision< T_Scalar > _phiInc; // polar angle [0, pi]

   public:
    dr_ScatteringByAHardSphere(const dr_ScatteringByAHardSphere &other);
    dr_ScatteringByAHardSphere(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const scalar::Precision< T_Scalar > z = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0., const scalar::Precision< T_Scalar > phiInc = 3.14159265358979323846 / 2.);
    ~dr_ScatteringByAHardSphere();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //***************************
  // Transmission scattering problem
  // We consider the diffraction of a plane monochromatic wave in the +z
  // direction by a sphere of radius R0.
  //***************************

  template< class T_Scalar >
  class transmissionScatteringByASphere final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _kint;
    const scalar::Precision< T_Scalar > _kout;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _rhoint;
    const scalar::Precision< T_Scalar > _rhoout;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const scalar::Precision< T_Scalar > _z;
    const unsigned int _nbrOfterm;

   public:
    transmissionScatteringByASphere(const transmissionScatteringByASphere &other);
    transmissionScatteringByASphere(const scalar::Precision< T_Scalar > kint, const scalar::Precision< T_Scalar > kout, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > rhoint, const scalar::Precision< T_Scalar > rhoout, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const scalar::Precision< T_Scalar > z = 0., const unsigned int nbrOfterm = 10);
    ~transmissionScatteringByASphere();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };


} // namespace gmshfem::analytics::helmholtz3D

#endif // H_GMSHFEM_HELMHOLTZ3D
