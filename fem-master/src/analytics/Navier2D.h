// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_NAVIER2D
#define H_GMSHFEM_NAVIER2D

#include "AnalyticalNode.h"

namespace gmshfem::analytics::navier2D
{


  //********************************
  // SoftPWavesScatteringByACylinder
  //********************************

  // Exact Navier solution for a cylindrical wall of radius _R0 (zero displacement on the boundary)
  // for a pressure incident wave : u_inc= ( i_kP*exp(-i_kP X),0 )^T.
  // The (theoric infinite) sum is truncated at _nbOfterm.

  template< class T_Scalar >
  class SoftPWavesScatteringByACylinder final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   protected:
    const scalar::Precision< T_Scalar > _kP;
    const scalar::Precision< T_Scalar > _kS;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    SoftPWavesScatteringByACylinder(const SoftPWavesScatteringByACylinder< T_Scalar > &other);
    SoftPWavesScatteringByACylinder(const scalar::Precision< T_Scalar > kP, const scalar::Precision< T_Scalar > kS, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    virtual ~SoftPWavesScatteringByACylinder();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //********************************
  // SoftSWavesScatteringByACylinder
  //********************************

  // Exact Navier solution for a cylindrical wall of radius _R0 (zero displacement on the boundary)
  // for a shear incident wave : u_inc= ( 0, -i_kS*exp(-i_kS X) )^T.
  // The (theoric infinite) sum is truncated at _nbOfterm.

  template< class T_Scalar >
  class SoftSWavesScatteringByACylinder final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   protected:
    const scalar::Precision< T_Scalar > _kP;
    const scalar::Precision< T_Scalar > _kS;
    const scalar::Precision< T_Scalar > _R0;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    SoftSWavesScatteringByACylinder(const SoftSWavesScatteringByACylinder< T_Scalar > &other);
    SoftSWavesScatteringByACylinder(const scalar::Precision< T_Scalar > kP, const scalar::Precision< T_Scalar > kS, const scalar::Precision< T_Scalar > R0, const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    virtual ~SoftSWavesScatteringByACylinder();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //*****************************************
  // SoftPWavesScatteringByACylinderWithLOABC
  //*****************************************

  // Truncated Navier solution for a cylindrical wall of radius _R0 (zero displacement on the boundary)
  // for a pressure incident wave : u_inc= ( i_kP*exp(-i_kP X),0 )^T.
  // LO-ABC onto the cylinder of radius _R1.
  // The (theoric infinite) sum is truncated at _nbOfterm.

  template< class T_Scalar >
  class SoftPWavesScatteringByACylinderWithLOABC final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   protected:
    const scalar::Precision< T_Scalar > _kP;
    const scalar::Precision< T_Scalar > _kS;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _R1;
    const scalar::Precision< T_Scalar > _lambda;
    const scalar::Precision< T_Scalar > _mu;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    SoftPWavesScatteringByACylinderWithLOABC(const SoftPWavesScatteringByACylinderWithLOABC< T_Scalar > &other);
    SoftPWavesScatteringByACylinderWithLOABC(const scalar::Precision< T_Scalar > kP, const scalar::Precision< T_Scalar > kS, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > R1, const scalar::Precision< T_Scalar > lambda, const scalar::Precision< T_Scalar > mu, const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    virtual ~SoftPWavesScatteringByACylinderWithLOABC();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //*****************************************
  // SoftSWavesScatteringByACylinderWithLOABC
  //*****************************************

  // Truncated Navier solution for a cylindrical wall of radius _R0 (zero displacement on the boundary)
  // for a shear incident wave : u_inc= ( 0, -i_kS*exp(-i_kS X) )^T.
  // LO-ABC onto the cylinder of radius _R1.
  // The (theoric infinite) sum is truncated at _nbOfterm.

  template< class T_Scalar >
  class SoftSWavesScatteringByACylinderWithLOABC final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    const scalar::Precision< T_Scalar > _kP;
    const scalar::Precision< T_Scalar > _kS;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _R1;
    const scalar::Precision< T_Scalar > _lambda;
    const scalar::Precision< T_Scalar > _mu;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    SoftSWavesScatteringByACylinderWithLOABC(const SoftSWavesScatteringByACylinderWithLOABC< T_Scalar > &other);
    SoftSWavesScatteringByACylinderWithLOABC(const scalar::Precision< T_Scalar > kP, const scalar::Precision< T_Scalar > kS, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > R1, const scalar::Precision< T_Scalar > lambda, const scalar::Precision< T_Scalar > mu, const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    virtual ~SoftSWavesScatteringByACylinderWithLOABC();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //*****************************************
  // SoftPWavesScatteringByACylinderWithHOABC
  //*****************************************

  // Truncated Navier solution for a cylindrical wall of radius _R0 (zero displacement on the boundary)
  // for a pressure incident wave : u_inc= ( i_kP*exp(-i_kP X),0 )^T.
  // HO-ABC onto the cylinder of radius _R1.
  // The (theoric infinite) sum is truncated at _nbOfterm.

  template< class T_Scalar >
  class SoftPWavesScatteringByACylinderWithHOABC final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   protected:
    const scalar::Precision< T_Scalar > _f;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _R1;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _lambda;
    const scalar::Precision< T_Scalar > _mu;
    const scalar::Precision< T_Scalar > _rho;
    const unsigned int _padeOrder;
    const scalar::Precision< T_Scalar > _padeAngle;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    SoftPWavesScatteringByACylinderWithHOABC(const SoftPWavesScatteringByACylinderWithHOABC< T_Scalar > &other);
    SoftPWavesScatteringByACylinderWithHOABC(const scalar::Precision< T_Scalar > f, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > R1, const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > lambda = 1., const scalar::Precision< T_Scalar > mu = 1., const scalar::Precision< T_Scalar > rho = 1., const unsigned int padeOrder = 4, const scalar::Precision< T_Scalar > padeAngle = 3.14159265359 / 4., const scalar::Precision< T_Scalar > thetaInc = 0.);
    virtual ~SoftPWavesScatteringByACylinderWithHOABC();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //*****************************************
  // SoftSWavesScatteringByACylinderWithHOABC
  //*****************************************

  // Truncated Navier solution for a cylindrical wall of radius _R0 (zero displacement on the boundary)
  // for a shear incident wave : u_inc= ( 0, -i_kS*exp(-i_kS X) )^T.
  // HO-ABC onto the cylinder of radius _R1.
  // The (theoric infinite) sum is truncated at _nbOfterm.

  template< class T_Scalar >
  class SoftSWavesScatteringByACylinderWithHOABC final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   protected:
    const scalar::Precision< T_Scalar > _f;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _R1;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _lambda;
    const scalar::Precision< T_Scalar > _mu;
    const scalar::Precision< T_Scalar > _rho;
    const unsigned int _padeOrder;
    const scalar::Precision< T_Scalar > _padeAngle;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    SoftSWavesScatteringByACylinderWithHOABC(const SoftSWavesScatteringByACylinderWithHOABC< T_Scalar > &other);
    SoftSWavesScatteringByACylinderWithHOABC(const scalar::Precision< T_Scalar > f, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > R1, const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > lambda = 1., const scalar::Precision< T_Scalar > mu = 1., const scalar::Precision< T_Scalar > rho = 1., const unsigned int padeOrder = 4, const scalar::Precision< T_Scalar > padeAngle = 3.14159265359 / 4., const scalar::Precision< T_Scalar > thetaInc = 0.);
    virtual ~SoftSWavesScatteringByACylinderWithHOABC();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };


} // namespace gmshfem::analytics::navier2D

#endif // H_GMSHFEM_NAVIER2D
