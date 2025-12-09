// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_HELMHOLTZ2D
#define H_GMSHFEM_HELMHOLTZ2D

#include "AnalyticalNode.h"

namespace gmshfem::analytics::helmholtz2D
{


  //
  //  Scattering by a soft or hard cylinder:
  //   * Wave sign convention : e^iwt
  //   * Plane wave : e^ikx
  //

  //**************************
  // ScatteringByASoftCylinder
  //**************************

  template< class T_Scalar >
  class ScatteringByASoftCylinder final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    ScatteringByASoftCylinder(const ScatteringByASoftCylinder &other);
    ScatteringByASoftCylinder(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    ~ScatteringByASoftCylinder();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //**************************
  // ScatteringByAHardCylinder
  //**************************

  template< class T_Scalar >
  class ScatteringByAHardCylinder final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    ScatteringByAHardCylinder(const ScatteringByAHardCylinder &other);
    ScatteringByAHardCylinder(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    ~ScatteringByAHardCylinder();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //*****************************
  // dr_ScatteringByAHardCylinder
  //*****************************

  template< class T_Scalar >
  class dr_ScatteringByAHardCylinder final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const unsigned int _nbrOfterm;
    const scalar::Precision< T_Scalar > _thetaInc;

   public:
    dr_ScatteringByAHardCylinder(const dr_ScatteringByAHardCylinder &other);
    dr_ScatteringByAHardCylinder(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int nbrOfterm = 10, const scalar::Precision< T_Scalar > thetaInc = 0.);
    ~dr_ScatteringByAHardCylinder();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //**************************
  // DuctModeSolution
  //**************************

  template< class T_Scalar >
  class DuctModeSolution final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k; // Free field wavenumber
    const scalar::Precision< T_Scalar > _M; // Mach number
    const scalar::Precision< T_Scalar > _H; // Duct height
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const unsigned int _n; // Mode number
    const bool _BC; // Boundary condition

   public:
    DuctModeSolution(const DuctModeSolution &other);
    DuctModeSolution(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > M, const scalar::Precision< T_Scalar > H, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int n = 0, const bool BC = 0);
    ~DuctModeSolution();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //**************************
  // DuctModeSolutionAiry
  //**************************

  template< class T_Scalar >
  class DuctModeSolutionAiry final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k; // Free field wavenumber
    const scalar::Precision< T_Scalar > _H; // Duct height
    const scalar::Precision< T_Scalar > _a; // linear profile parameters ax+b
    const scalar::Precision< T_Scalar > _b;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const unsigned int _n; // Mode number
    const bool _BC; // Boundary condition - 0 Nemann BC (hard wall), 1 Dirichlet BC (soft wall)

   public:
    DuctModeSolutionAiry(const DuctModeSolutionAiry &other);
    DuctModeSolutionAiry(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > H, const scalar::Precision< T_Scalar > a, const scalar::Precision< T_Scalar > b, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const unsigned int n = 0, const bool BC = 0);
    ~DuctModeSolutionAiry();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };


  //**************************
  // DuctModeSolution Multimodal
  //**************************

  template< class T_Scalar >
  class DuctModeSolutionMulti final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k; // Free field wavenumber
    const scalar::Precision< T_Scalar > _M; // Mach number
    const scalar::Precision< T_Scalar > _H; // Duct height
    const std::vector< scalar::Precision< T_Scalar > > _A; // Modal amplitudes
    const std::vector< unsigned int > _n; // Modal indices
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const bool _BC; // Boundary condition

   public:
    DuctModeSolutionMulti(const DuctModeSolutionMulti &other);
    DuctModeSolutionMulti(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > M, const scalar::Precision< T_Scalar > H, const std::vector< scalar::Precision< T_Scalar > > A, const std::vector< unsigned int > n, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const bool BC = 0);
    ~DuctModeSolutionMulti();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //**************************
  // MonopoleFreeField
  //**************************

  template< class T_Scalar >
  class MonopoleFreeField final : public function::AnalyticalOperation< T_Scalar, Degree::Degree0 >
  {
   private:
    const scalar::Precision< T_Scalar > _k; // Free field wavenumber
    const scalar::Precision< T_Scalar > _M; // Mach number
    const scalar::Precision< T_Scalar > _theta; // Flow angle
    const scalar::Precision< T_Scalar > _A; // Source amplitude
    const scalar::Precision< T_Scalar > _xs; // Source x-position
    const scalar::Precision< T_Scalar > _ys; // Source y-position
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;

   public:
    MonopoleFreeField(const MonopoleFreeField &other);
    MonopoleFreeField(const scalar::Precision< T_Scalar > k, const scalar::Precision< T_Scalar > M, const scalar::Precision< T_Scalar > theta, const scalar::Precision< T_Scalar > A, const scalar::Precision< T_Scalar > xs, const scalar::Precision< T_Scalar > ys, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0.);
    ~MonopoleFreeField();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree0 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };


} // namespace gmshfem::analytics::helmholtz2D

#endif // H_GMSHFEM_HELMHOLTZ2D
