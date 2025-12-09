// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MAXWELL3D
#define H_GMSHFEM_MAXWELL3D

#include "AnalyticalNode.h"

namespace gmshfem::analytics::maxwell3D
{


  //
  // Diffraction by a conducting sphere
  //   * Plane wave : (e^ikz,0,0)
  // Maxwell' equations :
  // curl H = -k1 E
  // curl E = k2 H
  // k² =-k1*k2

  //***************************
  // Diffraction by a conducting sphere, J=H^n
  //***************************

  template< class T_Scalar >
  class JElectricSurfaceCurrent final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    const scalar::Precision< T_Scalar > _k_int;
    const scalar::Precision< T_Scalar > _k_out;
    const scalar::Precision< std::complex< T_Scalar > > _k2_int;
    const scalar::Precision< std::complex< T_Scalar > > _k2_out;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const scalar::Precision< T_Scalar > _z;
    const unsigned int _nbrOfterm;

   public:
    JElectricSurfaceCurrent(const scalar::Precision< T_Scalar > k_int, const scalar::Precision< T_Scalar > k_out, const scalar::Precision< std::complex< T_Scalar > > k2_int, const scalar::Precision< std::complex< T_Scalar > > k2_out, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const scalar::Precision< T_Scalar > z = 0., const unsigned int nbrOfterm = 15);
    JElectricSurfaceCurrent(const JElectricSurfaceCurrent &other);
    ~JElectricSurfaceCurrent();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };

  //***************************
  // Diffraction by a conducting sphere, M=E^n
  //***************************

  template< class T_Scalar >
  class MMagneticSurfaceCurrent final : public function::AnalyticalOperation< T_Scalar, Degree::Degree1 >
  {
   private:
    const scalar::Precision< T_Scalar > _k_int;
    const scalar::Precision< T_Scalar > _k_out;
    const scalar::Precision< std::complex< T_Scalar > > _k2_int;
    const scalar::Precision< std::complex< T_Scalar > > _k2_out;
    const scalar::Precision< T_Scalar > _R0;
    const scalar::Precision< T_Scalar > _x;
    const scalar::Precision< T_Scalar > _y;
    const scalar::Precision< T_Scalar > _z;
    const unsigned int _nbrOfterm;

   public:
    MMagneticSurfaceCurrent(const scalar::Precision< T_Scalar > k_int, const scalar::Precision< T_Scalar > k_out, const scalar::Precision< std::complex< T_Scalar > > k2_int, const scalar::Precision< std::complex< T_Scalar > > k2_out, const scalar::Precision< T_Scalar > R0, const scalar::Precision< T_Scalar > x = 0., const scalar::Precision< T_Scalar > y = 0., const scalar::Precision< T_Scalar > z = 0., const unsigned int nbrOfterm = 15);
    MMagneticSurfaceCurrent(const MMagneticSurfaceCurrent &other);
    ~MMagneticSurfaceCurrent();

    void operator()(function::OutputVector< T_Scalar, Degree::Degree1 > &values, const std::vector< scalar::Precision< T_Scalar >, numa::allocator< scalar::Precision< T_Scalar > > > &points) const override;
  };


} // namespace gmshfem::analytics::maxwell3D

#endif // H_GMSHFEM_MAXWELL3D
