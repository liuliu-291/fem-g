// GetDP - Copyright (C) 1997-2020 P. Dular and C. Geuzaine, University of Liege
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/getdp/getdp/issues.

#include "Bessel.h"
#include <iostream>

#if defined(HAVE_NO_FORTRAN)

static void zbesj_(double*, double*, double*, int*, int*, double*,
    double*, int*, int*)
{
  std::cerr << "Bessel functions require Fortran compiler" << std::endl;
}

static void zbesk_(double*, double*, double*, int*, int*, double*,
    double*, int*, int*)
{
  std::cerr << "Bessel functions require Fortran compiler" << std::endl;
}

static void zbesy_(double*, double*, double*, int*, int*, double*,
    double*, int*, double*, double*, int*)
{
  std::cerr << "Bessel functions require Fortran compiler" << std::endl;
}

static void zbesh_(double*, double*, double*, int*, int*, int*,
    double*, double*, int*, int*)
{
  std::cerr << "Bessel functions require Fortran compiler" << std::endl;
}

static void zairy_(double*, double*, int*, int*, double*,
    double*, int*, int*)
{
  std::cerr << "Bessel functions require Fortran compiler" << std::endl;
}

#else

#if defined(HAVE_UNDERSCORE)
#define zbesj_ zbesj
#define zbesk_ zbesk
#define zbesy_ zbesy
#define zbesh_ zbesh
#define zairy_ zairy
#endif

extern "C" {
  void zbesj_(double*, double*, double*, int*, int*, double*,
      double*, int*, int*);
  void zbesk_(double*, double*, double*, int*, int*, double*,
      double*, int*, int*);

  void zbesy_(double*, double*, double*, int*, int*, double*,
      double*, int*, double*,
      double*, int*);
  void zbesh_(double*, double*, double*, int*, int*, int*, double*,
      double*, int*, int*);

  void zairy_(double*, double*, int*, int*, double*,
      double*, int*, int*);
}

#endif

static int BesselError(int ierr, const char *str)
{
  static int warn=0;

  switch(ierr){
    case 0 :
      return 0;
    case 1 :
      std::cerr << "Input error in " << str << std::endl;
      return BESSEL_ERROR_INPUT;
    case 2 :
      return BESSEL_OVERFLOW;
    case 3 :
      if(!warn){
        std::cout << "Half machine accuracy lost in " << str << " (large argument or order)" << std::endl;
        warn = 1;
      }
      return BESSEL_HALF_ACCURACY;
    case 4 :
      std::cerr << "Complete loss of significance in " << str << " (argument or order too large)" << std::endl;
      return BESSEL_NO_ACCURACY;
    case 5 :
      std::cerr << "Failed to converge in " << str << std::endl;
      return BESSEL_NO_CONVERGENCE;
    default:
      std::cout << "Unknown Bessel status in " << str << " (" << ierr << ")" << std::endl;
      return ierr;
  }
}

// First kind Bessel functions

int BesselJn(double n, int num, double x, double *val)
{
  int nz = 0, ierr = 0, kode = 1;
  double xi = 0.0;
  double* ji = new double[num];

  zbesj_(&x, &xi, &n, &kode, &num, val, ji, &nz, &ierr) ;

  delete[] ji;

  return BesselError(ierr, "BesselJn");
}

int BesselJnComplex(double n, int num, double xr, double xi, double *valr, double *vali)
{
  int nz = 0, ierr = 0, kode = 1;

  zbesj_(&xr, &xi, &n, &kode, &num, valr, vali, &nz, &ierr) ;

  return BesselError(ierr, "BesselJnComplex");
}

int BesselKnComplex(double n, int num, double xr, double xi, double *valr, double *vali)
{
  int nz = 0, ierr = 0, kode = 1;

  zbesk_(&xr, &xi, &n, &kode, &num, valr, vali, &nz, &ierr) ;

  return BesselError(ierr, "BesselKnComplex");
}

int BesselSphericalJn(double n, int num, double x, double *val)
{
  int ierr = BesselJn(n+0.5, num, x, val);
  double coef = sqrt(0.5*M_PI/x);
  for(int i = 0; i < num; i++){
    val[i] *= coef;
  }
  return BesselError(ierr, "BesselSphericalJn");
}

int BesselAltSphericalJn(double n, int num, double x, double *val)
{
  int ierr = BesselJn(n+0.5, num, x, val);
  double coef = sqrt(0.5*M_PI*x);
  for(int i = 0; i < num; i++){
    val[i] *= coef;
  }
  return BesselError(ierr, "BesselAltSphericalJn");
}

// Second kind Bessel functions

int BesselYn(double n, int num, double x, double *val)
{
  int nz = 0, ierr = 0, kode = 1;
  double xi = 0.0;
  double* yi = new double[num];
  double* auxyr = new double[num];
  double* auxyi = new double[num];

  zbesy_(&x, &xi, &n, &kode, &num, val, yi, &nz, auxyr, auxyi, &ierr);

  delete[] yi;
  delete[] auxyr;
  delete[] auxyi;

  return BesselError(ierr, "BesselYn");
}

int BesselYnComplex(double n, int num, double xr, double xi, double *valr, double *vali)
{
  int nz = 0, ierr = 0, kode = 1;
  double* auxyr = new double[num];
  double* auxyi = new double[num];

  zbesy_(&xr, &xi, &n, &kode, &num, valr, vali, &nz, auxyr, auxyi, &ierr) ;

  delete[] auxyr;
  delete[] auxyi;

  return BesselError(ierr, "BesselYnComplex");
}

int BesselSphericalYn(double n, int num, double x, double *val)
{
  int ierr = BesselYn(n+0.5, num, x, val);
  double coef = sqrt(0.5*M_PI/x);
  for(int i = 0; i < num; i++){
    val[i] *= coef;
  }
  return BesselError(ierr, "BesselSphericalYn");
}

int BesselAltSphericalYn(double n, int num, double x, double *val)
{
  int ierr = BesselYn(n+0.5, num, x, val);
  double coef = sqrt(0.5*M_PI*x);
  for(int i = 0; i < num; i++){
    val[i] *= coef;
  }
  return BesselError(ierr, "BesselAltSphericalYn");
}

// Hankel functions (type = 1 or 2)

int BesselHn(int type, double n, int num, double x, std::complex<double> *val)
{
  int nz = 0, ierr = 0, kode = 1;
  double* hr = new double[num];
  double* hi = new double[num];
  double xi = 0.0;

  zbesh_(&x, &xi, &n, &kode, &type, &num, hr, hi, &nz, &ierr);

  for(int i=0; i < num; i++){
    val[i] = std::complex<double>(hr[i], hi[i]);
  }

  delete[] hr;
  delete[] hi;

  return BesselError(ierr, "BesselHn");
}

int BesselSphericalHn(int type, double n, int num, double x, std::complex<double> *val)
{
  int ierr = BesselHn(type, n+0.5, num, x, val);
  double coef = sqrt(0.5*M_PI/x);
  for(int i = 0; i < num; i++){
    val[i] *= coef;
  }
  return BesselError(ierr, "BesselSphericalHn");
}

int BesselAltSphericalHn(int type, double n, int num, double x, std::complex<double> *val)
{
  int ierr = BesselHn(type, n+0.5, num, x, val);
  double coef = sqrt(0.5*M_PI*x);
  for(int i = 0; i < num; i++){
    val[i] *= coef;
  }
  return BesselError(ierr, "BesselAltSphericalHn");
}

// Airy Function and derivative (id=1)

int AiryComplex(double xr, double xi, int id, double *valr, double *vali)
{
  int nz = 0, ierr = 0, kode = 1;

  zairy_(&xr, &xi, &id, &kode, valr, vali, &nz, &ierr);
  //val = std::complex<double>(valr, vali);

  return BesselError(ierr, "AiryComplex");
}

// Utilities for backward compatibility

double Spherical_j_n(int n, double x)
{
  double res;
  BesselSphericalJn(n, 1, x, &res);
  return res;
}

double AltSpherical_j_n(int n, double x)
{
  double res;
  BesselAltSphericalJn(n, 1, x, &res);
  return res;
}

double Spherical_y_n(int n, double x)
{
  double res;
  BesselSphericalYn(n, 1, x, &res);
  return res;
}

double AltSpherical_y_n(int n, double x)
{
  double res;
  BesselAltSphericalYn(n, 1, x, &res);
  return res;
}
