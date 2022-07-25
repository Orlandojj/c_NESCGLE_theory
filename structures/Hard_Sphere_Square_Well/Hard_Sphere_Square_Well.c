#include "Hard_Sphere_Square_Well.h"

/*
Function: Auxiliary function that computes for the FT of the perturbation potential -\beta u_p(k)
System: Square Well
Inputs:
  -T* = Dimensionless temperature / Domain (0:\infty)
  -lambda range of the well in terms of the diameter sigma (1:infty)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  --\beta u_p(k): FT of the perturbation well
*/

double sw_m_beta_uk( const double T, const double lambda, const double k ){
  double c_aux;
  double k2,lk,l5,l3,coslk,sinlk,sink,cosk;
  k2 = k * k;
  if (k > 0.0750) {
    lk = lambda * k ;
    sink  = gsl_sf_sin( k );
    cosk  = gsl_sf_cos( k );
    sinlk = gsl_sf_sin( lk );
    coslk = gsl_sf_cos( lk );
    c_aux = ((cosk - lambda * coslk) / k) + ((sinlk - sink) / k2) ;
    c_aux = c_aux / k ;
  }
  else {
    l3 = lambda * lambda * lambda;
    l5 = l3 * lambda * lambda;
    c_aux = (1.0/3.0) * (l3 - 1.0) - (1.0/30.0) * (l5 - 1.0) * k2;
  }
  c_aux = 4.0 * M_PI * c_aux / T;
  return c_aux;
  }

/*
Function: FT Direct correlation function
System: Hard Sphere + Square Well
Approximation: Verlet-Weiss + Sharma-Sharma
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -T* = Dimensionless temperature / Domain (0:\infty)
  -lambda range of the well in terms of the diameter sigma / Domain (1:infty)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -c: FT direct correlation / type double / Range NA
*/
double ck_hssw_vwsh( const double phi, const double T, const double lambda, const double k ) {
  double c, chs;
  chs = ck_hs_vw(phi,k);
  c = chs + sw_m_beta_uk( T, lambda, k );
  return c;
}

/*
Function: Inverse of static structure factor
System: Hard Sphere
Approximation: Percus Yevick
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -is: Inverse of static structure factor / type double / Range (0:infty)
*/
double is_hssw_vwsh( const double phi, const double T, const double lambda, const double k ){
  double is, phi_vw, chs;
  phi_vw = phi*(1.0 - (phi / 16.0));
  chs = ck_hs_vw(phi,k);
  is = (phi_vw * chs) + (phi * sw_m_beta_uk(T,lambda,k));
  is = 1.0 - 6.0 * is / M_PI;
  return is;
}

/*
Function: Static structure factor
System: Hard Sphere
Approximation: Percus Yevick
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -s: Static structure factor / type double / Range [0:infty)
*/
double sk_hssw_vwsh( const double phi, const double T, const double lambda, const double k ){
  double s;
  s = 1.0 / is_hssw_vwsh(phi, T, lambda, k );
  return s;
}
