#include "./Hard_Disk.h"

double ck_hd_roth( const double phi, const double k ){
  double c;
  double j0,j1;
  double kh,k2;
  double etanp1,etanp2;
  double dum1, dum2, dum3;
  kh= k / 2.0; k2= k * k;
  etanp1 = 1.0 - phi;
  etanp2 = etanp1 * etanp1;
  j0 = gsl_sf_bessel_J0(kh); j1 = gsl_sf_bessel_J1(kh);
  dum1 = - ( ( 5.0 / 4.0 ) * etanp2 * k2 * j0 * j0 ) ;
  dum2 = ( 4.0 * ( ( phi - 20.0 ) * phi + 7.0 ) + ( 5.0 * etanp2 * k2 / 4.0 ) ) * j1 * j1 ;
  dum3 = ( 2.0 * ( phi - 13.0 ) * etanp1 * k * j0 * j1 );
  c = M_PI * ( dum1 + dum2 + dum3 ) / ( 6.0 * ( etanp1 * etanp2 ) * k2 ) ;
  return c;
}

double is_hd_roth( const double phi, const double k ){
  double is;
  double rho;
  rho = phi * 4.0 / M_PI;
  is = 1.0 - ( rho * ck_hd_roth( phi, k ) );
  return is;
}

double sk_hd_roth( const double phi, const double k ){
  double s;
  s = 1.0 / is_hd_roth( phi, k);
  return s;
}
