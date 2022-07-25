#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>

/*
Function: FT Direct correlation function
System: Hard Sphere
Approximation: Percus Yevick
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -c: FT direct correlation / type double / Range NA
*/
double ck_hs_py( const double  phi, const double  k) {
  double c;
  double k_func,ckdum,dum1,dum2,dum3,dum4,dumsin,dumcos;
  double k2,k3,k4,k6;
  dum1 = - pow( 1.0 + 2.0 * phi, 2 ) ;
  dum2 = 6.0 * phi * pow( ( 1.0 + ( phi / 2.0 ) ), 2 );
  dum3 = - 0.50 * phi * pow( ( 1.0 + ( 2.0 * phi ) ), 2 );
  dum4 = ( 1.0 - phi );
  dum4 = pow(dum4, 4);
  dum1 = dum1 / dum4;
  dum2 = dum2 / dum4;
  dum3 = dum3 / dum4;
  k2 = k * k;
  if (k > 0.0750 ) {
    dumsin = gsl_sf_sin( k );
    dumcos = gsl_sf_cos( k );
    k3 = k2 * k;
    k4 = k3 * k;
    k6 = k4 * k2;
    c = ( dum1 * ( dumsin - k * dumcos ) / ( k3 ) ) +
      ( dum2 * ( ( ( 2.0 * k ) * dumsin ) + ( ( - ( k2 ) + 2.0 ) * dumcos ) - 2.0 ) /
      ( k4 ) ) + ( dum3 * ( ( ( 4.0 * ( k3 ) - 24.0 * k ) * dumsin )
      + ( ( - ( k4 ) + 12.0 * ( k2 ) - 24.0 ) * dumcos ) + 24.0 ) / ( k6 ));
  }
  else{
    c = dum1 * ( (1.0 / 3.0) - ( k2 / 30.0 ) ) + dum2 * ( 0.250 -  ( k2 / 36.0 ) )
      + dum3 * ( ( 1.0 / 6.0 ) - ( k2 / 48.0 ) );
  }

  c = c * 4.0 * M_PI;
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
double is_hs_py( const double phi, const double k ){
  double is;
  is = ck_hs_py(phi,k);
  is = 1.0 - 6.0 * phi * is / M_PI;
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
double sk_hs_py( const double phi, const double k ){
  double s;
  s = 1.0 / is_hs_py(phi,k);
  return s;
}

/*
Function: FT Direct correlation function
System: Hard Sphere
Approximation: Percus Yevick + Verlet Weis
Ref: doi=10.1103/PhysRevA.5.939
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -c: FT direct correlation / type double / Range NA
*/
double ck_hs_vw( const double phi, const double k ){
  double c;
  double k_vw,phi_vw;
  phi_vw = phi*(1.0 - (phi / 16.0)); /* Density correction from Verlet-Weiss */
  k_vw = k * pow( ( phi_vw / phi ) ,  1.0 / 3.0  ); /* Wave vector correction from Verlet-Weiss */
  c = ck_hs_py(phi_vw,k_vw);
  return c;
}



/*
Function: Inverse of static structure factor
System: Hard Sphere
Approximation: Percus Yevick + Verlet Weis
Ref: doi=10.1103/PhysRevA.5.939
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -is: Inverse of static structure factor / type double / Range (0:infty)
*/
double is_hs_vw( const double phi, const double k ){
  double is;
  double k_vw,phi_vw;
  phi_vw = phi*(1.0 - (phi / 16.0)); /* Density correction from Verlet-Weiss */
  k_vw = k * pow( ( phi_vw / phi ) ,  1.0 / 3.0  ); /* Wave vector correction from Verlet-Weiss */
  is = is_hs_py(phi_vw,k_vw);
  return is;
}

/*
Function: Static structure factor
System: Hard Sphere
Approximation: Percus Yevick + Verlet Weis
Ref: doi=10.1103/PhysRevA.5.939
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -s: Static structure factor / type double / Range [0:infty)
*/
double sk_hs_vw( const double phi, const double k ){
  double s;
  s = 1.0 / is_hs_vw(phi,k);
  return s;
}
