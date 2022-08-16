#include "Hard_Sphere_Double_Yukawa.h"

/*
  For a Sharma-Sharma schme we want to compute for the FT of a double exponential potential, in particular, the 
  exponential potential reads:

  u(r) = u₀ exp[-z* (r-σ) / σ ] for r ≥ σ,

  whereas its FT is expressed as:

  βu(k)= 4πβuₐσ exp(zₐ) ∫𝛔∞ exp(-zₐr/σ)r sin(kr) dr - 4πβuᵣσ exp(zᵣ)∫₁∞ exp(-zᵣr/σ)r sin(kr) dr
  
  wherethe sub-index ᵣ indicates a repulsive potential, and ₐ an attractive potential. 
  By taking into account that the integral:
  
  ∫𝛔∞ exp(-zr/σ)r sin(kr) dr = 
  
  exp(-z) { [ -k ( z² + 2z + k²)cos(k) ] cos(k) - [z³+z²+(z-1)k²] sin(k) } / k(z² +k²)² ,

  thus, we need to ask for Ta=βuₐ, Tr=βuᵣ, zₐ, zᵣ and k to compute for βu(k)
*/


/* Function that computes the 3D FT of a repulsive exponential potential with 
range parameter z and scaled temperature T at a wave-vector magnitude k */
double
beta_u_k_yuk(const double T, const double z, const double k){
  double uk, k2, k3, z2, z3, sink_k, kcosk_k, denominator;
  uk=0.0; /* Default value */
  if ( T > 0.0 && z > 0.0 && k >= 0.0 ) {
    k2 = k * k;
    z2 = z * z;
    denominator = z2 + k2;
    denominator *= T;
    /* small wave-vector expansion */
    if ( k > 0.075 ) {
      sink_k = sin(k) / k;
      kcosk_k = cos(k);
    }
    else {
      sink_k  = 1.0 - ( k2 / 6.0 ) + ( k2 * k2 / 120.0 ) ;
      kcosk_k = 1.0 - ( k2 / 2.0 ) + ( k2 * k2 / 24.0  ) ;
    }
    
    uk =  4.0 * M_PI  * (  kcosk_k + z * sink_k ) / denominator;
  }
  return uk;
}

double
beta_u_k_double_yuk( const double Ta, const double Tr, const double za, const double zr, const double k  ){
    return beta_u_k_yuk( Tr, zr, k ) - beta_u_k_yuk( Ta, za, k ); /* Note, T=0 turns off the potential for either part */
}

/*
Function: FT Direct correlation function
System: Hard Sphere + double exponential
Approximation: Verlet-Weiss + Sharma-Sharma
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -Ta* = Dimensionless attractive temperature / Domain (0:\infty)
  -Tr* = Dimensionless repulsive temperature / Domain (0:\infty)
  -za  = range of attraction parameter / Domain (0:infty)
  -zr  = range of repulsion parameter / Domain (0:infty)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -c: FT direct correlation / type double / Range NA
*/
double ck_dble_yukawa_vwsh( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ) {
  double c, chs;
  chs = ck_hs_vw(phi,k);
  c = chs - beta_u_k_double_yuk( Ta, Tr, za, zr, k  );
  return c;
}

/*
Function: Inverse of static structure factor
System: Hard Sphere + double exponential
Approximation: Verlet-Weiss + Sharma-Sharma
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -Ta* = Dimensionless attractive temperature / Domain (0:\infty)
  -Tr* = Dimensionless repulsive temperature / Domain (0:\infty)
  -za  = range of attraction parameter / Domain (0:infty)
  -zr  = range of repulsion parameter / Domain (0:infty)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -is: Inverse of static structure factor / type double / Range (0:infty)
*/
double is_dble_yukawa_vwsh( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ){
  double is, phi_vw, chs;
  phi_vw = phi*(1.0 - (phi / 16.0));
  chs = ck_hs_vw(phi,k);
  is = (phi_vw * chs) - (phi * beta_u_k_double_yuk( Ta, Tr, za, zr, k  ));
  is = 1.0 - 6.0 * is / M_PI;
  return is;
}

/*
Function: Static structure factor
System: Hard Sphere + double exponential
Approximation: Verlet-Weiss + Sharma-Sharma
Inputs:
  -phi: Volume fraction / type double / Domain (0:1)
  -Ta* = Dimensionless attractive temperature / Domain (0:\infty)
  -Tr* = Dimensionless repulsive temperature / Domain (0:\infty)
  -za  = range of attraction parameter / Domain (0:infty)
  -zr  = range of repulsion parameter / Domain (0:infty)
  -k: wave vector magnitude / type double / Domain [0:infty)
Outputs:
  -s: Static structure factor / type double / Range [0:infty)
*/
double sk_dble_yukawa_vwsh( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ){
  double s;
  s = 1.0 / is_dble_yukawa_vwsh( phi, Ta, Tr, za, zr, k  );
  return s;
}

/* The wrong way of implementing Verlet Weiss for all the above functions... */
double ck_dble_yukawa_vwsh2( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ) {
  double c, chs, k_vw, phi_vw;
  phi_vw = phi*(1.0 - (phi / 16.0));
  k_vw = k * pow( ( phi_vw / phi ) ,  1.0 / 3.0  );
  chs = ck_hs_py(phi_vw,k_vw);
  c = chs - beta_u_k_double_yuk( Ta, Tr, za, zr, k_vw  );
  return c;
}
double is_dble_yukawa_vwsh2( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ){
  double is, phi_vw;
  phi_vw = phi*(1.0 - (phi / 16.0));
  is = phi_vw * ck_dble_yukawa_vwsh2(phi, Ta, Tr, za, zr, k  );
  is = 1.0 - 6.0 * is / M_PI;
  return is;
}
double sk_dble_yukawa_vwsh2( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ){
  double s;
  s = 1.0 / is_dble_yukawa_vwsh2( phi, Ta, Tr, za, zr, k  );
  return s;
}
