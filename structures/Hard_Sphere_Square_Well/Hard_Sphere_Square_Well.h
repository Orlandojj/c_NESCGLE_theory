#ifndef HARD_SPHERE_SQUARE_WELL_DOT_H    /* This is an "include guard" */
#define HARD_SPHERE_SQUARE_WELL_DOT_H    /* prevents the file from being included twice. */


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include "../Hard_Sphere/Hard_Sphere.h"

/* Sharma-Sharma w Verlet-Weiss approx Function */
/*
  phi = volume fraction
  T = Dimensionless temperature
  lambda = well length
  k =  wave vector magnitude
*/
double ck_hssw_vwsh( const double phi, const double T, const double lambda, const double k ); /* Direct correlation function c(k) */
double is_hssw_vwsh( const double phi, const double T, const double lambda, const double k ); /* Inverse of static structure factor S^-1(k) */
double sk_hssw_vwsh( const double phi, const double T, const double lambda, const double k ); /* Static structure factor S(k) */


#endif /* HARD_SPHERE_DOT_H */
