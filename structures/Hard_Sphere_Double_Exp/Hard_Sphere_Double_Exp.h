#ifndef HARD_SPHERE_DOUBLE_EXP_DOT_H    /* This is an "include guard" */
#define HARD_SPHERE_DOUBLE_EXP_DOT_H    /* prevents the file from being included twice. */


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../Hard_Sphere/Hard_Sphere.h"

/* Sharma-Sharma w Verlet-Weiss approx Function */
/*
  phi = volume fraction
  T = Dimensionless temperature
  lambda = well length
  k =  wave vector magnitude
*/
double ck_dble_exp_vwsh( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ); /* Direct correlation function c(k) */
double is_dble_exp_vwsh( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ); /* Inverse of static structure factor S^-1(k) */
double sk_dble_exp_vwsh( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ); /* Static structure factor S(k) */
/* For the version 2 of functions the verlet weiss correction is also considered within the functions */
double ck_dble_exp_vwsh2( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ); /* Direct correlation function c(k) */
double is_dble_exp_vwsh2( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ); /* Inverse of static structure factor S^-1(k) */
double sk_dble_exp_vwsh2( const double phi, const double Ta, const double Tr, const double za, const double zr, const double k  ); /* Static structure factor S(k) */


#endif /* HARD_SPHERE_DOUBLE_EXP_DOT_H */
