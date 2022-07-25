#ifndef HARD_DISK_DOT_H    /* This is an "include guard" */
#define HARD_DISK_DOT_H    /* prevents the file from being included twice. */

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>

/* Roth approx Functions */
double ck_hd_roth( const double phi, const double k );
double is_hd_roth( const double phi, const double k );
double sk_hd_roth( const double phi, const double k );

#endif /* HARD_SPHERE_DOT_H */
