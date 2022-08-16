/* This header serves to import all structures headers defined in each
system folder. For the imported impl */
#ifndef STRUCTURES_DOT_H    /* This is an "include guard" */
#define STRUCTURES_DOT_H    /* prevents the file from being included twice. */

#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "./Hard_Sphere/Hard_Sphere.h" /* Hard Sphere header */
#include "./Hard_Disk/Hard_Disk.h" /* Hard Disk header */
#include "./Hard_Sphere_Square_Well/Hard_Sphere_Square_Well.h" /* Hard Sphere Square Well header */
#include "./Hard_Sphere_Double_Exp/Hard_Sphere_Double_Exp.h" /* Hard Sphere Square Double Yukawa */
#include "./Hard_Sphere_Double_Yukawa/Hard_Sphere_Double_Yukawa.h" /* Hard Sphere Double Yukawa header */

typedef struct structure_grid{
	gsl_vector * k;
	gsl_vector * kw;
	gsl_vector * S;
} structure_grid;

typedef struct liquid_params{
	double rho;
	double phi;
	double dim;
	double Tem;
	double * up;
	int nup;
} liquid_params;

void structure_grid_ini(structure_grid * Sg, int knp);

liquid_params
liquid_params_ini_phi(double phi, double dim, int ip);


void s_name_constructor(const char * sys, const char * approx, const char * extension, 
const int id, const liquid_params lp, char ** name );

double
s_function_selector_mono_sph( const char * sys, const char * approx, const char * fun, const
double k, const liquid_params lp );

void
s_grid_save_file( const structure_grid sg, const char * folder, 
const char * prefix, const char * suffix );

void
gsl_vector_s_function_selector_mono_sph(gsl_vector * sk, const char * sys, 
const char * approx, const char * fun, const gsl_vector * k, liquid_params lp );

void liquid_params_free(liquid_params * lp);

liquid_params liquid_params_unit_phi(double dim, int nup);

liquid_params liquid_params_sum(liquid_params lp1, liquid_params lp2);

liquid_params liquid_params_dif(liquid_params lp1, liquid_params lp2);

liquid_params liquid_params_scale(liquid_params lp, double scale);

liquid_params liquid_params_div(liquid_params lp1, liquid_params lp2);

double liquid_params_norm(liquid_params lp);

void structure_grid_free(structure_grid * Sg);

void structure_grid_memcpy(structure_grid dest, const structure_grid source);

double radial_distribution_3D(double r, double rho, const structure_grid Sg);

void gsl_vector_radial_distribution_3D(gsl_vector * g, gsl_vector * r, double rho, const structure_grid Sg);


#endif /* STRUCTURES_DOT_H */
