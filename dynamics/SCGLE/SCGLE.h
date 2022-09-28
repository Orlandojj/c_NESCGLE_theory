/* This header serves to import all structures headers defined in each
system folder. For the imported impl */
#ifndef SCGLE_DOT_H    /* This is an "include guard" */
#define SCGLE_DOT_H    /* prevents the file from being included twice. */

/* Lib dependencies */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "../../structures/structures.h"
#include "../../math/math_aux.h"

/* Data structure for the parameters 
needed for the computation of the dynamics */

typedef struct dynamics_parameters{
	int st;
	const int it;
	double dtau;
	double kc;
	const double D0;
	const double tol;
	const int minimum_decimations_number;
	const int maximum_decimations_number;
} dynamics_parameters;

/* Data structure for the computation of the 
medium times dynamics */
typedef struct intermediate_times_variables{
	/* Internal computation vectors */
	gsl_vector * tau;
	gsl_vector * delta_z;
	gsl_vector * lambda_k;
	/* Internal computation matrixes */
	gsl_matrix * Fc;
	gsl_matrix * Fs;
} intermediate_times_variables;

typedef struct dynamics_save_variables{
	/* Input Vectors */
	gsl_vector * k;
	gsl_vector * S;
	/* Internal computation vectors */
	gsl_vector * lambda_k;
	/* Output Vectors */
	gsl_vector * tau;
	gsl_vector * delta_zeta_t;
	gsl_vector * delta_eta_t;
	gsl_vector * msd;
	gsl_vector * Dt;
	gsl_vector * tau_alpha;
	/* Output Matrixes */
	gsl_matrix * Fc;
	gsl_matrix * Fs;
	/* Output Files */
	FILE * F_dyn;
	FILE * F_taua;
	FILE * F_Fc;
	FILE * F_Fs;
	/* Output Scalars */
	double Dl;
	double eta;
	double gamma;
	double Delz_inf;
} dynamics_save_variables;

typedef struct dynamics_save_options{
	/* if int = 1 -> save/write option activated */
	int delta_zeta;
	int delta_eta;
	int msd;
	int Dt;
	int Fc;
	int Fs;
	int gamma;
	int Dl;
	int tau_alpha;
} dynamics_save_options;

void free_dynamics_save_variables( dynamics_save_variables * dsv );

void close_dynamics_save_variables( dynamics_save_options dso, dynamics_save_variables * dsv);

void free_close_dynamics_save_variables( dynamics_save_options dso, dynamics_save_variables *dsv );

double lambda_spherical_mono( double k, double kc );

void gsl_vector_lambda_spherical_mono (gsl_vector ** lamk, const gsl_vector * k, const double kc );

dynamics_parameters dynamics_parameters_manual_ini(int st, int it, double dtau, double kc, double D0, int mindn,int maxdn );

dynamics_parameters dynamics_parameters_auto_ini();

dynamics_parameters dynamics_parameters_auto_HD_ini();

dynamics_save_options dynamics_save_options_auto_ini();

dynamics_save_options dynamics_save_options_no_save_ini();

int dynamics_save_options_sum_tau( dynamics_save_options dso );
int dynamics_save_options_sum_tau_only( dynamics_save_options dso );
int dynamics_save_options_sum_k( dynamics_save_options dso );

void dynamics_save_variables_spherical_mono_ini( dynamics_save_variables * dsv, const dynamics_save_options dso, gsl_vector * k, gsl_vector * Sk, const int knp, const dynamics_parameters dp, const char * folder, const char * prefix, const char * suffix );

double gamma_spherical_mono(structure_grid Sg, gsl_vector * lamk, liquid_params lp);

void dynamics_spherical_mono( liquid_params lp, dynamics_parameters dp, structure_grid Sg,
	dynamics_save_variables * dsv, dynamics_save_options dso);

liquid_params arrest_lp_in_limits(char * sys, char * approx, liquid_params lp0, liquid_params lp1, structure_grid Sg, gsl_vector * lamk, double tol);

void gsl_vector_fc_non_ergodicity_parameter_spherical_mono(gsl_vector ** fc, structure_grid Sg, gsl_vector * lamk, double gamma);

void gsl_vector_fs_non_ergodicity_parameter_spherical_mono(gsl_vector ** fs, structure_grid Sg, gsl_vector * lamk, double gamma);

void dynamics_mono_spherical_standard_defined_structures( liquid_params lp, char * sys, char * approx, char * folder );


#endif /* HARD_SPHERE_DOT_H */
