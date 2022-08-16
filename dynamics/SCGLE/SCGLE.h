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

typedef struct dyn_params{
	int st;
	const int it;
	double dtau;
	double kc;
	const double D0;
	const double tol;
}dyn_params;

typedef struct it_vars{
	gsl_vector * tau;
	gsl_vector * delz;
	gsl_vector * msd;
	gsl_vector * Dt;
	gsl_vector * dele;
	gsl_vector * tau_alpha;
	gsl_vector * lamk;
	gsl_matrix * Fc;
	gsl_matrix * Fs;
	gsl_vector * lamk_w;
	gsl_matrix * Fc_w;
	gsl_matrix * Fs_w;
} it_vars;

typedef struct save_dyn_vars{
	gsl_vector * tau;
	gsl_vector * delz;
	gsl_vector * dele;
	gsl_matrix * Fc;
	gsl_matrix * Fs;
	gsl_vector * k;
	gsl_vector * S;
	FILE * F_dyn;
	FILE * F_taua;
	FILE * F_Fc;
	FILE * F_Fs;
}save_dyn_vars;

typedef struct save_dyn_op{
	int save_delz;
	int save_deleta;
	int save_F;
	int save_gamma;
	int save_Dl;
	int write_delz;
	int write_taua;
	int write_deleta;
	int write_F;
}save_dyn_op;

typedef struct dyn_scalar{
	double Dl;
	double eta;
	double gamma;
	double Delz_inf;
}dyn_scalar;

void save_dyn_vars_free( save_dyn_vars * sdv );

double lambda_sph_mono( double k, double kc );

void lambda_sph_mono_gsl_v (gsl_vector ** lamk, const gsl_vector * k, const double kc );

dyn_params dyn_params_ini(int st, int it, double dtau, double kc, double D0 );

dyn_params dyn_params_ini_std();

dyn_params dyn_params_ini_HD();

save_dyn_op save_dyn_op_ini();

save_dyn_op no_save_dyn_op_ini();

void save_dyn_vars_ini( save_dyn_vars * sdv, const save_dyn_op so, const int knp, const char * folder, const char * prefix, const char * suffix );

dyn_scalar dyn_scalar_ini();

double gamma_sph_mono(structure_grid Sg, gsl_vector * lamk, liquid_params lp);

void dynamics_sph_mono( liquid_params lp, dyn_params dp, structure_grid Sg,
	save_dyn_vars * dyn, save_dyn_op op, dyn_scalar * ds );

liquid_params arrest_lp_in_limits(char * sys, char * approx, liquid_params lp0, liquid_params lp1, structure_grid Sg, gsl_vector * lamk, double tol);

void gsl_vector_fc_non_ergo_param_mono_sph(gsl_vector ** fc, structure_grid Sg, gsl_vector * lamk, double gamma);

void gsl_vector_fs_non_ergo_param_mono_sph(gsl_vector ** fs, structure_grid Sg, gsl_vector * lamk, double gamma);



#endif /* HARD_SPHERE_DOT_H */
