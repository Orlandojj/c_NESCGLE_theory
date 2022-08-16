/* This header serves to import all structures headers defined in each
system folder. For the imported impl */
#ifndef NESCGLE_DOT_H    /* This is an "include guard" */
#define NESCGLE_DOT_H    /* prevents the file from being included twice. */

#include "../SCGLE/SCGLE.h"
#include <math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "../../structures/structures.h"

typedef struct inst_change_vars{
  gsl_vector * k;
  gsl_vector * kw;
  gsl_vector * Si;
  gsl_vector * Sf;
  gsl_vector * kwr;
  gsl_vector * Swri;
  gsl_vector * Swrf;
  liquid_params lpi;
  liquid_params lpf;
}inst_change_vars;

inst_change_vars inst_change_vars_ini(int knp, int knwr, liquid_params lpi, liquid_params lpf);
void inst_change_vars_free(inst_change_vars * icv);
void aux_ic_sph_mono(gsl_vector ** aux_a, gsl_vector * k, gsl_vector * Skf, double D0 );
void sku_inst_change_mono_sph(gsl_vector ** Sku, const gsl_vector * Ski, const gsl_vector * Skf, const gsl_vector * aux_a, const double u);
void inst_change_gamma_ua_sph_mono(inst_change_vars icv, gsl_vector * lamk, dyn_params dp, double * gamma_ua, double * ua );
void inst_change_mono_sph( inst_change_vars icv, char * folder, char * fname, dyn_params dp, save_dyn_op op, int write_S );
/*
void
inst_change_mono_sph( inst_change_vars icv, dyn_params dp, save_dyn_vars * dyn, save_dyn_op op )
*/
#endif /* NESCGLE_DOT_H */
