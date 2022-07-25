/* This header serves to import all structures headers defined in each
system folder. For the imported impl */
#ifndef MATH_AUX_DOT_H    /* This is an "include guard" */
#define MATH_AUX_DOT_H    /* prevents the file from being included twice. */
#include <gsl/gsl_integration.h>
#include <gsl/gsl_vector.h>
#include <string.h>

void gsl_vectors_composed_integration_space(const gsl_integration_fixed_workspace * w1,
const gsl_integration_fixed_workspace * w2, gsl_vector ** nodes, gsl_vector ** weights);

#endif /*  */
