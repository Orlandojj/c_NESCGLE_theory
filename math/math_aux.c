#include "math_aux.h"


/* Joins two integration workspaces into a gsl_vector adn frees the integration workspaces memory */
void
gsl_vectors_composed_integration_space(const gsl_integration_fixed_workspace * w1,
const gsl_integration_fixed_workspace * w2, gsl_vector ** nodes, gsl_vector ** weights){
  /* Copying data into the first elements */
  memcpy(nodes[0]->data,w1->x,sizeof(double)*w1->n);
  memcpy(weights[0]->data,w1->weights,sizeof(double)*w1->n);
  /* Copying data starting at w1->n */
  memcpy(&nodes[0]->data[w1->n],w2->x,sizeof(double)*w2->n);
  memcpy(&weights[0]->data[w1->n],w2->weights,sizeof(double)*w2->n);
  return;
}
