#include <stdio.h>
#include <math.h>
#include "./structures/structures.h"
#include "./dynamics/dynamics.h"
#include "./math/math_aux.h"
#include <time.h>


int main(void) {
  int i1,i2;
  double phi;
  clock_t start,end;
  start = clock();
  liquid_params lp=liquid_params_ini_phi(0.15,3,0);
  liquid_params final_lp=liquid_params_ini_phi(0.15,3,1);
  final_lp.up[0]=1.50;
  final_lp.Tem=1.50;
  printf("Program started\n");
  
  instant_change_dynamics_spherical_mono_standard_defined_structures( lp, final_lp, "HS", "HSSW", "VW", "SH", "./data/HSSW/" );
  /*(lp.up[0]=1.5;
  lp.Tem = 3.0;*/
  /*dynamics_mono_spherical_standard_defined_structures( lp, "HS", "VW", "./data/tests/" );*/
  /*final_lp.up[0] = 2.0;
  final_lp.up[1] = 1.0;
  final_lp.up[2] = 0.5;
  final_lp.Tem=0.80;}*/
  liquid_params_free(&final_lp);
  liquid_params_free(&lp);
  end=clock();
  double time_taken = ((double)(end-start))/ CLOCKS_PER_SEC;
  printf( "%s%f%s\n", "Program ended, time taken : ", time_taken, " sec");
  return 0;
}
