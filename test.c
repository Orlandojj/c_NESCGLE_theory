#include <stdio.h>
#include <math.h>
#include "./structures/structures.h"
#include "./dynamics/dynamics.h"
#include "./math/math_aux.h"
#include <time.h>

/* Function to test gamma quench */
void test_gamma_quench(){
  /* Liquid parameters variables */
  double phi=0.13, dim=3.0, Ti=4.0, Tf=0.5;
  char * sys="HSSW", * approx="SH", * fun="sk";
  double up[]={1.5};
  liquid_params lpi = liquid_params_ini_phi(phi, dim, 1);
  liquid_params lpf = liquid_params_ini_phi(phi, dim, 1);
  lpi.Tem = Ti; lpf.Tem = Tf;
  memcpy(lpi.up,up,sizeof(lpi.up)); 
  memcpy(lpf.up,up,sizeof(lpf.up));
  /* Integration and structure grids variables */
  const int knph1 = pow(2,8);
  const int knph2 = pow(2,8);
  const int knp = knph1+knph2;
  inst_change_vars icv = inst_change_vars_ini(knp, 0, lpi, lpf);
  gsl_integration_fixed_workspace * w1, *w2;
  const gsl_integration_fixed_type * int_T = gsl_integration_fixed_legendre;
  w1 = gsl_integration_fixed_alloc(int_T, knph1, 0.0, 10.0, 0.0, 0.0);
  w2 = gsl_integration_fixed_alloc(int_T, knph2, 10.0, 40.96, 0.0, 0.0);
  gsl_vectors_composed_integration_space(w1,w2, &icv.k, &icv.kw);
  gsl_vector_s_function_selector_mono_sph(icv.Si, sys, approx, fun,icv.k, icv.lpi );
  gsl_vector_s_function_selector_mono_sph(icv.Sf, sys, approx, fun,icv.k, icv.lpf );
  gsl_integration_fixed_free(w1); gsl_integration_fixed_free(w2);
  /* Defining dynamic variables */
  dyn_params dp = dyn_params_ini_std();
  double gamma_ua;
  double ua;
  int i1;
  gsl_vector * lamk = gsl_vector_alloc(knp); lambda_sph_mono_gsl_v( &lamk, icv.k, dp.kc );
  /* Computing gamma of ua and ua */
  inst_change_gamma_ua_sph_mono( icv, lamk, dp, &gamma_ua, &ua );
  printf("%s %1.9e %s %1.9e \n","gamma=",gamma_ua, " ua=", ua );
  /* Computing S(k;u=ua) */
  gsl_vector * aux_a=gsl_vector_alloc(knp);
  gsl_vector * sku=gsl_vector_alloc(knp);
  aux_ic_sph_mono(&aux_a,icv.k, icv.Sf, dp.D0);
  sku_inst_change_mono_sph(&sku,icv.Si, icv.Sf, aux_a, ua);
  FILE * F_test=fopen("test_k.dat", "w");
  for ( i1 = 0; i1 < knp; ++i1 ){
    fprintf(F_test,"%1.9e \t %1.9e \t %1.9e \t %1.9e \t %1.9e \n ", icv.k->data[i1], icv.Si->data[i1],
    icv.Sf->data[i1], sku->data[i1], lamk->data[i1] );
  }
  fclose(F_test);
  inst_change_vars_free(&icv);
  gsl_vector_free(aux_a); 
  gsl_vector_free(sku); 
  gsl_vector_free(lamk);
  liquid_params_free(&lpi);
  liquid_params_free(&lpf);
}

void dyn_test(){
  /* Liquid parameters */
  double phi=0.58;
  double T=1.0, lambda=1.5;
  double dim = 3.0;
  liquid_params lp = liquid_params_ini_phi(phi, dim, 0);
  /* Integration variables, structure grid definition and structure factor computation */
  const double k0=0.0; const double k1=10.0; const double k2=40.96;
  const int np  = pow(2,8);
  const int nph = np / 2;
  gsl_integration_fixed_workspace * w1, * w2;
  const gsl_integration_fixed_type * int_T = gsl_integration_fixed_legendre;
  w1 = gsl_integration_fixed_alloc(int_T, nph, k0, k1, 0.0, 0.0);
  w2 = gsl_integration_fixed_alloc(int_T, nph, k1, k2, 0.0, 0.0);
  structure_grid Sg={NULL,NULL,NULL}; structure_grid_ini(&Sg,np);
  gsl_vectors_composed_integration_space(w1,w2, &Sg.k, &Sg.kw);
  gsl_integration_fixed_free(w1); gsl_integration_fixed_free(w2);
  char * sys="HS";
  char * approx="VW";
  char * fun="sk";
  gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k, lp);
  /* Defining Dynamic variables */
  dyn_params dp = dyn_params_ini_std(); /* dyn_params_ini_HD(); */
  /* File handling variables */
  char * folder="./data/";
  char * fname = (char *)malloc(400*sizeof(char));
  s_name_constructor(sys,approx, "dat",1,lp, &fname );
  s_grid_save_file( Sg, "S_", "", fname );
  /* Initializing variables needed for the dynamics */
  printf("Initializing dynamic variable \n");
  double k_w[] = { 2.0, 4.0, 5.0, 2.0*M_PI, 7.18 }; const int nk_w = sizeof(k_w) / sizeof(k_w[0]);
  save_dyn_vars dyn; save_dyn_vars_ini(&dyn, nk_w, folder, "", fname );
  memcpy(dyn.k->data, k_w, sizeof(double)*nk_w);
  gsl_vector_s_function_selector_mono_sph(dyn.S, sys, approx, fun, dyn.k, lp); 
  save_dyn_op op = save_dyn_op_ini();
  dyn_scalar ds = dyn_scalar_ini();
  /* Computing the dynamics */
  printf("%s \n", "Initializing dynamics");
  dynamics_sph_mono( lp, dp, Sg, &dyn, op, &ds );
  printf("%s \t %1.9e \n", "Dl=", ds.Dl);
  printf("%s \t %1.9e \n", "gamma=", ds.gamma);
  /* Freeing memory */
  free(fname);
  structure_grid_free(&Sg);
  save_dyn_vars_free( &dyn );
  /*save_dyn_vars_free( &dyn );*/
  return;
}

void quench_test(){
  /* Initial state parameters */
  double dim = 3.0;
  double phi_i=0.3, T_i=2.0, lambda_i=1.5;
  liquid_params lp_i = liquid_params_ini_phi(phi_i, dim, 1); 
  lp_i.Tem=T_i; 
  lp_i.up[0]=lambda_i;
  char * sys_i="HSSW";
  char * approx_i="SH";
  /* Final state parameters */
  double phi_f=phi_i, T_f=0.1, lambda_f=lambda_i;
  liquid_params lp_f = liquid_params_ini_phi(phi_f, dim, 1); 
  lp_f.Tem=T_f; 
  lp_f.up[0]=lambda_f;
  char * sys_f=sys_i;
  char * approx_f=approx_i;
  /* Wave vector variables: number of points, segmentation, nodes, weights, integration variables and associated vectors */
  const int knp  = pow(2,10);
  const int knph = knp / 2;
  const double k0=0.0,k2=40.96,k1=k2/4.0; /* nodes limits */ 
  gsl_integration_fixed_workspace * w1, * w2;
  const gsl_integration_fixed_type * int_T = gsl_integration_fixed_legendre;
  w1 = gsl_integration_fixed_alloc(int_T, knph, k0, k1, 0.0, 0.0);
  w2 = gsl_integration_fixed_alloc(int_T, knph, k1, k2, 0.0, 0.0);
  double kwr[] = {2.0,4.0,5.0,2.0*M_PI,7.18};
  const int knwr = sizeof(kwr)/sizeof(kwr[0]);
  inst_change_vars icv = inst_change_vars_ini(knp, knwr, lp_i, lp_f);
  gsl_vectors_composed_integration_space(w1,w2, &icv.k, &icv.kw);\
  gsl_integration_fixed_free(w1); gsl_integration_fixed_free(w2);
  memcpy(icv.kwr->data,kwr,sizeof(kwr)); 
  /* Defining Dynamic variables */
  dyn_params dp = dyn_params_ini_std(); /* dyn_params_ini_HD(); */
  save_dyn_op op = save_dyn_op_ini();
  char * folder="./data/";
  char * prefix="";
  char * suffix="id_01.dat";
  char * fname= (char *)malloc(400*sizeof(char));
  /* Computing structure factor for initial and final state */
  gsl_vector_s_function_selector_mono_sph(icv.Si, sys_i, approx_i, "sk", icv.k, lp_i);
  gsl_vector_s_function_selector_mono_sph(icv.Sf, sys_f, approx_f, "sk", icv.k, lp_f);
  gsl_vector_s_function_selector_mono_sph(icv.Swri, sys_i, approx_i, "sk", icv.kwr, lp_i);
  gsl_vector_s_function_selector_mono_sph(icv.Swrf, sys_f, approx_f, "sk", icv.kwr, lp_f);
  s_name_constructor(sys_f,approx_f, "dat",20,lp_f, &fname );
  printf("name constructed \n");
  printf("%s \n", fname);
  printf("Initialization of quench \n");
  fflush(stdout);
  inst_change_mono_sph( icv, folder, fname, dp, op );
  free(fname);
  liquid_params_free(&lp_i);
  liquid_params_free(&lp_f);
  inst_change_vars_free(&icv);
  return ;
}

void HD_eq_dyn(){
  double phi=0.6, k, ck, sk_py, sk_vw;
  double T=1.0, lambda=1.5, sk_HSSW_vwsh;
  const int np  = pow(2,10);
  const int nph = np / 2;
  double * nodes, * weights;
  void * f_p;
  double a;
  gsl_integration_fixed_workspace * w1, * w2;
  const gsl_integration_fixed_type * int_T = gsl_integration_fixed_legendre;
  gsl_vector * v_n   = gsl_vector_alloc (nph);
  gsl_vector * v_w   = gsl_vector_alloc (nph);
  gsl_vector * v_F   = gsl_vector_alloc (nph);
  gsl_vector * v_dum = gsl_vector_alloc (nph);
  int i1,i2;
  FILE * sk_file;
  /* Defining Dynamic variables */
  double dim = 2.0;
  liquid_params lp = liquid_params_ini_phi(phi, dim, 0);
  dyn_params dp =  dyn_params_ini_HD();
  dp.kc=4.0;
  structure_grid Sg; structure_grid_ini(&Sg,np);
  char * prefix="";
  char * folder="./data/HD_ROTH_OLD2/";
  char * suffix="id_01.dat";
  printf("Initialization of save dynamics variables \n");
  char * fname= (char *)malloc(200*sizeof(char));
  char * sfname= (char *)malloc(200*sizeof(char));
  char * sys="HD";
  char * approx="ROTH";
  char * fun="sk";
  double sk;
  w1 = gsl_integration_fixed_alloc(int_T, nph, 0.0, 10.0, 0.0, 0.0);
  w2 = gsl_integration_fixed_alloc(int_T, nph, 10.0, 40.96, 0.0, 0.0);
  gsl_vectors_composed_integration_space(w1,w2, &Sg.k, &Sg.kw);
  FILE * scalar_dyns = fopen("./data/Dl.dat", "w");
  fprintf(scalar_dyns,"%s\n", "#|| 1 phi || 2 Dl(phi) ||#");
  for (i1=1; i1<73; ++i1){
    phi = i1*0.01; lp = liquid_params_ini_phi(phi, dim, 0);
    gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k, lp);
    printf("constructing name\n");
    s_name_constructor(sys,approx, "dat",0,lp, &fname );
    printf("name constructed \n");
    printf("%s \n", fname);
    s_grid_save_file( Sg, folder,prefix, fname );
    save_dyn_vars dyn; save_dyn_vars_ini(&dyn, 5, folder, prefix, fname); /*save_dyn_vars_ini( char * prefix, char * suffix );*/
    save_dyn_op op = save_dyn_op_ini();
    dyn_scalar ds = dyn_scalar_ini();
    printf("%s %1.9e %s %1.9e %s %1.9e \n","lp.phi=",lp.phi, "lp.rho=",lp.rho, "lp.dim=",lp.dim);
    printf("Initialization of dynamics \n");
    dyn.k->data[0] = 2.0;
    gsl_vector_set( dyn.S, 0, sk_hs_vw(phi, dyn.k->data[0] ) );
    dyn.k->data[1] = 4.0;
    gsl_vector_set( dyn.S, 1, sk_hs_vw(phi, dyn.k->data[1] ) );
    dyn.k->data[2] = 5.0;
    gsl_vector_set( dyn.S, 2, sk_hs_vw(phi, dyn.k->data[2] ) );
    dyn.k->data[3] = 2.0*M_PI;
    gsl_vector_set( dyn.S, 3, sk_hs_vw(phi, dyn.k->data[3] ) );
    dyn.k->data[4] = 7.18;
    gsl_vector_set( dyn.S, 4, sk_hs_vw(phi, dyn.k->data[4] ) );
    printf("%s %1.9e \n", "S(k=7.4)=",dyn.S->data[3]);
    printf("%s \n", "Initiating dynamics");
    dynamics_sph_mono( lp, dp, Sg, &dyn, op, &ds );
    printf("%s \t %1.9e \n", "Dl=", ds.Dl);
    printf("%s \t %1.9e \n", "gamma=", ds.gamma);
    fclose(dyn.F_dyn);
    fclose(dyn.F_Fc);
    fclose(dyn.F_Fs);
    fclose(dyn.F_taua);
    fprintf(scalar_dyns,"%1.9e \t %1.9e \n", lp.phi, ds.Dl);
    fflush(scalar_dyns);
  }
  fclose(scalar_dyns);
  return ;
}

void HD_arrest_finding(){
  /* Defining liquid parameters and structure parameters */
  double dim = 2.0;
  liquid_params lp_fluid = liquid_params_ini_phi(0.50, dim, 0);
  liquid_params lp_arrest = liquid_params_ini_phi(0.74, dim, 0);
  char * sys="HD";
  char * approx="ROTH";
  char * fun="sk";
  /* Defining dynamics parameters */
  dyn_params dp =  dyn_params_ini_HD(); /*dp.kc=4.0*M_PI;*/
  /* Defining structure grids through integration parameters */
  const int np  = pow(2,9);
  const int nph = np / 2;
  gsl_integration_fixed_workspace * w1, * w2;
  const gsl_integration_fixed_type * int_T = gsl_integration_fixed_legendre;
  structure_grid Sg,fdum;
  structure_grid_ini(&Sg,np);structure_grid_ini(&fdum,np);
  w1 = gsl_integration_fixed_alloc(int_T, nph, 0.0, 10.0, 0.0, 0.0);
  w2 = gsl_integration_fixed_alloc(int_T, nph, 10.0, 40.96, 0.0, 0.0);
  gsl_vectors_composed_integration_space(w1,w2, &Sg.k, &Sg.kw);
  gsl_integration_fixed_free(w1); gsl_integration_fixed_free(w2);
  structure_grid_memcpy(fdum,Sg);
  /* Auxiliary arrays for dynamics */
  gsl_vector * fc = gsl_vector_alloc(np);
  gsl_vector * fs = gsl_vector_alloc(np);
  gsl_vector * lamk = gsl_vector_alloc(np); lambda_sph_mono_gsl_v( &lamk, Sg.k, dp.kc );
  /* Computing the arrest liquid parameters */
  lp_arrest=arrest_lp_in_limits(sys,approx, lp_fluid, lp_arrest, Sg, lamk, 1E-10);
  printf("%s \t %1.9e \n", "The arrest area fraction is", lp_arrest.phi);
  /* Computing the structure factor of arrested state */
  gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k,  lp_arrest);
  double gamma = gamma_sph_mono( Sg, lamk, lp_arrest);
  printf("%s \t %1.9e \n", "gamma=", gamma);
  gsl_vector_fs_non_ergo_param_mono_sph(&fs, Sg, lamk, gamma);
  printf("fs computed \n");
  gsl_vector_fc_non_ergo_param_mono_sph(&fc, Sg, lamk, gamma);
  printf("fc computed \n");
  fflush(stdout);
  char * folder = "./data/";
  char * prefix = "fc_";
  char * fname = (char *)malloc(200*sizeof(char));
  s_name_constructor(sys,approx, "dat",1,lp_arrest, &fname );
  prefix="S_";
  s_grid_save_file( Sg, folder, prefix, fname );
  gsl_vector_memcpy(fdum.S,fc);
  printf("%s\n",fname);
  prefix="fc_";
  s_grid_save_file( fdum, folder, prefix, fname );
  prefix="fs_";
  gsl_vector_memcpy(fdum.S,fs);
  s_grid_save_file( fdum, folder, prefix, fname );
  free(fname);
  gsl_vector_free(lamk);
  gsl_vector_free(fc); 
  gsl_vector_free(fs); 
  structure_grid_free(&Sg);
  structure_grid_free(&fdum);
  return ;
}

void HD_eq_dyn_near_arrest(){
  const int np  = pow(2,10);
  const int nph = np / 2;
  int i1,i2;
  /* Defining Dynamic variables */
  double dim = 2.0;
  liquid_params lp_fluid = liquid_params_ini_phi(0.50, dim, 0);
  liquid_params lp_arrest = liquid_params_ini_phi(0.74, dim, 0);
  liquid_params lp_writte = liquid_params_ini_phi(0.50, dim, 0);
  dyn_params dp =  dyn_params_ini_HD(); /*dp.kc=1000.0*M_PI;*/
  structure_grid Sg; structure_grid_ini(&Sg,np);
  char * prefix="";
  char * folder="./data/";
  char * suffix="id_01.dat";
  char * fname= (char *)malloc(200*sizeof(char));
  char * sys="HD";
  char * approx="ROTH";
  char * fun="sk";
  /* setting integration workspace */
  gsl_integration_fixed_workspace * w1, * w2;
  const gsl_integration_fixed_type * int_T = gsl_integration_fixed_legendre;
  w1 = gsl_integration_fixed_alloc(int_T, nph, 0.0, 10.0, 0.0, 0.0);
  w2 = gsl_integration_fixed_alloc(int_T, nph, 10.0, 40.96, 0.0, 0.0);
  gsl_vectors_composed_integration_space(w1,w2, &Sg.k, &Sg.kw);
  gsl_vector * lamk = gsl_vector_alloc(np); lambda_sph_mono_gsl_v( &lamk, Sg.k, dp.kc );
  FILE * lam_file_test = fopen("lam.dat","w");
  for (i1=0; i1<Sg.k->size;i1++){
    fprintf(lam_file_test,"%1.9e\t%1.9e\n",Sg.k->data[i1],lamk->data[i1]);
  }
  fclose(lam_file_test);
  lp_arrest=arrest_lp_in_limits(sys,approx, lp_fluid, lp_arrest, Sg, lamk, 1E-10);
  printf("%s \t %1.9e \n", "The arrest area fraction is", lp_arrest.phi);
  printf("Computing dynamics relative to arrest state \n");
  /* Computing dynamics relative to dynamically arrested states, advancing logarithmically in *(10^0.1) 
  intervals starting from a 10^(-6) separation finishing in a 10^(-1) separation */
  save_dyn_vars dyn;
  save_dyn_op op = save_dyn_op_ini();
  dyn_scalar ds = dyn_scalar_ini();
  double interval=pow(10.0,0.1);
  double dphi = 1e-6;
  strcpy(fname,folder); strcat(fname,"Dl.dat");
  FILE * Dl_file=fopen(fname,"w");
  fprintf(Dl_file,"#|| ϕ || ϕ-ϕₐ || Dl # \n" );
  for (i1=1; i1<=59; i1++){
    lp_writte = liquid_params_ini_phi(0.01*i1, dim, 0);
    printf("constructing name\n");
    s_name_constructor(sys,approx, "dat",1,lp_writte, &fname );
    printf("name constructed \n");
    printf("%s \n", fname);
    lp_fluid = liquid_params_ini_phi(lp_arrest.phi - dphi, dim, 0);
    dphi = dphi*interval;
    gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k,  lp_fluid);
    s_grid_save_file( Sg, folder,prefix, fname );
    save_dyn_vars_ini(&dyn, 5, folder, prefix, fname );
    dyn.k->data[0]=2.0; dyn.k->data[1]=4.0; dyn.k->data[2]=5.0; dyn.k->data[3]=2.0*M_PI; dyn.k->data[4]=7.18;
    gsl_vector_s_function_selector_mono_sph(dyn.S, sys, approx, fun, dyn.k, lp_fluid);
    dynamics_sph_mono( lp_fluid, dp, Sg, &dyn, op, &ds );
    fprintf(Dl_file,"%1.9e\t%1.9e\t%1.9e \n",lp_fluid.phi,lp_arrest.phi-lp_fluid.phi,ds.Dl);
    fflush(Dl_file);
  }
  fclose(Dl_file);
  return ;
}

int main(void) {
  clock_t start,end;
  start = clock();
  printf("Program started\n");
  /*dyn_test();*/
  /*HD_arrest_finding();*/
  /*HD_eq_dyn_near_arrest();*/
  /*HD_eq_dyn();*/
  quench_test();
  /*test_gamma_quench();*/
  end=clock();
  double time_taken = ((double)(end-start))/ CLOCKS_PER_SEC;
  printf( "%s%f%s\n", "Program ended, time taken : ", time_taken, " sec");
  return 0;
}
