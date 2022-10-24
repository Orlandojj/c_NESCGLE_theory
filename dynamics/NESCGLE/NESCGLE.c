#include "NESCGLE.h"

instant_change_variables
instant_change_variables_ini(int knp, int knwr,  liquid_params lpi, liquid_params lpf){
  instant_change_variables icv;
  int i1;
  icv.k=gsl_vector_alloc(knp);
  icv.kw=gsl_vector_alloc(knp);
  icv.Si=gsl_vector_alloc(knp);
  icv.Sf=gsl_vector_alloc(knp);
  icv.kwr=NULL;
  icv.Swri=NULL;
  icv.Swrf=NULL;
  if (knwr>0) {
    icv.kwr=gsl_vector_alloc(knwr);
    icv.Swri=gsl_vector_alloc(knwr);
    icv.Swrf=gsl_vector_alloc(knwr);
  }
  icv.lpi=liquid_params_ini_phi(lpi.phi, lpi.dim, lpi.nup);
  icv.lpf=liquid_params_ini_phi(lpf.phi, lpf.dim, lpf.nup);
  icv.lpi.Tem=lpi.Tem;icv.lpf.Tem=lpf.Tem;
  if (lpi.nup>0){for (i1=0; i1<lpi.nup; ++i1){icv.lpi.up[i1] = lpi.up[i1];}}
  if (lpf.nup>0){for (i1=0; i1<lpf.nup; ++i1){icv.lpf.up[i1] = lpf.up[i1];}}
  return icv;
}

void instant_change_variables_free( instant_change_variables * icv ){
  if ( icv->k    != NULL ) { gsl_vector_free(icv->k);    }
  if ( icv->kw   != NULL ) { gsl_vector_free(icv->kw);   }
  if ( icv->Si   != NULL ) { gsl_vector_free(icv->Si);   }
  if ( icv->Sf   != NULL ) { gsl_vector_free(icv->Sf);   }
  if ( icv->kwr  != NULL ) { gsl_vector_free(icv->kwr);  }
  if ( icv->Swri != NULL ) { gsl_vector_free(icv->Swri); }
  if ( icv->Swrf != NULL ) { gsl_vector_free(icv->Swrf); }
  liquid_params_free(&icv->lpi); 
  liquid_params_free(&icv->lpf);
}

void aux_ic_sph_mono(gsl_vector ** aux_a, gsl_vector * k, gsl_vector * Skf, double D0 ){
  gsl_vector_memcpy(*aux_a,k);
  gsl_vector_mul(*aux_a,*aux_a);
  gsl_vector_div(*aux_a,Skf);
  gsl_vector_scale(*aux_a, 2.0 * D0 );
  return;
}

void
sku_inst_change_mono_sph(gsl_vector ** Sku, const gsl_vector * Ski, const gsl_vector * Skf, const gsl_vector * aux_a, const double u){
  int i1;
  int knp=Ski->size;
  double dum1,dum2;
  for (i1=0; i1 < knp; ++i1){
    dum1 = aux_a->data[i1]*u;
    if (dum1 < 500.0 ){
      gsl_vector_set(*Sku,i1,Skf->data[i1] + ( Ski->data[i1] - Skf->data[i1] ) * exp(-dum1));
    }
    else {
      gsl_vector_set(*Sku,i1,Skf->data[i1]);
    }
  }
  return;
}

int
assert_positive_Sk(gsl_vector * Sk){
  int assertion=1;
   if ( 0 >= gsl_vector_min(Sk) ) { assertion=0; }
   return assertion;
}

void inst_change_gamma_ua_sph_mono(instant_change_variables icv, gsl_vector * lamk, dynamics_parameters dp, double * gamma_ua, double * ua )
{
  double gamma_i, gamma_f,gamma_t;
  int assert_v,convergence;
  int knp = icv.k->size;
  structure_grid Sg = {NULL,NULL,NULL}; structure_grid_ini(&Sg,knp);
  gsl_vector * aux_a=gsl_vector_alloc(knp);
  double u_max,u_min,u_test,error;
  const double tol = 1E-10;
  gsl_vector_memcpy(Sg.k,icv.k);
  gsl_vector_memcpy(Sg.kw,icv.kw);
  /* Checking bounds */
  gsl_vector_memcpy(Sg.S,icv.Si);
  assert_v=assert_positive_Sk(Sg.S);
  if(assert_v==1){
    gamma_i = gamma_spherical_mono( Sg, lamk, icv.lpi);
    if (gamma_i<1E40){
      printf("Warning, initiating from an arrested state, returning u value to liquid if posible \n");}
  }
  else if (assert_v==0){
    gamma_i = 1E99;
  }
  else{
    structure_grid_free(&Sg); gsl_vector_free(aux_a);
    printf("Not a valid assertion for S(k)");
    exit(1);
  }
  aux_ic_sph_mono(&aux_a,icv.k, icv.Sf, dp.D0);
  u_min=0.0;
  u_max=1.0;
  convergence=0;
  while(convergence==0){
    sku_inst_change_mono_sph(&Sg.S,icv.Si, icv.Sf, aux_a, u_max);
    gamma_f = gamma_spherical_mono( Sg, lamk, icv.lpf);
    if ( gamma_f < 1E40 ){convergence = 1;}
    else if (u_max < 1E40){ u_min=u_max; u_max *= 2.0;}
    else {ua[0]=1E40;gamma_ua[0]=gamma_f;
    structure_grid_free(&Sg); gsl_vector_free(aux_a);return;}
  }
  convergence=0;
  while(convergence==0){
    u_test = ( u_max + u_min ) * 0.5;
    sku_inst_change_mono_sph(&Sg.S, icv.Si, icv.Sf, aux_a, u_test);
    gamma_t = gamma_spherical_mono( Sg, lamk, icv.lpf);
    if ( gamma_t < 1E40 ){ u_max = u_test;  gamma_f = gamma_t; }
    else { u_min = u_test; gamma_i = gamma_t; }
    error = (u_max-u_min) / u_max;
    if (error<tol){convergence=1; ua[0] = u_max; gamma_ua[0] = gamma_f; 
    structure_grid_free(&Sg); gsl_vector_free(aux_a); return;}
  }
  structure_grid_free(&Sg); gsl_vector_free(aux_a); return;
}

void instant_change_dynamics_spherical_mono( instant_change_variables icv, char * folder, char * fname, dynamics_parameters dp, dynamics_save_options op, int write_S ) {
  int knp=icv.k->size;
  int knpwr=icv.Swri->size;
  const double Dl_ratio_smoothness = 1.05; /* Value greater than 1, recommended lower than 1.10 */
  structure_grid Sg = {NULL,NULL,NULL}; 
  structure_grid_ini(&Sg,knp);  
  gsl_vector_memcpy(Sg.k,icv.k); 
  gsl_vector_memcpy(Sg.kw,icv.kw);
  structure_grid Sgw = {NULL, NULL, NULL};
  structure_grid_ini(&Sgw,knpwr);
  double Dl_initial, Dl_final, Dl_no_save;
  double u, t, ua, gamma_ua;
  gsl_vector * aux_a = gsl_vector_alloc(knp);aux_ic_sph_mono( &aux_a,icv.k, icv.Sf, dp.D0);
  gsl_vector * lamk = gsl_vector_alloc(knp); 
  gsl_vector_lambda_spherical_mono(&lamk, icv.k, dp.kc );
  gsl_vector * aux_a_w=gsl_vector_alloc(knpwr); aux_ic_sph_mono( &aux_a_w, icv.kwr, icv.Swrf, dp.D0);
  dynamics_save_variables dyn_save;
  dynamics_save_options no_save=dynamics_save_options_no_save_ini();
  int i1, convergence;
  double du_dt_1, du_dt_2, du, dt, t_pre, Dl_pre, t_save, u_save, du_save;
  double e_Dl;
  char * u_char=(char *)malloc(20*sizeof(char));
  char * u_char_S=(char *)malloc(20*sizeof(char));
  char * u_char_g=(char *)malloc(20*sizeof(char));
  /* Minimum mobility tolerance and file writting scale in terms of waiting time powers */
  const double tol_Dl=dp.tol, t_save_scale= pow(10.0,0.25);
  int write_dyn = 0;
  /* Radial distribution function grid initialization */
  int rnp = 2000;
  double rmax = 20.0;
  structure_grid gr = {NULL,NULL,NULL}; structure_grid_ini(&gr,rnp);
  for (i1=0; i1<rnp; i1++){ gr.k->data[i1] = rmax * ( (double) i1 / (double) rnp) ; gr.kw->data[i1] = 0.0;}
  if ( dynamics_save_options_sum_tau_only(op) > 0 ) { write_dyn = 1; }
  /* Save format and variables for mobility files */
  printf("Opening files\n");
  char * t_file_name = (char *)malloc(400*sizeof(char));
  strcpy(t_file_name,folder); 
  strcat(t_file_name,"b_"); 
  strcat(t_file_name, fname);
  FILE * t_file=fopen(t_file_name,"w");
  strcpy(t_file_name,folder); 
  strcat(t_file_name,"t_save_"); 
  strcat(t_file_name, fname);
  FILE * t_save_file;
  fprintf(t_file,"#1t,2b(t),3u(t),4uₐ-u(t)\n" );
  if (write_dyn == 1) {
    t_save_file=fopen(t_file_name,"w");
    fprintf(t_save_file,"#1quench_id,2t,3u(t),4b(t)\n" );
  }
  free(t_file_name);
  int i_save=0;
  /* Computing limit values of ua */
  printf("computing ua\n");
  inst_change_gamma_ua_sph_mono( icv, lamk, dp, &gamma_ua, &ua );
  gsl_vector_free(lamk);
  printf("%s \t %1.9e\n","ua=",ua);
  /* Computing initial dynamics */
  printf("computing initial dynamics\n");
  gsl_vector_memcpy(Sg.S,icv.Si); 
  sprintf(u_char,"%s%03d%s","u_",0,"_");
  sprintf(u_char_S,"%s%03d%s","u_",0,"_S_");
  sprintf(u_char_g,"%s%03d%s","u_",0,"_g_");
  dynamics_save_variables_spherical_mono_ini(&dyn_save, op,icv.kwr, icv.Swri, icv.k->size, dp, folder, u_char, fname); 
  gsl_vector_memcpy(dyn_save.k,icv.kwr);
  gsl_vector_memcpy(dyn_save.S,icv.Swri);
  gsl_vector_lambda_spherical_mono(&dyn_save.lambda_k,dyn_save.k,dp.kc); 
  if ( write_S == 1 ) {
    s_grid_save_file( Sg, folder, u_char_S, fname ); 
    gsl_vector_radial_distribution_3D(gr.S, gr.k, icv.lpi.rho, Sg); 
    s_grid_save_file( gr, folder, u_char_g, fname );
  }
  dynamics_spherical_mono( icv.lpi, dp, Sg, &dyn_save, op);
  Dl_initial = dyn_save.Dl;
  free_close_dynamics_save_variables(op,&dyn_save);
  /* Computing final state dynamics */
  printf("computing final dynamics\n");
  sprintf(u_char,"%s%03d%s","u_",999,"_");
  sprintf(u_char_S,"%s%03d%s","u_",999,"_S_"); 
  sprintf(u_char_g,"%s%03d%s","u_",999,"_g_");
  dynamics_save_variables_spherical_mono_ini(&dyn_save, op, icv.kwr, icv.Swrf, icv.k->size, dp, folder, u_char, fname);
  gsl_vector_memcpy(dyn_save.k,icv.kwr);
  if ( ua < 1E40){
    sku_inst_change_mono_sph(&Sg.S,icv.Si, icv.Sf, aux_a, ua);
    sku_inst_change_mono_sph(&dyn_save.S, icv.Swri, icv.Swrf, aux_a_w, ua);
    gsl_vector_lambda_spherical_mono(&dyn_save.lambda_k,dyn_save.k,dp.kc); 
    dynamics_spherical_mono( icv.lpf, dp, Sg, &dyn_save, op );
  }
  else{
    gsl_vector_memcpy(Sg.S,icv.Sf);
    gsl_vector_memcpy(dyn_save.S,icv.Swrf);
    gsl_vector_lambda_spherical_mono(&dyn_save.lambda_k,dyn_save.k,dp.kc); 
    dynamics_spherical_mono( icv.lpf, dp, Sg, &dyn_save, op );
  }
  if ( write_S == 1 ) {
    s_grid_save_file( Sg, folder, u_char_S, fname ); 
    gsl_vector_radial_distribution_3D(gr.S, gr.k, icv.lpi.rho, Sg); 
    s_grid_save_file( gr, folder, u_char_g, fname );
    }
  Dl_final = dyn_save.Dl;
  free_close_dynamics_save_variables(op,&dyn_save);
  /* Variables initialization for the waiting times loop */
  convergence = 0;
  du      = 1E-7 * Dl_initial; /* Initial u differential with dt≈du/Dl; where Dl is the long time diffusion coefficient of the initial state */
  t       = 0.0;
  u       = 0.0;
  t_save  = 1E-5;
  du_dt_1 = 0.0;
  i_save  = 0;
  if (write_dyn == 1) {fprintf( t_save_file,"%d%s%1.9e%s%1.9e%s%1.9e \n", i_save,",", 0.0,",", 0.0,",", Dl_initial );}
  fprintf( t_file,"%1.9e%s%1.9e%s%1.9e%s%1.9e \n", 0.0,",", Dl_initial, ",",0.0,",", ua );
  Dl_pre=Dl_initial;
  while ( convergence == 0 ){
    /* Computing the dynamics for the next u-value */
    u+=du;
    sku_inst_change_mono_sph(&Sg.S, icv.Si, icv.Sf, aux_a, u);
    dynamics_spherical_mono( icv.lpf, dp, Sg, &dyn_save, no_save );
    Dl_no_save = dyn_save.Dl;
    t_pre=t;
    dt = du / Dl_no_save;
    t += dt;
    
    /* Computing the decision to end the time loop */
    e_Dl = fabs(1.0 - ( Dl_no_save / Dl_final )); /* Error between current Dl and final Dl */
    if ( e_Dl < tol_Dl || u >= ua || Dl_no_save <= tol_Dl ){convergence=1;}
    if ( Dl_no_save > tol_Dl ) {
      printf("%s %1.9e \t %s %1.9e \t %s %1.9e\n","Dl rel error w final value: ",e_Dl,"u: ", u,"t: ",t);
      fprintf( t_file,"%1.9e%s%1.9e%s%1.9e%s%1.9e \n", t,",", Dl_no_save,",", u,",", ua-u );
      fflush(t_file);
      /* Computing the dynamics for the time-save value */
      
      while ( t_pre <= t_save && t > t_save && write_dyn == 1 ){
        i_save +=1;
        sprintf(u_char,"%s%03d%s","u_",i_save,"_");
        sprintf(u_char_S,"%s%03d%s","u_",i_save,"_S_");
        sprintf(u_char_g,"%s%03d%s","u_",i_save,"_g_"); 
        du_save = ( t_save - t_pre ) * 0.5 * ( Dl_no_save + Dl_pre );
        u_save  = u-du + du_save;
        sku_inst_change_mono_sph(&Sg.S, icv.Si, icv.Sf, aux_a, u_save);
        sku_inst_change_mono_sph(&Sgw.S, icv.Swri, icv.Swrf, aux_a_w, u_save);
        dynamics_save_variables_spherical_mono_ini( &dyn_save, op, Sgw.k, Sgw.S, Sg.k->size, dp, folder, u_char, fname);
        gsl_vector_memcpy(dyn_save.k,icv.kwr);
        s_grid_save_file( Sg, folder, u_char_S, fname );
        gsl_vector_radial_distribution_3D(gr.S, gr.k, icv.lpi.rho, Sg); 
        s_grid_save_file( gr, folder, u_char_g, fname );
        gsl_vector_lambda_spherical_mono(&dyn_save.lambda_k,dyn_save.k,dp.kc); 
        dynamics_spherical_mono( icv.lpf, dp, Sg, &dyn_save, op );
        fprintf(t_save_file,"%d%s%1.9e%s%1.9e%s%1.9e \n", i_save,",", t_save,",", u_save,",", dyn_save.Dl );
        fflush(t_save_file);
        t_save *= t_save_scale;
        free_close_dynamics_save_variables(op, &dyn_save);
      }
    /* Computing the change in du in terms of change in (du/dt)  */
    if( Dl_pre / Dl_no_save > Dl_ratio_smoothness ) {du *= 0.05 / ( ( Dl_pre / Dl_no_save ) - 1.0) ;}
    else if( Dl_no_save / Dl_pre > Dl_ratio_smoothness ) {du *= 0.05 / ( ( Dl_no_save / Dl_pre ) - 1.0);}
    else if ( u + 2.0*du < ua && 2.0 * du < 0.05 * ua ) { du *= 2.0; }
    }
    Dl_pre = Dl_no_save;
    
  }
  fprintf(t_file,"%1.9e%s%1.9e%s%1.9e%s%1.9e \n", t,",", Dl_final,",", ua,",", 0.0 );
  if (write_dyn == 1) {fclose(t_save_file);}
  fclose(t_file);
  gsl_vector_free(aux_a);
  gsl_vector_free(aux_a_w);
  free(u_char);
  free(u_char_S);
  free(u_char_g);
  structure_grid_free(&Sg);
  structure_grid_free(&Sgw);
  structure_grid_free(&gr);
}

void instant_change_dynamics_spherical_mono_standard_defined_structures( liquid_params initial_lp, liquid_params final_lp, char * initial_sys, char * final_sys, char * initial_approx, char * final_approx, char * folder ) {
  /* Integration variables, structure grid definition and structure factor computation */
  const double k0=0.0; const double k1=10.0; const double k2=40.96;
	const int np  = pow(2,8);
	const int nph = np / 2;
	char * fun="sk";
  /* Computing integration nodes and weights */
	gsl_integration_fixed_workspace * w1, * w2;
	const gsl_integration_fixed_type * int_T = gsl_integration_fixed_legendre;
	w1 = gsl_integration_fixed_alloc(int_T, nph, k0, k1, 0.0, 0.0);
	w2 = gsl_integration_fixed_alloc(int_T, nph, k1, k2, 0.0, 0.0);
  /* Initializing and computing the instant change variables needed for the computation of instant change dynamics */
  const int knwr = 5;
	instant_change_variables icv = instant_change_variables_ini(np,knwr,initial_lp,final_lp);
  gsl_vectors_composed_integration_space(w1,w2, &icv.k, &icv.kw);
	gsl_integration_fixed_free(w1); 
  gsl_integration_fixed_free(w2);
	gsl_vector_s_function_selector_mono_sph(icv.Si, initial_sys, initial_approx, fun, icv.k, initial_lp);
  gsl_vector_s_function_selector_mono_sph(icv.Sf, final_sys, final_approx, fun, icv.k, final_lp);
  icv.kwr->data[0] = 2.0;
  icv.kwr->data[1] = 4.0;
  icv.kwr->data[2] = 5.0;
  icv.kwr->data[3] = 2.0 * M_PI;
  icv.kwr->data[4] = 7.18;
  gsl_vector_s_function_selector_mono_sph(icv.Swri, initial_sys, initial_approx, fun, icv.kwr, initial_lp);
  gsl_vector_s_function_selector_mono_sph(icv.Swrf, final_sys, final_approx, fun, icv.kwr, final_lp);
  /* Setting the name structure for save files */
  char * fname = (char *)malloc(400*sizeof(char));
	s_name_constructor( final_sys, final_approx, "csv", 1, final_lp, &fname );
  /* Setting dynamics parameters and dynamics save options */
  dynamics_parameters dp = dynamics_parameters_auto_ini();
  dynamics_save_options dso = dynamics_save_options_auto_ini();
	/* Initializing dynamics parameters and dynamics save options */
  instant_change_dynamics_spherical_mono( icv, folder, fname, dp, dso, 1 );
  instant_change_variables_free( &icv );
  free(fname);
  return;
}