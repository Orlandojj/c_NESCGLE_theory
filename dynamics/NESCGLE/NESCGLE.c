#include "NESCGLE.h"

inst_change_vars
inst_change_vars_ini(int knp, int knwr,  liquid_params lpi, liquid_params lpf){
  inst_change_vars icv;
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
  if (lpi.nup>0){memcpy(icv.lpi.up, lpi.up,sizeof(lpi.up));} 
  if (lpf.nup>0){memcpy(icv.lpf.up, lpf.up,sizeof(lpf.up));}
  return icv;
}

void inst_change_vars_free(inst_change_vars * icv){
  gsl_vector_free(icv->k);
  gsl_vector_free(icv->kw);
  gsl_vector_free(icv->Si);
  gsl_vector_free(icv->Sf);
  if ( icv->kwr != NULL ) {gsl_vector_free(icv->kwr);}
  if ( icv->Swri != NULL ) {gsl_vector_free(icv->Swri);}
  if ( icv->Swrf != NULL ) {gsl_vector_free(icv->Swrf);}
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

void
inst_change_gamma_ua_sph_mono(inst_change_vars icv,
gsl_vector * lamk, dyn_params dp, double * gamma_ua, double * ua )
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
    gamma_i = gamma_sph_mono( Sg, lamk, icv.lpi);
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
    gamma_f = gamma_sph_mono( Sg, lamk, icv.lpf);
    if ( gamma_f < 1E40 ){convergence = 1;}
    else if (u_max < 1E40){ u_min=u_max; u_max *= 2.0;}
    else {ua[0]=1E40;gamma_ua[0]=gamma_f;
    structure_grid_free(&Sg); gsl_vector_free(aux_a);return;}
  }
  convergence=0;
  while(convergence==0){
    u_test = ( u_max + u_min ) * 0.5;
    sku_inst_change_mono_sph(&Sg.S, icv.Si, icv.Sf, aux_a, u_test);
    gamma_t = gamma_sph_mono( Sg, lamk, icv.lpf);
    if ( gamma_t < 1E40 ){ u_max = u_test;  gamma_f = gamma_t; }
    else { u_min = u_test; gamma_i = gamma_t; }
    error = (u_max-u_min) / u_max;
    if (error<tol){convergence=1; ua[0] = u_max; gamma_ua[0] = gamma_f; 
    structure_grid_free(&Sg); gsl_vector_free(aux_a); return;}
  }
  structure_grid_free(&Sg); gsl_vector_free(aux_a); return;
}

void
inst_change_mono_sph( inst_change_vars icv, char * folder, char * fname, dyn_params dp, save_dyn_op op, int write_S ) {
  int knp=icv.k->size;
  int knpwr=icv.Swri->size;
  structure_grid Sg = {NULL,NULL,NULL}; structure_grid_ini(&Sg,knp);  
  gsl_vector_memcpy(Sg.k,icv.k); gsl_vector_memcpy(Sg.kw,icv.kw);
 
  dyn_scalar ds = dyn_scalar_ini();
  dyn_scalar dsf = dyn_scalar_ini();
  dyn_scalar ds_save = dyn_scalar_ini();
  double u, t, ua, gamma_ua;
  gsl_vector * aux_a = gsl_vector_alloc(knp);aux_ic_sph_mono( &aux_a,icv.k, icv.Sf, dp.D0);
  gsl_vector * lamk = gsl_vector_alloc(knp); lambda_sph_mono_gsl_v(&lamk, icv.k, dp.kc );
  gsl_vector * aux_a_w=gsl_vector_alloc(knpwr); aux_ic_sph_mono( &aux_a_w, icv.kwr, icv.Swrf, dp.D0);
  save_dyn_vars dyn_save;
  save_dyn_op no_save=no_save_dyn_op_ini();
  int i1, convergence;
  double du_dt_1, du_dt_2, du, dt, t_pre, Dl_pre, t_save, u_save, du_save;
  double e_Dl;
  char * u_char=(char *)malloc(20*sizeof(char));
  char * u_char_S=(char *)malloc(20*sizeof(char));
  char * u_char_g=(char *)malloc(20*sizeof(char));
  /* Minimum mobility tolerance and file writting scale in terms of waiting time powers */
  const double tol_Dl=dp.tol, t_save_scale= pow(10.0,0.1);
  int write_dyn = 0;
  /* Radial distribution function grid initialization */
  int rnp = 1000;
  double rmax = 10.0;
  structure_grid gr = {NULL,NULL,NULL}; structure_grid_ini(&gr,rnp);
  for (i1=0; i1<rnp; i1++){ gr.k->data[i1] = rmax * ( (double) i1 / (double) rnp) ; gr.kw->data[i1] = 0.0;}
  if ( op.write_deleta == 1 || op.write_delz == 1 || op.write_F == 1 || op.write_taua == 1 || write_S == 1 ) { write_dyn = 1; }
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
  fprintf(t_file,"#|| 1 t || 2 b(t) || 3 u(t) || 4 uₐ-u(t) ||# \n" );
  if (write_dyn == 1) {
    t_save_file=fopen(t_file_name,"w");
    fprintf(t_save_file,"#|| 1 quench id || 2 t || 3 u(t) || 4 b(t) ||# \n" );
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
  save_dyn_vars_ini(&dyn_save, op, icv.kwr->size, folder, u_char, fname); 
  gsl_vector_memcpy(dyn_save.k,icv.kwr);
  gsl_vector_memcpy(dyn_save.S,icv.Swri);
  if ( write_S == 1 ) {
    s_grid_save_file( Sg, folder, u_char_S, fname ); 
    gsl_vector_radial_distribution_3D(gr.S, gr.k, icv.lpi.rho, Sg); 
    s_grid_save_file( gr, folder, u_char_g, fname );
  }
  dynamics_sph_mono( icv.lpi, dp, Sg, &dyn_save, op, &ds );
  save_dyn_vars_free(&dyn_save);
  /* Computing final state dynamics */
  printf("computing final dynamics\n");
  sprintf(u_char,"%s%03d%s","u_",999,"_");
  sprintf(u_char_S,"%s%03d%s","u_",999,"_S_"); 
  sprintf(u_char_g,"%s%03d%s","u_",999,"_g_");
  save_dyn_vars_ini(&dyn_save, op, icv.kwr->size, folder, u_char, fname);
  gsl_vector_memcpy(dyn_save.k,icv.kwr);
  if ( ua < 1E40){
    sku_inst_change_mono_sph(&Sg.S,icv.Si, icv.Sf, aux_a, ua);
    sku_inst_change_mono_sph(&dyn_save.S, icv.Swri, icv.Swrf, aux_a_w, ua);
    dynamics_sph_mono( icv.lpf, dp, Sg, &dyn_save, op, &dsf );
  }
  else{
    gsl_vector_memcpy(Sg.S,icv.Sf);
    gsl_vector_memcpy(dyn_save.S,icv.Swrf);
    dynamics_sph_mono( icv.lpf, dp, Sg, &dyn_save, op, &dsf );
  }
  if ( write_S == 1 ) {
    s_grid_save_file( Sg, folder, u_char_S, fname ); 
    gsl_vector_radial_distribution_3D(gr.S, gr.k, icv.lpi.rho, Sg); 
    s_grid_save_file( gr, folder, u_char_g, fname );
    }
  save_dyn_vars_free(&dyn_save);
  /* Variables initialization for the waiting times loop */
  convergence = 0;
  du      = 1E-7 * ds.Dl; /* Initial u differential with dt≈du/Dl; where Dl is the long time diffusion coefficient of the initial state */
  t       = 0.0;
  u       = 0.0;
  t_save  = 1E-3;
  du_dt_1 = 0.0;
  i_save  = 0;
  if (write_dyn == 1) {fprintf( t_save_file,"%d \t %1.9e \t %1.9e \t %1.9e \n", i_save, 0.0, 0.0, ds.Dl );}
  fprintf( t_file,"%1.9e \t %1.9e \t %1.9e \t %1.9e \n", 0.0, ds.Dl, 0.0, ua );
  while ( convergence == 0 ){
    /* Computing the dynamics for the next u-value */
    Dl_pre = ds.Dl;
    u+=du;
    sku_inst_change_mono_sph(&Sg.S, icv.Si, icv.Sf, aux_a, u);
    dynamics_sph_mono( icv.lpf, dp, Sg, &dyn_save, no_save, &ds );
    t_pre=t;
    dt = du / ds.Dl;
    t += dt;
    /* Computing the change in du in terms of change in (du/dt)  */
    if( Dl_pre / ds.Dl > 1.05) {du *= 0.05 / ( ( Dl_pre / ds.Dl ) - 1.0) ;}
    else if( ds.Dl / Dl_pre > 1.05) {du *= 0.05 / ( ( ds.Dl / Dl_pre ) - 1.0);}
    else if ( u + 2.0*du < ua && 2.0 * du < 0.05 * ua ) { du *= 2.0; }
    /* Computing the decision to end the time loop */
    e_Dl = fabs(1.0 - ( ds.Dl / dsf.Dl )); /* Error between current Dl and final Dl */
    if ( e_Dl < tol_Dl || u >= ua || ds.Dl <=tol_Dl ){convergence=1;}
    if ( ds.Dl > tol_Dl ) {
      printf("%s %1.9e \t %s %1.9e \t %s %1.9e\n","Dl rel error w final value: ",e_Dl,"u: ", u,"t: ",t);
      fprintf( t_file,"%1.9e \t %1.9e \t %1.9e \t %1.9e \n", t, ds.Dl, u, ua-u );
      fflush(t_file);
      /* Computing the dynamics for the time-save value */
      
      while ( t_pre <= t_save && t > t_save && write_dyn == 1 ){
        i_save +=1;
        sprintf(u_char,"%s%03d%s","u_",i_save,"_");
        sprintf(u_char_S,"%s%03d%s","u_",i_save,"_S_");
        sprintf(u_char_g,"%s%03d%s","u_",i_save,"_g_"); 
        save_dyn_vars_ini( &dyn_save, op, icv.kwr->size, folder, u_char, fname);
        gsl_vector_memcpy(dyn_save.k,icv.kwr);
        du_save = ( t_save - t_pre ) * 0.5 * ( ds.Dl + Dl_pre );
        u_save  = u-du + du_save;
        sku_inst_change_mono_sph(&Sg.S, icv.Si, icv.Sf, aux_a, u_save);
        sku_inst_change_mono_sph(&dyn_save.S, icv.Swri, icv.Swrf, aux_a_w, u_save);
        s_grid_save_file( Sg, folder, u_char_S, fname );
        gsl_vector_radial_distribution_3D(gr.S, gr.k, icv.lpi.rho, Sg); 
        s_grid_save_file( gr, folder, u_char_g, fname );
        dynamics_sph_mono( icv.lpf, dp, Sg, &dyn_save, op, &ds_save );
        fprintf(t_save_file,"%d \t %1.9e \t %1.9e \t %1.9e \n", i_save, t_save, u_save, ds_save.Dl );
        fflush(t_save_file);
        t_save *= t_save_scale;
        save_dyn_vars_free(&dyn_save);
      }
    }
  }
  fprintf(t_file,"%1.9e \t %1.9e \t %1.9e \t %1.9e \n", t, dsf.Dl, ua, 0.0 );
  if (write_dyn == 1) {fclose(t_save_file);}
  fclose(t_file);
  gsl_vector_free(aux_a);
  gsl_vector_free(aux_a_w);
  free(u_char);
  free(u_char_S);
  free(u_char_g);
  structure_grid_free(&Sg);
  structure_grid_free(&gr);
}
