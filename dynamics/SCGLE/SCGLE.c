#include "SCGLE.h"


/* Initialization Methods for data structures */

dynamics_parameters dynamics_parameters_manual_ini(int st, int it, double dtau, double kc, double D0, int mindn,int maxdn ){
	dynamics_parameters dp={st,it,dtau, kc, D0, mindn,maxdn};
	return dp;
}

dynamics_parameters dynamics_parameters_auto_ini(){
	dynamics_parameters dp={10,128,1e-7,2.0*M_PI*1.305,1.0,1E-7,20,60};
	return dp;
}

dynamics_parameters dynamics_parameters_auto_ini_HD(){
	dynamics_parameters dp={10,512,1e-7,4.0*M_PI,1.0,1E-7,20,60};
	return dp;
}

dynamics_save_options dynamics_save_options_auto_ini(){
	dynamics_save_options dso;
	dso.gamma=1;
	dso.Dl=1;
	dso.msd = 1;
	dso.Dt = 1;
	dso.delta_zeta=1;
	dso.tau_alpha=1;
	dso.delta_eta=1;
	dso.Fc=1;
	dso.Fs=1;
	return dso;
}

dynamics_save_options dynamics_save_options_no_save_ini(){
	dynamics_save_options dso;
	dso.delta_eta=0;
	dso.delta_zeta=0;
	dso.Fc=0;
	dso.Fs=0;
	dso.tau_alpha=0;
	dso.msd = 0;
	dso.Dt = 0;
	dso.gamma=1;
	dso.Dl=1;
	return dso;
}

/* Auxiliary function to know the amount of save variables to write  */
int dynamics_save_options_sum_tau(dynamics_save_options dso) {
	int sum;
	sum = dso.delta_eta + dso.delta_zeta + dso.Fc + dso.Fs + dso.msd;
	return sum;
}
int dynamics_save_options_sum_tau_only(dynamics_save_options dso) {
	int sum;
	sum = dso.delta_eta + dso.delta_zeta + dso.msd;
	return sum;
}
/* Auxiliary function to know the amount of save variables to write  */
int dynamics_save_options_sum_k(dynamics_save_options dso) {
	int sum;
	sum = dso.Fc + dso.Fs + dso.tau_alpha;
	return sum;
}

void dynamics_save_variables_spherical_mono_ini( dynamics_save_variables * dsv, const dynamics_save_options dso, gsl_vector * k, gsl_vector * Sk, const int knp, const dynamics_parameters dp, const char * folder, const char * prefix, const char * suffix ){
	const int tfb = 2 * (sizeof(folder) + sizeof(prefix) + sizeof(suffix) + 100*sizeof(char));
	char * dyn_name  = (char *)malloc(tfb);
	char * taua_name = (char *)malloc(tfb);
	char * fc_name   = (char *)malloc(tfb);
	char * fs_name   = (char *)malloc(tfb);
	strcpy(dyn_name, folder); 
	strcat(dyn_name, prefix);
	strcat(dyn_name, "dyn_"); 
	strcat(dyn_name, suffix);
	strcpy(taua_name, folder); 
	strcat(taua_name, prefix);
	strcat(taua_name, "taua_"); 
	strcat(taua_name, suffix);
	strcpy(fc_name, folder); 
	strcat(fc_name, prefix);
	strcat(fc_name, "Fc_");   
	strcat(fc_name, suffix);
	strcpy(fs_name, folder); 
	strcat(fs_name, prefix);
	strcat(fs_name, "Fs_");   
	strcat(fs_name, suffix);
	printf("%s \n", dyn_name);
	dsv->k=NULL;
	dsv->S=NULL;
	dsv->lambda_k=NULL;
	dsv->tau=NULL;
	dsv->delta_zeta_t=NULL;
	dsv->delta_eta_t=NULL;
	dsv->msd=NULL;
	dsv->Dt=NULL;
	dsv->tau_alpha=NULL;
	dsv->Fc=NULL;
	dsv->Fs=NULL;
	if ( dynamics_save_options_sum_k(dso) > 0 ) {
		dsv->k=gsl_vector_alloc(k->size);
		dsv->S=gsl_vector_alloc(k->size);
		dsv->lambda_k = gsl_vector_alloc(k->size);
		gsl_vector_memcpy(dsv->k, k);
		gsl_vector_memcpy(dsv->S, Sk);
	}
	if ( dynamics_save_options_sum_tau(dso)>0 ) { 
		dsv->tau = gsl_vector_alloc(dp.it);
		if ( dso.delta_zeta == 1 ) { dsv->delta_zeta_t = gsl_vector_alloc(dp.it);}
		if ( dso.delta_eta == 1 ) { dsv->delta_eta_t = gsl_vector_alloc(dp.it);}
		if ( dso.msd == 1 ) { dsv->msd = gsl_vector_alloc(dp.it);}
		if ( dso.Dt == 1 ) { dsv->Dt = gsl_vector_alloc(dp.it); }
		}
	if ( dso.Fc == 1 || dso.Fs == 1 ) { 
		dsv->F_Fc   = fopen( fc_name, "w" ); 
		dsv->Fc = gsl_matrix_alloc( dp.it, k->size );
		dsv->F_Fs   = fopen( fs_name, "w" ); 
		dsv->Fs = gsl_matrix_alloc( dp.it, k->size );
		}
	if ( dso.tau_alpha == 1 ){ 
		dsv->F_taua = fopen(taua_name, "w"); 
		dsv->tau_alpha = gsl_vector_alloc(knp);
		gsl_vector_set_all(dsv->tau_alpha,dp.dtau);
	} 
	if ( dynamics_save_options_sum_tau_only( dso ) > 0 ) { dsv->F_dyn  = fopen( dyn_name, "w" );}
	free( dyn_name );
	free( taua_name );
	free( fc_name );
	free( fs_name );
	return;
}

/*Function that computes the wave vector dependent part of the memory function
	lambda(k,kc)
	"k" -> wave vector magnitude
	"kc" ->  Adjustable parameter, found that for HS systems kc ≈ 2 π * 1.305*/
double lambda_spherical_mono( const double k, const double kc ) {
		double l;
		l = k / kc;
		l = l * l;
		l = 1.0 + l;
		l = 1.0 / l;
		return l;
	}

void gsl_vector_lambda_spherical_mono (gsl_vector ** lamk, const gsl_vector * k, const double kc ){
	;
	if( k!=NULL ) {
		int knp=k->size;
		int i1;
		for (i1=0; i1<knp; ++i1){gsl_vector_set(*lamk,i1, lambda_spherical_mono( k->data[i1], kc ));}
	}
	return;
}

intermediate_times_variables intermediate_times_variables_spherical_mono_ini( const dynamics_parameters dp, const structure_grid Sg){
	int knp = Sg.k->size;
	int nt=dp.it;
	gsl_vector * tau = gsl_vector_alloc(nt);
	gsl_vector * idelz = gsl_vector_alloc(nt);
	gsl_vector * lamk = gsl_vector_alloc(knp); 
	gsl_vector_lambda_spherical_mono( &lamk, Sg.k, dp.kc );
	gsl_matrix * iFc = gsl_matrix_alloc(nt, knp);
	gsl_matrix * iFs = gsl_matrix_alloc(nt, knp);
	int i1;
	double lam_val;
	for(i1=0; i1<nt ; ++i1){ gsl_vector_set(tau,i1,( i1 + 1 ) * dp.dtau); }
	intermediate_times_variables itv = { tau, idelz, lamk, iFc, iFs };
	return itv;
}




void free_dynamics_save_variables( dynamics_save_variables * dsv )
{
	if( dsv->tau !=NULL ) { gsl_vector_free( dsv->tau ); }
	if( dsv->delta_zeta_t !=NULL ) { gsl_vector_free( dsv->delta_zeta_t ); }
	if( dsv->delta_eta_t !=NULL ) { gsl_vector_free( dsv->delta_eta_t ); }
	if( dsv->Fc !=NULL ) { gsl_matrix_free( dsv->Fc ); }
	if( dsv->Fs !=NULL ) { gsl_matrix_free( dsv->Fs ); }
	if( dsv->k !=NULL ) { gsl_vector_free( dsv->k ); }
	if( dsv->S !=NULL ) { gsl_vector_free( dsv->S ); }
	if( dsv->msd != NULL ) { gsl_vector_free( dsv->msd ) ;}
	if( dsv->Dt != NULL ) { gsl_vector_free( dsv->Dt ); }
}

void close_dynamics_save_variables( dynamics_save_options dso, dynamics_save_variables * dsv ){
	if ( dynamics_save_options_sum_tau_only(dso) > 0 ) {fclose(dsv->F_dyn);}
	if ( dso.Fc == 1 ) { fclose(dsv->F_Fc); } 
	if ( dso.Fs == 1) { fclose(dsv->F_Fs);}
	if ( dso.tau_alpha == 1 ) {fclose(dsv->F_taua);}
	return;
}

void free_close_dynamics_save_variables( dynamics_save_options dso, dynamics_save_variables *dsv ) {
	free_dynamics_save_variables( dsv );
	close_dynamics_save_variables( dso, dsv );
	return;
}

void free_intermediate_times_variables(intermediate_times_variables * itv){
	if(itv->tau != NULL) 		{ gsl_vector_free (itv->tau); }
	if(itv->delta_z != NULL) 		{ gsl_vector_free (itv->delta_z); }
	if(itv->lambda_k !=NULL) 		{ gsl_vector_free (itv->lambda_k); }
	if(itv->Fc !=NULL) 			{ gsl_matrix_free (itv->Fc); }
	if(itv->Fs !=NULL) 			{ gsl_matrix_free (itv->Fs); }
}

/*Function that computes the integral for the memory function for a spherical monocomponent system
The function computes

Δζ*(t) = C(d) ∫₀∞ [1-S⁻¹(k)]² Fc(k,t) Fs(k,t) kᵈ⁺¹  dk

/* Inputs type							Variable name		Notes

	liquid_params										lp				Composed of double parameters

	structure_grid								 	Sg				Composed of 3 gsl_vector

	gsl_vector *										Fc				Vector of Fc(k,t)

	gsl_vector *										Fs				Vector of Fs(k,t)


Notes:
	-Currently only works for lp.dim=2 and lp.dim=3

*/
double
delta_zeta_spherical_mono( const liquid_params lp, const structure_grid Sg, const gsl_vector * Fc, const gsl_vector * Fs ){
	int knp = (Sg.k)->size;
 	gsl_vector * integrand = gsl_vector_alloc(knp);
 	double dum1;
 	int i1;
 	double z;
 	for (i1=0; i1<knp; i1++){
		dum1 = 1.0 - ( 1.0 / Sg.S->data[i1]);
		dum1 = dum1 * dum1;
	 	integrand->data[i1] = pow(gsl_vector_get (Sg.k, i1), lp.dim+1.0) * dum1 *
													Fc->data[i1] * Fs->data[i1] * Sg.kw->data[i1];
 	}
	z = 2.0 * lp.dim * (lp.rho) * pow(M_PI,lp.dim-1.0);
 	z = gsl_vector_sum(integrand) / z;
	/* HD adjustment for the memory */
	if (lp.dim==2.0) {
		z = z* (2.0 - 2.075 * lp.phi);
	}
	/* Deallocation of vectors memory */
	gsl_vector_free(integrand);
	return z;
}

double delta_zeta_spherical_mono_test_for_k_points( const liquid_params lp, const structure_grid Sg ){
	int knp = (Sg.k)->size;
 	gsl_vector * integrand = gsl_vector_alloc(knp);
 	double dum1;
 	int i1;
 	double z;
 	for (i1=0; i1<knp; i1++){
		dum1 = 1.0 - ( 1.0 / Sg.S->data[i1]);
		dum1 = dum1 * dum1;
	 	integrand->data[i1] = pow(gsl_vector_get (Sg.k, i1), lp.dim+1.0) * dum1 *
													Sg.S->data[i1] * Sg.kw->data[i1];
 	}
	z = 2.0 * lp.dim * (lp.rho) * pow(M_PI,lp.dim-1.0);
 	z = gsl_vector_sum(integrand) / z;
	/* HD adjustment for the memory */
	if (lp.dim==2.0) {
		z = z* (2.0 - 2.075 * lp.phi);
	}
	/* Deallocation of vectors memory */
	gsl_vector_free(integrand);
	return z;
}

/*
Function that computes the viscosity for a spherical monocomponent system

Δη(t) = (k_BT/60π²)∫dk(0, ∞) (k/S(k;t))⁴ (dS(k;t)/dk)² F²(k,τ;t)

the discretization becomes

Δη(t) = (k_BT/60π²)Σ(0,n) kw (k/S(k;t))⁴ (dS(k;t)/dk)² F²(k,τ;t)

for a simple approx for the derivative dS/dk ≈ ΔS/Δk

Δη(t) = (k_BT/60π²)Σ(0,n) (k/S(k;t))⁴ (ΔS(k;t)*F(k,τ;t))²/Δk

*/

double delta_eta_spherical_mono( const liquid_params lp, const structure_grid Sg, const gsl_vector * Fc ){
	int knp = (Sg.k)->size;
 	gsl_vector * integrand = gsl_vector_alloc(knp);
 	double dk, dS, k, Sk, Fk;
 	int i1;
 	double eta;
 	// simple approximation for the derivative
 	/*for (i1=0; i1<knp; i1++){
 		// dk
		if(i1>0) dk = Sg.k->data[i1]- Sg.k->data[i1-1];
		else dk = Sg.k->data[0];
		// dS
		dS = Sg.S->data[i1] - Sg.S->data[i1-1];
		// k, S(k) and F(k,τ)
		k = Sg.k->data[i1];
		Sk = Sg.S->data[i1];
		Fk = Fc->data[i1];
		if (Fk<1.0E-150) {Fk=0.0;}
		integrand->data[i1] = pow(k/Sk, 4) * pow(dS*Fk, 2) / dk;
 	}*/
 	// approxiating the local derivative as a function of three points
 	double dS_dk, f0, f1, f2, x0, x1, x2;
 	integrand->data[0] = 0.0;//pow((Sg.k->data[0])/(Sg.S->data[0]), 4) 
						//* pow((Sg.k->data[1]- Sg.k->data[0]) * Fc->data[i1], 2) / Sg.k->data[0] ;
 	for (i1=1; i1<knp-1; i1++){
 		// dk
		dk = Sg.k->data[0];
		// derivative dS/dk for three points (x_0, f(x, 0)), (x_1, f(x, 1)) and (x_2, f(x, 2))
		// f'(xj) = f(x0) * (2xj - x1 - x2)/((x0 - x1)*(x0 - x2))
		//			 +	f(x1) * (2xj - x0 - x2)/((x1 - x0)*(x1 - x2))
		//			 +	f(x2) * (2xj - x0 - x1)/((x2 - x0)*(x2 - x1))
		// for xj = x1
		// f'(x1) = f(x0) * (x1 - x2)/((x0 - x1)*(x0 - x2))
		//			 +	f(x1) * (2x1 - x0 - x2)/((x1 - x0)*(x1 - x2))
		//			 +	f(x2) * (x1 - x0)/((x2 - x0)*(x2 - x1))
		f0 = Sg.S->data[i1-1]; f1 = Sg.S->data[i1]; f2 = Sg.S->data[i1+1];
		x0 = Sg.k->data[i1-1]; x1 = Sg.k->data[i1]; x2 = Sg.k->data[i1+1];
		dS_dk = f0 * (x1 - x2)/((x0 - x1)*(x0 - x2))
				+ f1 * (2*x1 - x0 - x2)/((x1 - x0)*(x1 - x2))
				+ f2 * (x1 - x0)/((x2 - x0)*(x2 - x1));
		// k, S(k) and F(k,τ)
		k = Sg.k->data[i1];
		Sk = Sg.S->data[i1];
		Fk = Fc->data[i1];
		if (Fk<1.0E-150) {Fk=0.0;}
		// Δη(t) = (k_BT/60π²)Σ(0,n) dk (k/S(k;t))⁴ (dS(k;t)/dk)² F²(k,τ;t)
		integrand->data[i1] = pow(k/Sk, 4) * pow(dS_dk*Fk, 2) * Sg.kw->data[i1];
 	}
 	integrand->data[knp-1] = 0.0;

	eta = 60.0 * M_PI * M_PI;
 	eta = gsl_vector_sum(integrand) / eta;
	/* Deallocation of vectors memory */
	gsl_vector_free(integrand);
	return eta;
}

/* Auxiliary function that computes the integral ∫₀ᵗ Δζ(t')dt' */
double
int_Delz( const gsl_vector * Delz, const double prev_val, const double dtau, const int ini, const int end ){
	double integral;
	int i1;
	integral = 0;
	for (i1=ini; i1<end; ++i1){
		integral = integral + Delz->data[i1];
	}
	integral = ( dtau * integral ) + prev_val;
	return integral;
}

double Dl_sph_mono( double int_delz ){
	double Dl;
	Dl = 1.0 + int_delz;
	Dl = 1.0 / Dl;
	return Dl;
}

/* Function that computes for the correlation time dependent diffusion coefficient
for the it- time through the convolution D(t)= 1 - ∫₀ᵗ	 Δζ(t-t')D(t')dt' */

double Dt(int it, dynamics_parameters dp, gsl_vector * delta_z, gsl_vector * Dt ){
	int i1;
	double DDz, sum, dum1;
	sum = 0.0;
	for (i1=0; i1 < it; ++i1){
		DDz = delta_z->data[it-i1] - delta_z->data[it-i1-1];
		sum = sum + ( Dt->data[i1] * DDz );
	}
	dum1 = ( Dt->data[it-1] / dp.dtau ) - sum;
	dum1 = dum1 / ( delta_z->data[0] + ( 1.0 / dp.dtau ) );
	return dum1;
}
/* Function that computes for the mean-squared displacement
for the it- time through the convolution w(t)= t - ∫₀ᵗ	[ d ( w(t-t') ) / dt ] [ Δζ(t')dt' ]  */
double msdt(int it, dynamics_parameters dp, gsl_vector * delz, gsl_vector * msd){
	int i1;
	double h=dp.dtau;
	double Dmsd, sum, dum1;
	sum = -( delz->data[0] + ( 1.0 / h ) ) * msd->data[it-1];
	/*sum = msd->data[it-2] * ( h + delz->data[0] ) ;*/
	for (i1=1; i1 < it; ++i1){
		Dmsd = msd->data[it-i1] - msd->data[it-i1-1];
		sum = sum + ( delz->data[i1] * Dmsd );
	}
	dum1 = h * ( 1.0 - sum - ( delz->data[it] * msd->data[0] ) );
	dum1 = dum1 / (1.0 + ( h * delz->data[0] ) ) ;
	/*printf("%1.9e \t %1.9e \t %d \n", dum1, msd->data[it-1], it );*/
	return dum1;
}

/*Function that computes the small time limit of the intermediate scattering function
	Fc(k,t)
	"k" -> wave vector magnitude
	"t" -> correlation time value
	"s" -> static structure factor value */
double ini_Fc_sph_mono( double k, double t, double S, double D0 ){
	double f;
	f =  exp (- D0 * t * k * k / S) * S;
	return f;
}

/*Function that computes the small time limit of the self part of the intermediate scattering function
	Fs(k,t)
	"k" -> wave vector magnitude
	"t" -> correlation time value */
double ini_Fs_sph_mono( double k, double t, double D0 ){
	double f;
	f = exp (- D0 * t * k * k ) ;
	return f;
}

double tau_alpha_sph_mono( double Fs_pre, double tau_pre, double Fs, double tau, double D0, double taua_pre ){
	double taua;
	double m,c;
	const double F_taua = exp( - D0 );
	if ( Fs < F_taua && Fs_pre >= F_taua  ){
		m = ( Fs - Fs_pre ) / ( tau - tau_pre );
		c = Fs - (m*tau );
		taua = ( F_taua - c ) / m;
	}
	else if( Fs > F_taua ){
			taua = tau;
	}
	else{
		taua = taua_pre;
	}
	return taua;
}

/* Function that computes for the small time limits of the SCGLE formalism */
/* Inputs type							Variable name		Notes

	liquid_params										lp				Composed of double parameters

	structure_grid								 	Sg				Composed of 3 gsl_vector

	dyn_params											dp				Composed of integers and doubles that helps in computing the dynamics

																						Composed of dynamics variables computed for an
	itv_vars	*											itv				equally spaced time grid, you need to pass an address
																						as this function pretends to modify such values
*/
void short_times_dynamics_spherical_mono(liquid_params lp, structure_grid Sg, dynamics_parameters dp, dynamics_save_variables * dsv,
intermediate_times_variables * itv, dynamics_save_options dso ){
	int knp = (Sg.k)->size;
	int knp_w; if ( dso.Fc == 1 || dso.Fs == 1 ) { knp_w = dsv->k->size; }
	int i1,i2;
	double t,k_val,S_val,z,e;
	gsl_vector_view Fc_v;
	gsl_vector_view Fs_v;
	double sumDelz,dum1;
	sumDelz=0.0;
	for( i2 = 0; i2 < dp.st; ++i2 ){
		t = gsl_vector_get(itv->tau,i2);
		for( i1 = 0; i1 < knp; ++i1 ){
			k_val = gsl_vector_get (Sg.k, i1);
			S_val = gsl_vector_get (Sg.S, i1);
			gsl_matrix_set ( itv->Fc, i2, i1, ini_Fc_sph_mono( k_val, t, S_val, dp.D0 ) );
			gsl_matrix_set ( itv->Fs, i2, i1, ini_Fs_sph_mono( k_val, t, dp.D0 ) );
		}
		Fc_v = gsl_matrix_row( itv->Fc, i2);
		Fs_v = gsl_matrix_row( itv->Fs, i2);
		z = delta_zeta_spherical_mono( lp, Sg,&Fc_v.vector,&Fs_v.vector );
		gsl_vector_set ( itv->delta_z, i2, z );
		/* Dynamics save variables computation */
		/* Fc and Fs */
		if ( dso.Fc == 1 || dso.Fs == 1  ) {
			for( i1 = 0; i1 < knp_w; ++i1 ){
				k_val = gsl_vector_get (dsv->k, i1);
				S_val = gsl_vector_get (dsv->S, i1);
				gsl_matrix_set ( dsv->Fc, i2, i1, ini_Fc_sph_mono( k_val, t, S_val, dp.D0 ) ); 
				gsl_matrix_set ( dsv->Fs, i2, i1, ini_Fs_sph_mono( k_val, t, dp.D0 ) );
				}
		}
		/* Computation for tau_alpha */
		if ( dso.tau_alpha == 1 ) {
			for( i1 = 0; i1 < knp; ++i1 ){
				if (i2>0) {
					dsv->tau_alpha->data[i1] = tau_alpha_sph_mono( gsl_matrix_get(itv->Fs,i2-1,i1), itv->tau->data[i2-1], gsl_matrix_get(itv->Fs,i2,i1), itv->tau->data[i2], dp.D0, dsv->tau_alpha->data[i1] ) ;
				}
			}
		}
		/* Computation of Deta(t) */
		if ( dso.delta_eta == 1 ){
		e = delta_eta_spherical_mono( lp, Sg,&Fc_v.vector);
		gsl_vector_set ( dsv->delta_eta_t, i2, e );
		}
		/* Computation of D(t) */
		sumDelz = sumDelz + z;
		dum1 = dp.D0 * ( 1.0 - dp.dtau * sumDelz );
		if ( dso.Dt == 1 ) {gsl_vector_set ( dsv->Dt, i2, dum1 );};
		/* Computation of w(t) */
		if ( dso.msd == 1 ) { 
			if (i2==0) { dsv->msd->data[0] = dum1 * dp.dtau; }
			else { 
				for( i1 = 1; i1< dp.st; ++i1 ){
				dsv->msd->data[i1] = dsv->msd->data[i1-1] + dsv->Dt->data[i1] * dp.dtau;
				} 
			} 
		}
	}	
}



/* Function that computes for auxiliary vectors as and ac employed in medium_t_dynamics_sph_mono */
void alphas( structure_grid Sg, dynamics_parameters dp, gsl_vector * lamk,
double mem1, gsl_vector ** ac, gsl_vector ** as){
	int knp = (Sg.k)->size;
	int i1;
	double k_val,S_val,lam_val;
	double ac_val, as_val, adum;
	for( i1 = 0; i1 < knp; ++i1 ){
		k_val   = gsl_vector_get (Sg.k, i1);
		lam_val = gsl_vector_get (lamk, i1);
		S_val   = gsl_vector_get (Sg.S, i1);
		adum    = ( dp.D0 * k_val * k_val );
		as_val  = ( 1.0 / dp.dtau ) + ( lam_val * mem1 );
		ac_val  = as_val + ( adum / S_val );
		as_val  = as_val + adum;
		as_val  = 1.0 / as_val;
		ac_val  = 1.0 / ac_val;
		gsl_vector_set (*as, i1, as_val);
		gsl_vector_set (*ac, i1, ac_val);
	}
}

/* Function that computes for the dynamics variables employing the SCGLE formalism */
/* Inputs type							Variable name		Notes
	liquid_params										lp				Composed of double parameters
	structure_grid								 	Sg				Composed of 3 gsl_vector
	dyn_params											dp				Composed of integers and doubles that helps in computing the dynamics
	itv_vars												itv				Composed of dynamics variables computed for an equally spaced time grid
The SCGLE formalism consists in finding a solution to the next coupled equations:
	d[Fc(t)] /dt = λΔζ(t)S - k²D₀S⁻¹Fc(t) - λ d[ ∫₀ᵗ Δζ(t-t') Fc(t') dt' ]/dt
	d[Fs(t)] /dt = λΔζ(t) -   k²D₀Fs(t)   - λ d[ ∫₀ᵗ Δζ(t-t') Fs(t') dt' ]/dt
	Δζ(t) =  c ∫ [1-S⁻¹]² kᵈ⁺¹ Fc Fs dk,
where
	c = V(d) / (2π)ᵈ ρ,
with V(d) being the volume of a unitary d-dimensional sphere
Variables:
	k    = Sg.k
	dk  = Sg.kw
	S(k) = Sg.S
	λ(k) = lamk
	d    = lp.d
	ρ    = lp.rho
	γ    = Function Output
Notes:
	- Currently only works for d=2 and d=3
	*/
void intermediate_times_dynamics_spherical_mono_for_writting(const liquid_params lp, dynamics_parameters dp, intermediate_times_variables * itv, dynamics_save_variables * dsv){
		int knp = dsv->k->size;
		int i1,i2,i3,i4;
		double t,k_val,S_val,z;
		double dum1,deldelz_val,delz_val,dele_val;
		double delz_test;
		structure_grid Sg_dummy = { dsv->k, NULL, dsv->S }; 
		gsl_vector * ac   		= gsl_vector_alloc(knp);
		gsl_vector * as   		= gsl_vector_alloc(knp);
		gsl_vector * lamc 		= gsl_vector_alloc(knp);
		gsl_vector * lams 		= gsl_vector_alloc(knp);
		gsl_vector * dFc1 		= gsl_vector_alloc(knp);
		gsl_vector * dFs1 		= gsl_vector_alloc(knp);
		gsl_vector * Fcdum 		= gsl_vector_alloc(knp);
		gsl_vector * Fsdum 		= gsl_vector_alloc(knp);
		gsl_vector * Fcdum2 	= gsl_vector_alloc(knp);
		gsl_vector * Fsdum2 	= gsl_vector_alloc(knp);
		gsl_vector * deldelz	= gsl_vector_alloc(dp.it);
		gsl_vector_view Fc_v,Fc1;
		gsl_vector_view Fs_v,Fs1;
		double dum_val1,dum_val2;
		double Fc_val,Fs_val;
		/* Initialization of variables */
		alphas( Sg_dummy, dp, itv->lambda_k, itv->delta_z->data[0], &ac, &as);
		gsl_vector_memcpy(lamc,ac); 
		gsl_vector_mul(lamc, dsv->lambda_k); /* lamc = ac * lamk */
		Fc1 = gsl_matrix_row( dsv->Fc, 0); 
		gsl_vector_memcpy(lams,as); 
		gsl_vector_mul(lams, dsv->lambda_k); /* lams = as * lamk */
		Fs1 = gsl_matrix_row( dsv->Fs, 0); /* Fc1=Fc(0); Fs1=Fs(0) */
		gsl_vector_set_all(dFs1,1.0); 
		gsl_vector_sub(dFs1,&Fs1.vector); /* dFs1 = 1-Fs(0) */
		gsl_vector_memcpy(dFc1,dsv->S); 
		gsl_vector_sub(dFc1,&Fc1.vector); /* dFc1 = S -Fc(0) */
			/* deldelz(i) = delz(i) - delz(i-1) */
		gsl_vector_set_zero(deldelz);
		for(i1=1; i1<dp.it; ++i1){
			dum1 = gsl_vector_get(itv->delta_z,i1) - gsl_vector_get(itv->delta_z,i1-1);
			gsl_vector_set(deldelz,i1,dum1);
		}
		
		/* Computing the dynamic variables for new times */
		for(i1=dp.st; i1<dp.it; ++i1){
			gsl_vector_set_zero(Fcdum); /* Fcdum=0 */
			gsl_vector_set_zero(Fsdum); /* Fsdum=0 */
			for(i2=1; i2<i1; ++i2){
				Fc_v=gsl_matrix_row( dsv->Fc, i2 );
				Fs_v=gsl_matrix_row( dsv->Fs, i2 );
				deldelz_val = deldelz->data[i1-i2];
				for(i3=0; i3<knp; ++i3){
					/*  */
					Fc_val = gsl_vector_get(&Fc_v.vector,i3); /* Fc(i2) */
					Fs_val = gsl_vector_get(&Fs_v.vector,i3); /* Fs(i2) */
					Fcdum->data[i3] = Fcdum->data[i3] - ( Fc_val * deldelz_val ) ;
					Fsdum->data[i3] = Fsdum->data[i3] - ( Fs_val * deldelz_val ) ;
				}
			}
			Fc_v = gsl_matrix_row( dsv->Fc, i1-1 ); 
			Fs_v = gsl_matrix_row( dsv->Fs, i1-1 );
			delz_val = itv->delta_z->data[i1-1];
			for(i3=0; i3<knp; ++i3){
				Fc_val = gsl_vector_get(&Fc1.vector,i3); /* F(0) */
				dum_val1 = ( Fcdum->data[i3] + ( Fc_val * delz_val ) ) * dsv->lambda_k->data[i3];
				dum_val1 = dum_val1  + gsl_vector_get(&Fc_v.vector,i3) / dp.dtau;
				dum_val1 = dum_val1 * ac->data[i3];
				gsl_vector_set( Fcdum, i3, dum_val1 );
				Fcdum2->data[i3] = Fcdum->data[i3] + ( dFc1->data[i3] * itv->delta_z->data[i1] * lamc->data[i3] );
				gsl_matrix_set(dsv->Fc,i1,i3, gsl_vector_get(Fcdum2,i3) );

				Fs_val = gsl_vector_get(&Fs1.vector,i3); /* F(0) */
				dum_val1 = ( Fsdum->data[i3] + ( Fs_val * delz_val)  ) * dsv->lambda_k->data[i3];
				dum_val1 = dum_val1  + gsl_vector_get(&Fs_v.vector,i3) / dp.dtau;
				dum_val1 = dum_val1 * as->data[i3];
				gsl_vector_set( Fsdum, i3, dum_val1 );
				Fsdum2->data[i3] = Fsdum->data[i3] + ( dFs1->data[i3] * itv->delta_z->data[i1] * lams->data[i3] );
				gsl_matrix_set(dsv->Fs,i1,i3, gsl_vector_get(Fsdum2,i3) );
			}
		}
	/* Freeing memory */
	gsl_vector_free(ac);
	gsl_vector_free(as);
	gsl_vector_free(lamc);
	gsl_vector_free(lams);
	gsl_vector_free(dFc1);
	gsl_vector_free(dFs1);
	gsl_vector_free(Fcdum);
	gsl_vector_free(Fsdum);
	gsl_vector_free(Fcdum2);
	gsl_vector_free(Fsdum2);
	gsl_vector_free(deldelz);
	return;
}

/* Function that computes for the dynamics variables employing the SCGLE formalism */
/* Inputs type							Variable name		Notes
	liquid_params										lp				Composed of double parameters
	structure_grid								 	Sg				Composed of 3 gsl_vector
	dyn_params											dp				Composed of integers and doubles that helps in computing the dynamics
	itv_vars												itv				Composed of dynamics variables computed for an equally spaced time grid
The SCGLE formalism consists in finding a solution to the next coupled equations:
	d[Fc(t)] /dt = λΔζ(t)S - k²D₀S⁻¹Fc(t) - λ d[ ∫₀ᵗ Δζ(t-t') Fc(t') dt' ]/dt
	d[Fs(t)] /dt = λΔζ(t) -   k²D₀Fs(t)   - λ d[ ∫₀ᵗ Δζ(t-t') Fs(t') dt' ]/dt
	Δζ(t) =  c ∫ [1-S⁻¹]² kᵈ⁺¹ Fc Fs dk,
where
	c = V(d) / (2π)ᵈ ρ,
with V(d) being the volume of a unitary d-dimensional sphere
Variables:
	k    = Sg.k
	dk  = Sg.kw
	S(k) = Sg.S
	λ(k) = lamk
	d    = lp.d
	ρ    = lp.rho
	γ    = Function Output
Notes:
	- Currently only works for d=2 and d=3
	*/
void intermediate_times_dynamics_spherical_mono( const liquid_params lp, const structure_grid Sg, dynamics_parameters dp, intermediate_times_variables* itv, dynamics_save_variables * dsv, dynamics_save_options dso ){
	const double tol=1e-10;
	int knp = (Sg.k)->size;
	int i1,i2,i3,i4;
	double t,k_val,S_val,z,error;
	double dum1,deldelz_val,delz_val,m,c,dele_val;
	double delz_test;
	gsl_vector * ac = gsl_vector_alloc(knp);
	gsl_vector * as = gsl_vector_alloc(knp);
	gsl_vector * lamc = gsl_vector_alloc(knp);
	gsl_vector * lams = gsl_vector_alloc(knp);
	gsl_vector * dFc1 = gsl_vector_alloc(knp);
	gsl_vector * dFs1 = gsl_vector_alloc(knp);
	gsl_vector * Fcdum = gsl_vector_alloc(knp);
	gsl_vector * Fsdum = gsl_vector_alloc(knp);
	gsl_vector * Fcdum2 = gsl_vector_alloc(knp);
	gsl_vector * Fsdum2 = gsl_vector_alloc(knp);
	gsl_vector * deldelz = gsl_vector_alloc(dp.it);
	gsl_vector_view Fc_v,Fc1;
	gsl_vector_view Fs_v,Fs1;
	double dum_val1,dum_val2;
	double Fc_val,Fs_val;
	/* Initialization of variables */
	alphas( Sg, dp, itv->lambda_k, itv->delta_z->data[0], &ac, &as);
	gsl_vector_memcpy(lamc,ac); gsl_vector_mul(lamc, itv->lambda_k); /* lamc = ac * lamk */
	gsl_vector_memcpy(lams,as); gsl_vector_mul(lams, itv->lambda_k); /* lams = as * lamk */
	Fc1 = gsl_matrix_row( itv->Fc, 0); Fs1 = gsl_matrix_row( itv->Fs, 0); /* Fc1=Fc(0); Fs1=Fs(0) */
	gsl_vector_set_all(dFs1,1.0); gsl_vector_sub(dFs1,&Fs1.vector); /* dFs1 = 1-Fs(0) */
	gsl_vector_memcpy(dFc1,Sg.S); gsl_vector_sub(dFc1,&Fc1.vector); /* dFc1 = S -Fc(0) */
		/* deldelz(i) = delz(i) - delz(i-1) */
	gsl_vector_set_zero(deldelz);
	for(i1=1; i1<dp.st; ++i1){
		dum1 = gsl_vector_get(itv->delta_z,i1) - gsl_vector_get(itv->delta_z,i1-1);
		gsl_vector_set(deldelz,i1,dum1);
	}
	/* Computing the dynamic variables for new times */
	for(i1=dp.st; i1<dp.it; ++i1){
		gsl_vector_set_zero(Fcdum); /* Fcdum=0 */
		gsl_vector_set_zero(Fsdum); /* Fsdum=0 */
		for(i2=1; i2<i1; ++i2){
			Fc_v=gsl_matrix_row( itv->Fc, i2 );
			Fs_v=gsl_matrix_row( itv->Fs, i2 );
			deldelz_val = deldelz->data[i1-i2];
			for(i3=0; i3<knp; ++i3){
				/*  */
				Fc_val = gsl_vector_get(&Fc_v.vector,i3); /* Fc(i2) */
				Fs_val = gsl_vector_get(&Fs_v.vector,i3); /* Fs(i2) */
				Fcdum->data[i3] = Fcdum->data[i3] - ( Fc_val * deldelz_val ) ;
				Fsdum->data[i3] = Fsdum->data[i3] - ( Fs_val * deldelz_val );
			}
		}
		Fc_v=gsl_matrix_row( itv->Fc, i1-1 ); Fs_v=gsl_matrix_row( itv->Fs, i1-1 );
		delz_val = itv->delta_z->data[i1-1];
		for(i3=0; i3<knp; ++i3){
			Fc_val = gsl_vector_get(&Fc1.vector,i3); /* F(0) */
			dum_val1 = ( Fcdum->data[i3] + ( Fc_val * delz_val ) ) * itv->lambda_k->data[i3];
			dum_val1 = dum_val1  + gsl_vector_get(&Fc_v.vector,i3) / dp.dtau;
			dum_val1 = dum_val1 * ac->data[i3];
			gsl_vector_set( Fcdum, i3, dum_val1 );

			Fs_val = gsl_vector_get(&Fs1.vector,i3); /* F(0) */
			dum_val1 = ( Fsdum->data[i3] + ( Fs_val * delz_val)  ) * itv->lambda_k->data[i3];
			dum_val1 = dum_val1  + gsl_vector_get(&Fs_v.vector,i3) / dp.dtau;
			dum_val1 = dum_val1 * as->data[i3];
			gsl_vector_set(Fsdum, i3, dum_val1);
		}
		/* Iteration */
		error = 1.0; /* setting a new error value to initiate the iteration */
		m = deldelz->data[i1-1]/dp.dtau;
		c = itv->delta_z->data[i1-1] - m * itv->tau->data[i1-1];
		delz_test = m * itv->tau->data[i1] + c; /* gsl_vector_get(itv->delz,i1-1) - gsl_vector_get(itv->delz,i1-1); /* setting a test value for Δζ(i1) */
		while ( error > tol ) {
			/* Computing a new value for Fc(i1) and Fs(i1) in terms of the test value of Δζ(i1) */
			for (i3=0; i3<knp; ++i3){
				Fcdum2->data[i3] = Fcdum->data[i3] + ( dFc1->data[i3] * delz_test * lamc->data[i3] );
				Fsdum2->data[i3] = Fsdum->data[i3] + ( dFs1->data[i3] * delz_test * lams->data[i3] );
			}
			/* Computing Δζ with the computed values for Fc(i1) and and Fs(i1) */
			delz_val = delta_zeta_spherical_mono( lp, Sg, Fcdum2, Fsdum2 );

			/* Computing the relative error between the test Δζ value and the computed Δζ  */
			error = fabs(1.0 - (delz_test / delz_val ));
			/* Updating a new value to test */
			delz_test = delz_val;
		}
		/*printf( "%s %d %s %1.9e \n", "time=", i1, "delta_zeta=", delz_val);*/
		/* Saving the data for the i1-time */
		gsl_vector_set(itv->delta_z,i1,delz_val);
		for ( i2=0; i2 < knp; ++i2){
			gsl_matrix_set(itv->Fc,i1,i2, gsl_vector_get(Fcdum2,i2) );
			gsl_matrix_set(itv->Fs,i1,i2, gsl_vector_get(Fsdum2,i2) );
		}
		deldelz->data[i1] = gsl_vector_get(itv->delta_z,i1) - gsl_vector_get(itv->delta_z,i1-1);
		if (dso.delta_eta==1){
				dele_val = delta_eta_spherical_mono( lp, Sg, Fcdum2 );
				gsl_vector_set(dsv->delta_eta_t,i1,dele_val);
		}
	}
	/* Saving for writting options once the intermediate times have been computed */
	if ( dynamics_save_options_sum_tau_only(dso) > 0 ) {
		/* update tau */
		gsl_vector_memcpy(dsv->tau,itv->tau);
		/* update Delta Zeta */
		if ( dso.delta_zeta == 1 ) { gsl_vector_memcpy(dsv->delta_zeta_t,itv->delta_z); }
		
		/* update D(t) */
		if ( dso.Dt == 1 ) { 
			for(i1=dp.st; i1<dp.it; ++i1) { dsv->Dt->data[i1]  = Dt(i1, dp, itv->delta_z, dsv->Dt ); } 
			}
		
		/* update w(t) */
		if ( dso.msd == 1 ) { 
			for(i1=dp.st; i1<dp.it; ++i1) { dsv->msd->data[i1] = msdt(i1, dp, itv->delta_z, dsv->msd); }
			}
	}
	/* update tau_alpha */
	if ( dso.tau_alpha==1 ){
		for(i1=dp.st; i1<dp.it; ++i1) {
			for ( i2=0; i2 < knp; ++i2){ 
				dsv->tau_alpha->data[i2] = tau_alpha_sph_mono( 
					gsl_matrix_get(itv->Fs,i1-1,i2), itv->tau->data[i1-1],
					gsl_matrix_get(itv->Fs,i1,i2), itv->tau->data[i1],
					dp.D0,dsv->tau_alpha->data[i2] ) ;
			}
		}	
	}
	/* update Fc and Fs */
	if( dso.Fc == 1 || dso.Fs == 1 ){
		intermediate_times_dynamics_spherical_mono_for_writting(lp, dp, itv, dsv);
	}
	/* Freeing memory */
	gsl_vector_free(ac);
	gsl_vector_free(as);
	gsl_vector_free(lamc);
	gsl_vector_free(lams);
	gsl_vector_free(dFc1);
	gsl_vector_free(dFs1);
	gsl_vector_free(Fcdum);
	gsl_vector_free(Fsdum);
	gsl_vector_free(Fcdum2);
	gsl_vector_free(Fsdum2);
	gsl_vector_free(deldelz);
}



void save_half_intermediate_times_variables_spherical_mono( intermediate_times_variables * itv, dynamics_parameters * dp, dynamics_save_options dso, dynamics_save_variables * dsv ){
	int i1,i2,isave;
	int nt =  dp->it;
	const int nth = nt / 2;
	const int knp = itv->Fc->size2;
	int knp_w; if ( dso.Fc == 1 || dso.Fs == 1 ) { knp_w = dsv->Fc->size2; }
	/* Saving for the intermediate times variables */
	dp->dtau = 2.0 * dp->dtau;
	gsl_vector_scale(itv->tau,2.0);
	for(i1=0; i1<nth; ++i1){
		isave=1+i1*2;
		gsl_vector_set( itv->delta_z, i1,
			0.5 * (gsl_vector_get( itv->delta_z,isave)+ gsl_vector_get(itv->delta_z,isave-1)) );
		for(i2=0; i2<knp; ++i2){
			gsl_matrix_set(itv->Fc,i1,i2,0.5*(gsl_matrix_get(itv->Fc,isave,i2)+gsl_matrix_get(itv->Fc,isave-1,i2)));
			gsl_matrix_set(itv->Fs,i1,i2,0.5*(gsl_matrix_get(itv->Fs,isave,i2)+gsl_matrix_get(itv->Fs,isave-1,i2)));
		}
	}
	/* Saving for the dynamics_save_variables; 
		Note: dsv.tau and dsv.delta_zeta_t are updated in the intermediate times computation  */
	/* w(t) */
	if ( dso.msd == 1 ) {
		for(i1=0; i1<nth; ++i1){
			isave=1+i1*2;
			gsl_vector_set( dsv->msd, i1, 0.5 * (gsl_vector_get(dsv->msd,isave) + 
				gsl_vector_get(dsv->msd,isave-1)));
		}
	}
	/* D(t) */
	if ( dso.Dt == 1 ) {
		for(i1=0; i1<nth; ++i1){
			isave=1+i1*2;
			gsl_vector_set( dsv->Dt, i1, 0.5 * ( gsl_vector_get( dsv->Dt, isave ) +
				gsl_vector_get( dsv->Dt, isave - 1 ) ) );
		}
	}
	/* Delta eta (t) */
	if ( dso.delta_eta == 1 ) {
		for(i1=0; i1<nth; ++i1){
			isave=1+i1*2;
			gsl_vector_set( dsv->delta_eta_t, i1, 0.5 * (gsl_vector_get( 
				dsv->delta_eta_t, isave ) + gsl_vector_get( dsv->delta_eta_t, isave ) ) );
		}
	}
	/* Fc(k,t) y Fs(k,t) */
	if ( dso.Fc == 1 || dso.Fs == 1 ) {
		for( i1=0; i1<nth; ++i1 ){
			isave=1+i1*2;
			for(i2=0; i2<knp_w; ++i2){
				gsl_matrix_set( dsv->Fc, i1, i2, 0.5 * ( gsl_matrix_get( dsv->Fc, isave, i2 )
					+ gsl_matrix_get( dsv->Fc, isave - 1, i2 ) ) );
				gsl_matrix_set( dsv->Fs, i1, i2, 0.5 * ( gsl_matrix_get( dsv->Fs, isave, i2 )
					+ gsl_matrix_get( dsv->Fs, isave - 1, i2 ) ) );
				}
			}
		}
	return ;
}


/* Function that computes for gamma iteratively */
/* Inputs type							Variable name								Notes
	structure grid								 Sg							Composed of 3 gsl_vector
	gsl_vector										lamk
	liquid_params										lp					Composed of double parameters
The function to solve is:
	γ =  c ∫ { [S(k)-1]λ(k) }² kᵈ⁺¹ / {[λ(k)S(k) + γ k²] [ λ(k) + γ k² ]} dk,
where
	c = V(d) / (2π)ᵈ ρ,
with V(d) being the volume of a unitary d-dimensional sphere
Variables:
	k    = Sg.k
	dk  = Sg.kw
	S(k) = Sg.S
	λ(k) = lamk
	d    = lp.d
	ρ    = lp.rho
	γ    = Function Output
Notes:
	- Currently only works for d=2 and d=3
	- γ is only computed for values < 1E40
	- If γ = 1E40 the value is expected to diverge
	- If γ = 1E99 no convergence was found in 10,000 iteration steps
	*/
double gamma_spherical_mono(structure_grid Sg, gsl_vector * lamk, liquid_params lp){
	const double tol=1e-8;
	const double gamma_max=1e40;
	int knp = Sg.k->size;
	int i1,i2;
	int convergence;
	double dimen_dum;
	gsl_vector * k2   = gsl_vector_alloc(knp);
	gsl_vector * lams = gsl_vector_alloc(knp);
	gsl_vector * dum1 = gsl_vector_alloc(knp);
	double gamma, gamma_test,error;
	double dum_val1,integral_val, dimen_2_correction;
	/* dimen_dum = V(d) / (2π)ᵈ ρ; ONLY WORKS FOR d=2 and d=3 */
	dimen_dum = lp.dim * pow(2.0 * M_PI,lp.dim) * lp.rho;
	if (lp.dim==3.0){dimen_dum *= 1.0 / ( 4.0 * M_PI );}
	if (lp.dim==2.0){dimen_dum *= 1.0 / ( (2.0 - 2.075 * lp.phi) * 2.0 * M_PI);}
	/*if(lp.dim==2){dimen_dum *=(2.0 - 2.075 * lp.phi);} /* HD correction */
	/* k2=k² */
	gsl_vector_memcpy(k2,Sg.k);
	gsl_vector_mul(k2,k2);
	/* lams = λ(k) S(k) */
	gsl_vector_memcpy(lams,lamk);
	gsl_vector_mul(lams,Sg.S);
	for (i1=0; i1<knp; ++i1){
		/* dum_val1 = {[S(k)-1]λ(k)}²kᵈ⁺¹dk */
		dum_val1 = (Sg.S->data[i1] - 1.0) * lamk->data[i1];
		dum_val1 = dum_val1 * dum_val1 * Sg.kw->data[i1] * pow( Sg.k->data[i1], lp.dim + 1.0 );
		/* Assign dum_val1 to dum1 gsl vector */
		gsl_vector_set(dum1, i1, dum_val1 );
	}
	/* Initialization of loop variable and proposal of γ */
	convergence = 0;
	gamma = 1.0e-7;
	i2=0;
	while ( convergence == 0 ){
		/* Initialization of integral variable */
		integral_val = 0.0;
		for (i1=0; i1<knp; ++i1){
			/* dum_val1 = [λ(k) + S(k) γ k²] [ λ(k) + γ k² ] */
			dum_val1 = gamma * k2->data[i1] ;
			dum_val1 = ( lams->data[i1] + dum_val1 ) * ( lamk->data[i1] + dum_val1 );
			/* Computing the integral ∫ [S(k)-1]λ(k) kᵈ⁺¹ / {[λ(k)S(k) + γk²] [ λ(k) + γk² ]} dk */
			integral_val = integral_val + ( dum1->data[i1] / dum_val1 );
		}
		/* setting gamma_test to the computed gamma */
		gamma_test = dimen_dum / integral_val;
		/* computing the error between the proposed value of γ and the computed γ */
		error = fabs( 1.0 - ( gamma / gamma_test ) );
		/* setting the proposal to the just computed value */
		gamma = gamma_test;
		/* Loop exit conditions */
		if ( error <= tol ){
			convergence = 1;
		}
		if ( gamma > gamma_max ){
			convergence = 1;
			gamma = gamma_max;
		}
		if ( i2 > 10000 ){
			convergence = 1;
			gamma = 1E99;
		}

	}
	gsl_vector_free(k2);
	gsl_vector_free(lams);
	gsl_vector_free(dum1);
	return gamma;
}

/* 

Function that searches for the liquid parameters condition for dynamical arrest 
when two limits of liquid_parameters are given. The function needs to compute for S(k)
so it uses the structure module to compute the structure factor on the search of the 
liquid parameter values while moving linearly between all the parameters lp0 and lp1.

inputs: 
	char * sys -> string that indicates the system employed over gsl_vector_s_function_selector_mono_sph (defined in structures.c)
	char * approx -> string that indicates the approx employed over gsl_vector_s_function_selector_mono_sph (defined in structures.c)
	liquid_params lp0, lp1 -> Limits to search for arrest
	structure_grid Sg -> Grid in which to evaluate the structure, nodes and weights need to be precomputed
	gsl_vector * lamk -> Array of SCGLE interpolation function λ(k) evaluated on Sg.nodes
	double tol -> tolerance parameter for the lp norm of the difference between arrest and fluid parameters. Accept values R ∈ [1E-14:1E-1]
output:
	liquid_params lp_sup -> liquid parameters for which the system is found to be arrested with a tolerance=

*/
liquid_params arrest_lp_in_limits(char * sys, char * approx, liquid_params lp0, liquid_params lp1, structure_grid Sg, gsl_vector * lamk, double tol){
	liquid_params lp_test,lp_inf,lp_sup,lp_dif;
	liquid_params lp_unit=liquid_params_unit_phi(lp0.dim,lp0.nup);
	double gamma_t, rel_error;
	double const gamma_max=1E40;
	int const max_iter=10000;
	int i1;
	char * fun="sk";
	printf("Asserting limits \n");
	gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k, lp0);
  	gamma_t=gamma_spherical_mono(Sg, lamk, lp0);
	if ( gamma_t >= gamma_max ) { 
		lp_inf=lp0; 
		lp_sup=lp1;
		gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k, lp1);
		gamma_t = gamma_spherical_mono(Sg, lamk, lp1);
		if ( gamma_t >= gamma_max ) {
			printf("Warning, no arrest found for the given limits, returning lp0 \n");
			return lp0;
		}
	}
	else{
		lp_inf=lp1;
		lp_sup=lp0;
		gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k, lp1);
		gamma_t = gamma_spherical_mono(Sg, lamk, lp1);
		if ( gamma_t < gamma_max ){
			printf("Warning, no fluid found for the given limits, returning lp0 \n");
			return lp0;
		}
	}
	/* Bisection until tol or max number of iterations is reached */
	printf("Limits asserted, starting bisection \n");
	i1=0; rel_error=1.0;
	while ( i1 < max_iter && rel_error > tol ){
		i1+=1;
		lp_test = liquid_params_scale( liquid_params_sum(lp_inf,lp_sup), 0.5) ;
		gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k, lp_test);
		gamma_t = gamma_spherical_mono(Sg, lamk, lp_test);
		if ( gamma_t >= gamma_max){ lp_inf = lp_test;}
		else{ lp_sup = lp_test;}
		rel_error = liquid_params_norm( liquid_params_dif( lp_unit, liquid_params_div( lp_sup, lp_inf ) ) );
		/*printf("%s%1.9e%s%1.9e%s%1.9e\n","phi_test=",lp_test.phi,"  gamma=",gamma_t,"  error=",rel_error); /* Printing iterations */
	}
	/*free(fun);free(lp_test.up);free(lp_dif.up);free(lp_unit.up);*/
	return lp_sup;
}

void intermediate_times_variables_writting(dynamics_save_options dso, dynamics_save_variables * dsv, int ini){
	int i1,i2;
	if ( ini == 0 ) {
		if ( dynamics_save_options_sum_tau_only( dso ) > 0 ) {
			fprintf(dsv->F_dyn, "# || 1 τ ||" );
			i1=1;
			if( dso.delta_zeta == 1 ) { 
				i1+=1; 
				fprintf(dsv->F_dyn, " %1.1d %s", i1, "Δζ(τ) ||" ); 
			}
			if( dso.Dt == 1 ) { 
				i1+=1; 
				fprintf(dsv->F_dyn, " %1.1d %s", i1, "D(τ) ||" ); 
			}
			if( dso.msd == 1 ) { 
				i1+=1; 
				fprintf(dsv->F_dyn, " %1.1d %s", i1, "w(τ) ||" ); 
			}
			if( dso.delta_eta == 1 ) { 
				i1+=1; 
				fprintf(dsv->F_dyn, " %1.1d %s", i1, "Δη(τ) ||" ); 
			}
			fprintf(dsv->F_dyn, "# \n");
		}

		if ( dso.Fc == 1 ){
			fprintf( dsv->F_Fc, "# || 1 τ || ");
			for ( i1=0; i1 < dsv->k->size; ++i1 ){
				fprintf( dsv->F_Fc, "%d %s %3.3f %s", i1 + 2 ,"k=", dsv->k->data[i1]," || " );
			}
			fprintf(dsv->F_Fc,"# \n");
		}
		if ( dso.Fs == 1 ){
			fprintf(dsv->F_Fs, "# || 1 τ || ");
			for ( i1=0; i1 < dsv->k->size; ++i1 ){
				fprintf(dsv->F_Fs, "%d %s %3.3f %s", i1 + 2, "k=", dsv->k->data[i1], " || " );
			}
			fprintf(dsv->F_Fs,"# \n");
		}
	}

	if ( dynamics_save_options_sum_tau_only( dso ) > 0 ) {
		for ( i1 = ini; i1 < dsv->tau->size; ++i1 ) {  
			fprintf(dsv->F_dyn,"%1.9e \t ", dsv->tau->data[i1] );
			if ( dso.delta_zeta == 1 ) { fprintf(dsv->F_dyn,"%1.9e \t ", dsv->delta_zeta_t->data[i1] ); }
			if ( dso.Dt == 1 ) { fprintf(dsv->F_dyn,"%1.9e \t ", dsv->Dt->data[i1] ); }
			if ( dso.msd == 1 ) { fprintf(dsv->F_dyn,"%1.9e \t ", dsv->msd->data[i1] ); }
			if ( dso.delta_eta == 1 ) { fprintf(dsv->F_dyn,"%1.9e \t ", dsv->delta_eta_t->data[i1] ); }
			fprintf(dsv->F_dyn," \n");
		}
	}

	if ( dso.Fc == 1 ) {
		for( i1=ini; i1 < dsv->Fc->size1; ++i1 ){
			fprintf(dsv->F_Fc, "%1.9e \t", dsv->tau->data[i1] );
			for( i2=0; i2 < dsv->Fc->size2; ++i2 ){
				fprintf(dsv->F_Fc, "%1.9e\t",gsl_matrix_get(dsv->Fc,i1,i2) );
			}
		fprintf(dsv->F_Fc, "\n");
		}

	}
	if ( dso.Fs == 1 ) {
		for( i1=ini; i1 < dsv->Fs->size1; ++i1 ){
			fprintf(dsv->F_Fs, "%1.9e \t", dsv->tau->data[i1] );
			for( i2=0; i2 < dsv->Fs->size2; ++i2 ){
				fprintf(dsv->F_Fs, "%1.9e\t",gsl_matrix_get(dsv->Fs, i1, i2) );
			}
		fprintf(dsv->F_Fs, "\n");
		}
	}
	return;
}


void gsl_vector_fc_non_ergo_param_mono_sph(gsl_vector ** fc,structure_grid Sg, gsl_vector * lamk, double gamma){
	gsl_vector_memcpy(fc[0],Sg.S);
	gsl_vector * dummy = gsl_vector_alloc(Sg.k->size);
	gsl_vector_memcpy(dummy,Sg.k);
	gsl_vector_mul(dummy,dummy);
 	gsl_vector_scale(dummy,gamma);
	gsl_vector_div(dummy,Sg.S);
	gsl_vector_div(dummy,lamk);
	gsl_vector_add_constant(dummy,1.0);
	gsl_vector_div(fc[0],dummy);
	gsl_vector_free(dummy);
	return;
}

void gsl_vector_fs_non_ergo_param_mono_sph(gsl_vector ** fs, structure_grid Sg, gsl_vector * lamk, double gamma){
	gsl_vector_set_all(fs[0],1.0);
	gsl_vector * dummy = gsl_vector_alloc(Sg.k->size);
	gsl_vector_memcpy(dummy,Sg.k);
	gsl_vector_mul(dummy,dummy);
 	gsl_vector_scale(dummy,gamma);
	gsl_vector_div(dummy,lamk);
	gsl_vector_add_constant(dummy,1.0);
	gsl_vector_div(fs[0],dummy);
	gsl_vector_free(dummy);
	return;
}

void writting_taua(dynamics_save_options dso, dynamics_save_variables * dsv,
structure_grid Sg){
	int i1,i2;
	if ( dso.tau_alpha == 1 ){
		fprintf( dsv->F_taua, "# || 1 k || 2 τα(k) ||");
		for( i1 = 0; i1 < Sg.k->size; ++i1 ){
			fprintf(dsv->F_taua,"%1.9e \t %1.9e \n", Sg.k->data[i1], dsv->tau_alpha->data[i1]);
		}
	}
	return;
}

/* Function that computes for the dynamics employing the SCGLE formalism iteratively  \
until a convergence in Dl is found for a monocomponent spherical system */
/* Inputs type							Variable name		Notes
	liquid_params										lp				Composed of double parameters
	dyn_params											dp				Composed of integers and doubles that helps in computing the dynamics
	structure_grid								 	Sg				Composed of 3 gsl_vector
	save_dyn_vars_ini								dyn				Composed of variables to be saved either to drive or as an output
	save_dyn_op_ini									op				Composed of options to know what to save or not
	dyn_scalar											ds				Composed of double variables which save scalar information of the dynamics
	itv_vars												itv				Composed of dynamics variables computed for an equally spaced time grid


Notes:
	- Currently only works for d=2 and d=3
	- γ is only computed for values < 1E40
	- If γ = 1E40 the value is expected to diverge
	- If γ = 1E99 no convergence was found in 10,000 iteration steps
	*/
void dynamics_spherical_mono( liquid_params lp, dynamics_parameters dp, structure_grid Sg, dynamics_save_variables * dsv, dynamics_save_options dso ){
	int nt = dp.it;
	int knp = Sg.k->size;
	/* Initialization of intermediate times variables for computing and writting */
	intermediate_times_variables itv = intermediate_times_variables_spherical_mono_ini(dp,Sg);
	gsl_vector_view delz_int;
	dynamics_parameters dpd=dp;
	int i1, i2, convergence;
	const double Dl_tol = dp.tol;
	double dtau_decim;
	double intdelz,Dl,Dl_test;
	double error;
	double gamma_val;
	FILE * Fdyn;
	FILE * Flam;
	/* Initialization of dynamics variables */
		/* Computing gamma */
	gamma_val = gamma_spherical_mono(Sg, itv.lambda_k, lp);
	if( dso.gamma == 1 ){
		dsv->gamma = gamma_val;
		printf("%s %1.9e \n", "Gamma=",dsv->gamma);
	}
	/* Setting initial times */
	short_times_dynamics_spherical_mono( lp, Sg, dp, dsv, &itv, dso);
	/* Initialization of intermediate times */
	intermediate_times_dynamics_spherical_mono(lp, Sg, dp, &itv, dsv, dso);
	intermediate_times_variables_writting(dso, dsv, 0);
	/* Computing long time diffusion coefficient */
	intdelz = int_Delz( itv.delta_z, 0.0, dpd.dtau,0,nt );
	Dl_test = Dl_sph_mono( intdelz );
	/* Decimation */
	convergence = 0;
	i1=0;
	dpd.st=dp.it/2;
	while (convergence == 0){
		i1=i1+1;
		save_half_intermediate_times_variables_spherical_mono( &itv,  &dpd, dso, dsv);
		intermediate_times_dynamics_spherical_mono( lp, Sg, dpd, &itv, dsv, dso );
		intermediate_times_variables_writting(dso, dsv, dp.it/2);
		/* Computing long time diffusion coefficient */
		intdelz = int_Delz( itv.delta_z, intdelz, dpd.dtau,dp.it/2,nt );
		Dl = Dl_sph_mono( intdelz );
		error = fabs( 1.0 - ( Dl_test / Dl ) );
		Dl_test = Dl;
		if ( i1 >= dp.minimum_decimations_number ){
			if ( error < Dl_tol || Dl < Dl_tol || i1 >= dp.maximum_decimations_number ){ convergence = 1; } /* Convergence condition */
		}
	}
	writting_taua(dso, dsv, Sg);
	if(dso.Dl==1){dsv->Dl=Dl;}
	printf("%s \t %1.9e\n","Dynamics Dl=", Dl);
	/* Freeing memory and closing files */
	free_intermediate_times_variables(&itv);
	return;
}

void dynamics_mono_spherical_standard_defined_structures( liquid_params lp, char * sys, char * approx, char * folder ) {
	/* Integration variables, structure grid definition and structure factor computation */
  	const double k0=0.0; const double k1=10.0; const double k2=40.96;
	const int np  = pow(2,8);
	const int nph = np / 2;
	char * fun="sk";
	gsl_integration_fixed_workspace * w1, * w2;
	const gsl_integration_fixed_type * int_T = gsl_integration_fixed_legendre;
	w1 = gsl_integration_fixed_alloc(int_T, nph, k0, k1, 0.0, 0.0);
	w2 = gsl_integration_fixed_alloc(int_T, nph, k1, k2, 0.0, 0.0);
	structure_grid Sg={NULL,NULL,NULL}; 
	structure_grid_ini(&Sg,np);
	gsl_vectors_composed_integration_space(w1,w2, &Sg.k, &Sg.kw);
	gsl_integration_fixed_free(w1); gsl_integration_fixed_free(w2);
	gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k, lp);
	/* Defining Dynamic variables */
	dynamics_parameters dp = dynamics_parameters_auto_ini(); /* dyn_params_ini_HD(); */
	/* Initializing variables needed for the dynamics */
		/* File handling variables */
	char * fname = (char *)malloc(400*sizeof(char));
	s_name_constructor(sys,approx, "dat",1,lp, &fname );
	s_grid_save_file( Sg, folder, "S_", fname );
		/**/
	printf("Initializing dynamics variables \n");
	dynamics_save_options dso = dynamics_save_options_auto_ini();
	int nkw = 5;
	gsl_vector * kw = gsl_vector_alloc(nkw);
	gsl_vector * Sw = gsl_vector_alloc(nkw);
	kw->data[0] = 2.0; 
	kw->data[1] = 4.0; 
	kw->data[2] = 5.0;
	kw->data[3] = 2.0 * M_PI; 
	kw->data[4] = 7.18;
	gsl_vector_s_function_selector_mono_sph(Sw, sys, approx, fun, kw, lp); 
	dynamics_save_variables dyn;
	dynamics_save_variables_spherical_mono_ini(&dyn, dso, kw, Sw, Sg.k->size, dp, folder, "", fname );
	/* Computing the dynamics */
	printf("%s \n", "Starting computation of dynamics");
	dynamics_spherical_mono( lp, dp, Sg, &dyn, dso);
	printf("%s \t %1.9e \n", "Dl=", dyn.Dl);
	printf("%s \t %1.9e \n", "gamma=", dyn.gamma);
	/* Freeing memory */
	free(fname);
	structure_grid_free(&Sg);
	free_close_dynamics_save_variables( dso, &dyn );
	gsl_vector_free(kw);
	gsl_vector_free(Sw); 
	return;
}
