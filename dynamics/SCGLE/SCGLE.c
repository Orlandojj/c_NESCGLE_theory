#include "SCGLE.h"

/*Function that computes the wave vector dependent part of the memory function
	lambda(k,kc)
	"k" -> wave vector magnitude
	"kc" ->  Adjustable parameter, found that for HS systems kc ≈ 2 π * 1.305*/
double
lambda_sph_mono( const double k, const double kc ) {
		double l;
		l = k / kc;
		l = l * l;
		l = 1.0 + l;
		l = 1.0 / l;
		return l;
	}

void lambda_sph_mono_gsl_v (gsl_vector ** lamk, const gsl_vector * k, const double kc ){
	;
	if( k!=NULL ) {
		int knp=k->size;
		int i1;
		for (i1=0; i1<knp; ++i1){gsl_vector_set(*lamk,i1, lambda_sph_mono( k->data[i1], kc ));}
	}
	return;
}

dyn_params
dyn_params_ini(int st, int it, double dtau, double kc, double D0 ){
	dyn_params dp={st,it,dtau, kc, D0};
	return dp;
}

dyn_params
dyn_params_ini_std(){
	dyn_params dp={10,256,1e-7,2.0*M_PI*1.305,1.0,1E-6};
	return dp;
}

dyn_params
dyn_params_ini_HD(){
	dyn_params dp={10,512,1e-7,4.0*M_PI,1.0,1E-6};
	return dp;
}

dyn_scalar
dyn_scalar_ini(){
	dyn_scalar ds = {0,0,0,0};
	return ds;
}

save_dyn_op
save_dyn_op_ini(){
	save_dyn_op so;
	so.save_delz=0;
	so.save_deleta=0;
	so.save_F=0;
	so.save_gamma=0;
	so.save_Dl=1;
	so.write_delz=1;
	so.write_taua=1;
	so.write_deleta=1;
	so.write_F=1;
	return so;
}

save_dyn_op
no_save_dyn_op_ini(){
	save_dyn_op so;
	so.save_delz=0;
	so.save_deleta=0;
	so.save_F=0;
	so.save_gamma=0;
	so.save_Dl=1;
	so.write_delz=0;
	so.write_taua=0;
	so.write_deleta=0;
	so.write_F=0;
	return so;
}

void
save_dyn_vars_ini( save_dyn_vars * sdv, const int knp, const char * folder, const char * prefix, const char * suffix ){
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
	sdv->k=gsl_vector_alloc(knp);
	sdv->S=gsl_vector_alloc(knp);
	sdv->tau=NULL;
	sdv->delz=NULL;
	sdv->dele=NULL;
	sdv->Fc=NULL;
	sdv->Fs=NULL;
	sdv->F_dyn  = fopen(dyn_name, "w");
	sdv->F_Fc   = fopen(fc_name, "w");
	sdv->F_Fs   = fopen(fs_name, "w");
	sdv->F_taua = fopen(taua_name, "w");
	free(dyn_name);
	free(taua_name);
	free(fc_name);
	free(fs_name);
	return;
}

void save_dyn_vars_free( save_dyn_vars * sdv )
{
	if( sdv->tau !=NULL ) { gsl_vector_free(sdv->tau); }
	if( sdv->delz !=NULL ) { gsl_vector_free(sdv->delz); }
	if( sdv->dele !=NULL ) { gsl_vector_free(sdv->dele); }
	if( sdv->Fc !=NULL ) { gsl_matrix_free(sdv->Fc); }
	if( sdv->Fs !=NULL ) { gsl_matrix_free(sdv->Fs); }
	if( sdv->k !=NULL ) { gsl_vector_free(sdv->k); }
	if( sdv->S !=NULL ) { gsl_vector_free(sdv->S); }
}

void
save_dyn_vars_close( save_dyn_op sdo, save_dyn_vars * sdv){
	if ( sdo.write_deleta == 1 || sdo.write_delz == 1 ) {fclose(sdv->F_dyn);}
	if ( sdo.write_F == 1 ) {fclose(sdv->F_Fc); fclose(sdv->F_Fs);}
	if ( sdo.write_taua == 1 ) {fclose(sdv->F_taua);}
	return;
}

it_vars 
it_vars_ini( const dyn_params dp, const structure_grid Sg, const structure_grid Sg_w){
	int knp = Sg.k->size;
	int knp_w; if ( Sg_w.k != NULL ) { knp_w = Sg_w.k->size; }
	int nt=dp.it;
	gsl_vector * tau = gsl_vector_alloc(nt);
	gsl_vector * idelz = gsl_vector_alloc(nt);
	gsl_vector * msd = gsl_vector_calloc(nt); /* Not necesary so initializing it at 0 */
	gsl_vector * Dt = gsl_vector_calloc(nt); /* Not necesary so initializing it at 0 */
	gsl_vector * idele = gsl_vector_calloc(nt);/* Not necesary so initializing it at 0 */
	gsl_vector * tau_alpha = gsl_vector_calloc(knp); /* Not necesary so initializing it at 0 */
	gsl_vector * lamk = gsl_vector_alloc(knp); lambda_sph_mono_gsl_v( &lamk, Sg.k, dp.kc );
	gsl_matrix * iFc = gsl_matrix_alloc(nt, knp);
	gsl_matrix * iFs = gsl_matrix_alloc(nt, knp);
	gsl_vector * lamk_w = NULL; if ( Sg_w.k != NULL ) { lamk_w = gsl_vector_alloc(knp_w); lambda_sph_mono_gsl_v(&lamk_w, Sg_w.k, dp.kc ); }
	gsl_matrix * iFc_w = NULL; if ( Sg_w.k != NULL ) { iFc_w = gsl_matrix_alloc(nt, knp_w); }
	gsl_matrix * iFs_w = NULL; if ( Sg_w.k != NULL ) { iFs_w = gsl_matrix_alloc(nt, knp_w); }
	int i1;
	double lam_val;
	for(i1=0; i1<nt ; ++i1){ gsl_vector_set(tau,i1,( i1 + 1 ) * dp.dtau); }
	gsl_vector_set_all(tau_alpha,dp.dtau);
	it_vars itv={ tau, idelz, msd, Dt, idele, tau_alpha, lamk, iFc, iFs, lamk_w, iFc_w, iFs_w };
	return itv;
}


void
it_vars_free(it_vars * itv){
	if(itv->tau != NULL) 		{ gsl_vector_free (itv->tau); }
	if(itv->delz != NULL) 		{ gsl_vector_free (itv->delz); }
	if(itv->msd !=NULL) 		{ gsl_vector_free (itv->msd); }
	if(itv->Dt !=NULL) 			{ gsl_vector_free (itv->Dt); }
	if(itv->tau_alpha !=NULL) 	{ gsl_vector_free (itv->tau_alpha); }
	if(itv->dele !=NULL) 		{ gsl_vector_free (itv->dele); }
	if(itv->lamk !=NULL) 		{ gsl_vector_free (itv->lamk); }
	if(itv->Fc !=NULL) 			{ gsl_matrix_free (itv->Fc); }
	if(itv->Fs !=NULL) 			{ gsl_matrix_free (itv->Fs); }
	if(itv->lamk_w !=NULL)		{ gsl_vector_free (itv->lamk_w); }
	if(itv->Fc_w !=NULL) 		{ gsl_matrix_free (itv->Fc_w); }
	if(itv->Fs_w !=NULL) 		{ gsl_matrix_free (itv->Fs_w); }
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
del_z_sph_mono( const liquid_params lp, const structure_grid Sg, const gsl_vector * Fc, const gsl_vector * Fs ){
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

double
Dl_sph_mono( double int_delz ){
	double Dl;
	Dl = 1.0 + int_delz;
	Dl = 1.0 / Dl;
	return Dl;
}

/* Function that computes for the correlation time dependent diffusion coefficient
for the it- time through the convolution D(t)= 1 - ∫₀ᵗ	 Δζ(t-t')D(t')dt' */

double Dt(int it, dyn_params dp, gsl_vector * delz, gsl_vector * Dt ){
	int i1;
	double DDz, sum, dum1;
	sum = 0.0;
	for (i1=0; i1 < it; ++i1){
		DDz = delz->data[it-i1] - delz->data[it-i1-1];
		sum = sum + ( Dt->data[i1] * DDz );
	}
	dum1 = ( Dt->data[it-1] / dp.dtau ) - sum;
	dum1 = dum1 / ( delz->data[0] + ( 1.0 / dp.dtau ) );
	return dum1;
}
/* Function that computes for the mean-squared displacement
for the it- time through the convolution w(t)= t - ∫₀ᵗ	[ d ( w(t-t') ) / dt ] [ Δζ(t')dt' ]  */
double
msdt(int it, dyn_params dp, gsl_vector * delz, gsl_vector * msd){
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
double
ini_Fc_sph_mono( double k, double t, double S, double D0 ){
	double f;
	f =  exp (- D0 * t * k * k / S) * S;
	return f;
}

/*Function that computes the small time limit of the self part of the intermediate scattering function
	Fs(k,t)
	"k" -> wave vector magnitude
	"t" -> correlation time value */
double
ini_Fs_sph_mono( double k, double t, double D0 ){
	double f;
	f = exp (- D0 * t * k * k ) ;
	return f;
}

double
tau_alpha_sph_mono( double Fs_pre, double tau_pre, double Fs, double tau, double D0, double taua_pre ){
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
void
small_t_dynamics_sph_mono(liquid_params lp, structure_grid Sg, dyn_params dp,
it_vars * itv, save_dyn_op sdo, structure_grid Sg_w ){
	int knp = (Sg.k)->size;
	int knp_w; if ( sdo.write_F == 1 ) { knp_w = (Sg_w.k)->size; }
	int i1,i2;
	double t,k_val,S_val,z;
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
			if (i2>0) {
				itv->tau_alpha->data[i1] = tau_alpha_sph_mono( gsl_matrix_get(itv->Fs,i2-1,i1),
				itv->tau->data[i2-1], gsl_matrix_get(itv->Fs,i2,i1), itv->tau->data[i2], dp.D0,
				itv->tau_alpha->data[i1] ) ;
			}
		}
		if ( sdo.write_F == 1 ) {
			for( i1 = 0; i1 < knp_w; ++i1 ){
				k_val = gsl_vector_get (Sg_w.k, i1);
				S_val = gsl_vector_get (Sg_w.S, i1);
				gsl_matrix_set ( itv->Fc_w, i2, i1, ini_Fc_sph_mono( k_val, t, S_val, dp.D0 ) );
				gsl_matrix_set ( itv->Fs_w, i2, i1, ini_Fs_sph_mono( k_val, t, dp.D0 ) );
			}
		}
		Fc_v = gsl_matrix_row( itv->Fc, i2);
		Fs_v = gsl_matrix_row( itv->Fs, i2);
		z = del_z_sph_mono( lp, Sg,&Fc_v.vector,&Fs_v.vector );
		gsl_vector_set ( itv->delz, i2, z );
		sumDelz = sumDelz + z;
		dum1 = dp.D0 * ( 1.0 - dp.dtau * sumDelz );
		gsl_vector_set ( itv->Dt, i2, dum1 );
	}
	itv->msd->data[0] = itv->Dt->data[0] * dp.dtau;
	for( i1 = 1; i1< dp.st; ++i1 ){
		itv->msd->data[i1] = itv->msd->data[i1-1] + itv->Dt->data[i1] * dp.dtau;
	}
}



/* Function that computes for auxiliary vectors as and ac employed in medium_t_dynamics_sph_mono */
void
alphas( structure_grid Sg, dyn_params dp, gsl_vector * lamk,
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
void
medium_t_dynamics_sph_mono_writting(liquid_params lp, structure_grid Sg, dyn_params dp, it_vars * itv){
		int knp = (Sg.k)->size;
		int i1,i2,i3,i4;
		double t,k_val,S_val,z;
		double dum1,deldelz_val,delz_val;
		double delz_test;
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
		alphas( Sg, dp, itv->lamk_w, itv->delz->data[0], &ac, &as);
		gsl_vector_memcpy(lamc,ac); gsl_vector_mul(lamc, itv->lamk_w); /* lamc = ac * lamk */
		gsl_vector_memcpy(lams,as); gsl_vector_mul(lams, itv->lamk_w); /* lams = as * lamk */
		Fc1 = gsl_matrix_row( itv->Fc_w, 0); Fs1 = gsl_matrix_row( itv->Fs_w, 0); /* Fc1=Fc(0); Fs1=Fs(0) */
		gsl_vector_set_all(dFs1,1.0); gsl_vector_sub(dFs1,&Fs1.vector); /* dFs1 = 1-Fs(0) */
		gsl_vector_memcpy(dFc1,Sg.S); gsl_vector_sub(dFc1,&Fc1.vector); /* dFc1 = S -Fc(0) */
			/* deldelz(i) = delz(i) - delz(i-1) */
		gsl_vector_set_zero(deldelz);
		for(i1=1; i1<dp.it; ++i1){
			dum1 = gsl_vector_get(itv->delz,i1) - gsl_vector_get(itv->delz,i1-1);
			gsl_vector_set(deldelz,i1,dum1);
		}
		/* Computing the dynamic variables for new times */
		for(i1=dp.st; i1<dp.it; ++i1){
			gsl_vector_set_zero(Fcdum); /* Fcdum=0 */
			gsl_vector_set_zero(Fsdum); /* Fsdum=0 */
			for(i2=1; i2<i1; ++i2){
				Fc_v=gsl_matrix_row( itv->Fc_w, i2 );
				Fs_v=gsl_matrix_row( itv->Fs_w, i2 );
				deldelz_val = deldelz->data[i1-i2];
				for(i3=0; i3<knp; ++i3){
					/*  */
					Fc_val = gsl_vector_get(&Fc_v.vector,i3); /* Fc(i2) */
					Fs_val = gsl_vector_get(&Fs_v.vector,i3); /* Fs(i2) */
					Fcdum->data[i3] = Fcdum->data[i3] - ( Fc_val * deldelz_val ) ;
					Fsdum->data[i3] = Fsdum->data[i3] - ( Fs_val * deldelz_val ) ;
				}
			}
			Fc_v = gsl_matrix_row( itv->Fc_w, i1-1 ); Fs_v=gsl_matrix_row( itv->Fs_w, i1-1 );
			delz_val = itv->delz->data[i1-1];
			for(i3=0; i3<knp; ++i3){
				Fc_val = gsl_vector_get(&Fc1.vector,i3); /* F(0) */
				dum_val1 = ( Fcdum->data[i3] + ( Fc_val * delz_val ) ) * itv->lamk_w->data[i3];
				dum_val1 = dum_val1  + gsl_vector_get(&Fc_v.vector,i3) / dp.dtau;
				dum_val1 = dum_val1 * ac->data[i3];
				gsl_vector_set( Fcdum, i3, dum_val1 );
				Fcdum2->data[i3] = Fcdum->data[i3] + ( dFc1->data[i3] * itv->delz->data[i1] * lamc->data[i3] );
				gsl_matrix_set(itv->Fc_w,i1,i3, gsl_vector_get(Fcdum2,i3) );

				Fs_val = gsl_vector_get(&Fs1.vector,i3); /* F(0) */
				dum_val1 = ( Fsdum->data[i3] + ( Fs_val * delz_val)  ) * itv->lamk_w->data[i3];
				dum_val1 = dum_val1  + gsl_vector_get(&Fs_v.vector,i3) / dp.dtau;
				dum_val1 = dum_val1 * as->data[i3];
				gsl_vector_set( Fsdum, i3, dum_val1 );
				Fsdum2->data[i3] = Fsdum->data[i3] + ( dFs1->data[i3] * itv->delz->data[i1] * lams->data[i3] );
				gsl_matrix_set(itv->Fs_w,i1,i3, gsl_vector_get(Fsdum2,i3) );
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
void
medium_t_dynamics_sph_mono( liquid_params lp, structure_grid Sg, dyn_params dp,
it_vars * itv, save_dyn_op sdo, structure_grid Sg_w ){
	const double tol=1e-10;
	int knp = (Sg.k)->size;
	int i1,i2,i3,i4;
	double t,k_val,S_val,z,error;
	double dum1,deldelz_val,delz_val,m,c;
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
	alphas( Sg, dp, itv->lamk, itv->delz->data[0], &ac, &as);
	gsl_vector_memcpy(lamc,ac); gsl_vector_mul(lamc, itv->lamk); /* lamc = ac * lamk */
	gsl_vector_memcpy(lams,as); gsl_vector_mul(lams, itv->lamk); /* lams = as * lamk */
	Fc1 = gsl_matrix_row( itv->Fc, 0); Fs1 = gsl_matrix_row( itv->Fs, 0); /* Fc1=Fc(0); Fs1=Fs(0) */
	gsl_vector_set_all(dFs1,1.0); gsl_vector_sub(dFs1,&Fs1.vector); /* dFs1 = 1-Fs(0) */
	gsl_vector_memcpy(dFc1,Sg.S); gsl_vector_sub(dFc1,&Fc1.vector); /* dFc1 = S -Fc(0) */
		/* deldelz(i) = delz(i) - delz(i-1) */
	gsl_vector_set_zero(deldelz);
	for(i1=1; i1<dp.st; ++i1){
		dum1 = gsl_vector_get(itv->delz,i1) - gsl_vector_get(itv->delz,i1-1);
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
		delz_val = itv->delz->data[i1-1];
		for(i3=0; i3<knp; ++i3){
			Fc_val = gsl_vector_get(&Fc1.vector,i3); /* F(0) */
			dum_val1 = ( Fcdum->data[i3] + ( Fc_val * delz_val ) ) * itv->lamk->data[i3];
			dum_val1 = dum_val1  + gsl_vector_get(&Fc_v.vector,i3) / dp.dtau;
			dum_val1 = dum_val1 * ac->data[i3];
			gsl_vector_set( Fcdum, i3, dum_val1 );

			Fs_val = gsl_vector_get(&Fs1.vector,i3); /* F(0) */
			dum_val1 = ( Fsdum->data[i3] + ( Fs_val * delz_val)  ) * itv->lamk->data[i3];
			dum_val1 = dum_val1  + gsl_vector_get(&Fs_v.vector,i3) / dp.dtau;
			dum_val1 = dum_val1 * as->data[i3];
			gsl_vector_set(Fsdum, i3, dum_val1);
		}

		/*if (dp.st>200 ){
		/*if (i1 == dp.it-1 ){*/
			/*FILE * filetest=fopen("./F_test.dat", "w");
			for (i2 = 0; i2<knp; ++i2){
				fprintf(filetest, "%1.9e \t %1.9e \t %1.9e \t %1.9e \t %1.9e \n",Sg.k->data[i2],
				lamc->data[i2], lams->data[i2], ac->data[i2], as->data[i2]);
			}
			printf("%1.9e \n", dp.dtau);
			fclose(filetest);
			exit(1);
		}*/
		/* Iteration */
		error = 1.0; /* setting a new error value to initiate the iteration */
		m = deldelz->data[i1-1]/dp.dtau;
		c = itv->delz->data[i1-1] - m * itv->tau->data[i1-1];
		delz_test = m * itv->tau->data[i1] + c; /* gsl_vector_get(itv->delz,i1-1) - gsl_vector_get(itv->delz,i1-1); /* setting a test value for Δζ(i1) */
		while ( error > tol ) {
			/* Computing a new value for Fc(i1) and Fs(i1) in terms of the test value of Δζ(i1) */
			for (i3=0; i3<knp; ++i3){
				Fcdum2->data[i3] = Fcdum->data[i3] + ( dFc1->data[i3] * delz_test * lamc->data[i3] );
				Fsdum2->data[i3] = Fsdum->data[i3] + ( dFs1->data[i3] * delz_test * lams->data[i3] );
			}
			/* Computing Δζ with the computed values for Fc(i1) and and Fs(i1) */
			delz_val = del_z_sph_mono( lp, Sg, Fcdum2, Fsdum2 );

			/* Computing the relative error between the test Δζ value and the computed Δζ  */
			error = fabs(1.0 - (delz_test / delz_val ));
			/* Updating a new value to test */
			delz_test = delz_val;
		}
		/* Saving the data for the i1-time */
		gsl_vector_set(itv->delz,i1,delz_val);
		for ( i2=0; i2 < knp; ++i2){
			gsl_matrix_set(itv->Fc,i1,i2, gsl_vector_get(Fcdum2,i2) );
			gsl_matrix_set(itv->Fs,i1,i2, gsl_vector_get(Fsdum2,i2) );

			itv->tau_alpha->data[i2] = tau_alpha_sph_mono( gsl_matrix_get(itv->Fs,i1-1,i2),
			itv->tau->data[i1-1], gsl_matrix_get(itv->Fs,i1,i2), itv->tau->data[i1], dp.D0,
			itv->tau_alpha->data[i2] ) ;
		}

		deldelz->data[i1] = gsl_vector_get(itv->delz,i1) - gsl_vector_get(itv->delz,i1-1);
		itv->Dt->data[i1]  = Dt(i1, dp, itv->delz, itv->Dt );
		itv->msd->data[i1] = msdt(i1, dp, itv->delz, itv->msd);
		/*if(dp.st>200 ){
			printf("%1.9e\n", itv->delz->data[i1]-1041.4876421622946);
			itv->delz->data[i1]= 1041.4876421622946;
			itv->msd->data[i1] = msdt(i1, dp, itv->delz, itv->msd);
			exit(1);
		}*/
	}
	if( sdo.write_F == 1 ){
		medium_t_dynamics_sph_mono_writting(lp, Sg_w, dp, itv);
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



void
save_half_itv_sph_mono( it_vars * itv, dyn_params * dp, save_dyn_op sdo ){
	int i1,i2,isave;
	int nt =  dp->it;
	const int nth = nt / 2;
	const int knp = itv->Fc->size2;
	int knp_w; if ( sdo.write_F==1) {knp_w = itv->Fc_w->size2;}
	dp->dtau = 2.0 * dp->dtau;
	gsl_vector_scale(itv->tau,2.0);
	for(i1=0; i1<nth; ++i1){
		isave=1+i1*2;
		gsl_vector_set(itv->delz,i1,0.5 * (gsl_vector_get(itv->delz,isave)+gsl_vector_get(itv->delz,isave-1)));
		gsl_vector_set(itv->msd,i1,0.5 * (gsl_vector_get(itv->msd,isave)+gsl_vector_get(itv->msd,isave-1)));
		gsl_vector_set(itv->Dt,i1,0.5* (gsl_vector_get(itv->Dt,isave)+gsl_vector_get(itv->Dt,isave-1)));
		gsl_vector_set(itv->dele,i1,0.5*(gsl_vector_get(itv->dele,isave)+gsl_vector_get(itv->dele,isave)));
		for(i2=0; i2<knp; ++i2){
			gsl_matrix_set(itv->Fc,i1,i2,0.5*(gsl_matrix_get(itv->Fc,isave,i2)+gsl_matrix_get(itv->Fc,isave-1,i2)));
			gsl_matrix_set(itv->Fs,i1,i2,0.5*(gsl_matrix_get(itv->Fs,isave,i2)+gsl_matrix_get(itv->Fs,isave-1,i2)));
		}
		if ( sdo.write_F==1) {
			for(i2=0; i2<knp_w; ++i2){
				gsl_matrix_set(itv->Fc_w,i1,i2,0.5*(gsl_matrix_get(itv->Fc_w,isave,i2)+gsl_matrix_get(itv->Fc_w,isave-1,i2)));
				gsl_matrix_set(itv->Fs_w,i1,i2,0.5*(gsl_matrix_get(itv->Fs_w,isave,i2)+gsl_matrix_get(itv->Fs_w,isave-1,i2)));
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
double 
gamma_sph_mono(structure_grid Sg, gsl_vector * lamk, liquid_params lp){
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
  	gamma_t=gamma_sph_mono(Sg, lamk, lp0);
	if ( gamma_t >= gamma_max ) { 
		lp_inf=lp0; 
		lp_sup=lp1;
		gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k, lp1);
		gamma_t = gamma_sph_mono(Sg, lamk, lp1);
		if ( gamma_t >= gamma_max ) {
			printf("Warning, no arrest found for the given limits, returning lp0 \n");
			return lp0;
		}
	}
	else{
		lp_inf=lp1;
		lp_sup=lp0;
		gsl_vector_s_function_selector_mono_sph(Sg.S, sys, approx, fun, Sg.k, lp1);
		gamma_t = gamma_sph_mono(Sg, lamk, lp1);
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
		gamma_t = gamma_sph_mono(Sg, lamk, lp_test);
		if ( gamma_t >= gamma_max){ lp_inf = lp_test;}
		else{ lp_sup = lp_test;}
		rel_error = liquid_params_norm( liquid_params_dif( lp_unit, liquid_params_div( lp_sup, lp_inf ) ) );
		/*printf("%s%1.9e%s%1.9e%s%1.9e\n","phi_test=",lp_test.phi,"  gamma=",gamma_t,"  error=",rel_error); /* Printing iterations */
	}
	/*free(fun);free(lp_test.up);free(lp_dif.up);free(lp_unit.up);*/
	return lp_sup;
}

void
writting_medium_times_itvs(save_dyn_op sdo, save_dyn_vars * sdv,
structure_grid Sg_w, it_vars itv, int ini){
	int i1,i2;
	if(sdo.write_delz==1){
		if (ini==0){
			fprintf(sdv->F_dyn, "# || 1 τ || 2 Δζ(τ) || 3 D(τ) || 4 w(τ) || 5 Δη(τ)  || # \n" );
		}

		for( i1 = ini; i1 < itv.tau->size; ++i1 ){
			fprintf(sdv->F_dyn,"%1.9e \t %1.9e \t %1.9e \t %1.9e \t %1.9e \n",
			itv.tau->data[i1], itv.delz->data[i1], itv.Dt->data[i1],
			itv.msd->data[i1], itv.dele->data[i1] );
		}
	}
	if(sdo.write_F==1){
		if (ini==0){
			fprintf(sdv->F_Fc, "# || 1 τ || ");
			fprintf(sdv->F_Fs, "# || 1 τ || ");
			for (i1=0;i1<Sg_w.k->size;++i1){
				fprintf(sdv->F_Fc,"%d %s %3.3f %s" , i1+2,"k=",Sg_w.k->data[i1]," || " );
				fprintf(sdv->F_Fs,"%d %s %3.3f %s" , i1+2,"k=",Sg_w.k->data[i1]," || " );
			}
			fprintf(sdv->F_Fc,"# \n");
			fprintf(sdv->F_Fs,"# \n");
		}
		for( i1=ini; i1 < itv.Fc_w->size1; ++i1 ){
			fprintf(sdv->F_Fc, "%1.9e \t", itv.tau->data[i1] );
			fprintf(sdv->F_Fs, "%1.9e \t", itv.tau->data[i1] );
			for( i2=0; i2 < itv.Fc_w->size2; ++i2 ){
				fprintf(sdv->F_Fc, "%1.9e\t",gsl_matrix_get(itv.Fc_w,i1,i2) );
				fprintf(sdv->F_Fs, "%1.9e\t",gsl_matrix_get(itv.Fs_w,i1,i2) );
			}
		fprintf(sdv->F_Fc, "\n");
		fprintf(sdv->F_Fs, "\n");
		}
	}
	return;
}


void
gsl_vector_fc_non_ergo_param_mono_sph(gsl_vector ** fc,structure_grid Sg, gsl_vector * lamk, double gamma){
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

void
gsl_vector_fs_non_ergo_param_mono_sph(gsl_vector ** fs, structure_grid Sg, gsl_vector * lamk, double gamma){
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

void
writting_taua(save_dyn_op sdo, save_dyn_vars * sdv,
structure_grid Sg, it_vars itv){
	int i1,i2;
	if(sdo.write_taua==1){
		fprintf(sdv->F_taua,"# || 1 k || 2 τα(k) ||");
		for( i1 = 0; i1 < Sg.k->size; ++i1 ){
			fprintf(sdv->F_taua,"%1.9e \t %1.9e \n", Sg.k->data[i1], itv.tau_alpha->data[i1]);
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
void
dynamics_sph_mono( liquid_params lp, dyn_params dp, structure_grid Sg,
save_dyn_vars * dyn, save_dyn_op sdo, dyn_scalar * ds ){
	int nt = dp.it;
	int knp = Sg.k->size;
	structure_grid Sg_w={NULL,NULL,NULL};
	/* Initialization of writting grid */
	if (sdo.write_F == 1 ){
		int knp_w = dyn->k->size;
		structure_grid_ini(&Sg_w,knp_w);
		gsl_vector_memcpy(Sg_w.k,dyn->k);
		gsl_vector_memcpy(Sg_w.S,dyn->S);
	} 
	/* Initialization of intermediate times variables for computing and writting */
	it_vars itv = it_vars_ini(dp,Sg,Sg_w);
	gsl_vector_view delz_int;
	dyn_params dpd=dp;
	int i1, i2, convergence;
	const double Dl_tol = dp.tol;
	double dtau_decim;
	double intdelz,Dl,Dl_test;
	double error;
	FILE * Fdyn;
	FILE * Flam;
	/* Initialization of dynamics variables */
		/* Computing gamma */
	if(sdo.save_gamma==1){
		ds->gamma = gamma_sph_mono(Sg, itv.lamk, lp);
		printf("%s %1.9e \n", "Gamma=",ds->gamma);
	}
	/* Setting initial times */
	small_t_dynamics_sph_mono( lp, Sg, dp, &itv, sdo, Sg_w);
	/* Initialization of intermediate times */
	medium_t_dynamics_sph_mono(lp, Sg, dp, &itv, sdo, Sg_w);
	writting_medium_times_itvs(sdo, dyn, Sg_w, itv, 0);
	/* Computing long time diffusion coefficient */
	intdelz = int_Delz( itv.delz, 0.0, dpd.dtau,0,nt );
	Dl_test = Dl_sph_mono( intdelz );
	/* Decimation */
	convergence = 0;
	i1=0;
	dpd.st=dp.it/2;
	while (convergence == 0){
		i1=i1+1;
		save_half_itv_sph_mono( &itv,  &dpd, sdo);
		medium_t_dynamics_sph_mono( lp, Sg, dpd, &itv, sdo, Sg_w );
		writting_medium_times_itvs(sdo, dyn, Sg_w, itv, dp.it/2);
		/* Computing long time diffusion coefficient */
		intdelz = int_Delz( itv.delz, intdelz, dpd.dtau,dp.it/2,nt );
		Dl = Dl_sph_mono( intdelz );
		error = fabs( 1.0 - ( Dl_test / Dl ) );
		Dl_test = Dl;
		if ( error < Dl_tol || Dl < Dl_tol ){ convergence = 1; } /* Convergence condition */
	}
	writting_taua(sdo, dyn, Sg, itv);
	if(sdo.save_Dl==1){ds->Dl=Dl;}
	printf("%s \t %1.9e\n","Dynamics Dl=", Dl);
	/* Freeing memory and closing files */
	it_vars_free(&itv);
	structure_grid_free(&Sg_w);
	save_dyn_vars_close( sdo, dyn);
	return;
}
