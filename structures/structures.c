#include "structures.h"

void
structure_grid_ini(structure_grid * Sg, int knp){
	Sg->k  = gsl_vector_alloc(knp);
	Sg->kw = gsl_vector_alloc(knp);
	Sg->S  = gsl_vector_alloc(knp);
	return;
}

void
structure_grid_memcpy(structure_grid dest, const structure_grid source){
	gsl_vector_memcpy(dest.k,source.k);
	gsl_vector_memcpy(dest.kw,source.kw);
	gsl_vector_memcpy(dest.S,source.S);
	return;
}

void
structure_grid_free(structure_grid * Sg){
	if ( Sg->k != NULL ) { gsl_vector_free(Sg->k); }
	if ( Sg->kw != NULL ) { gsl_vector_free(Sg->kw); }
	if ( Sg->S != NULL ) { gsl_vector_free(Sg->S); }
	return;
}

liquid_params 
liquid_params_ini_phi(double phi, double dim, int nup){
	double * up=NULL;
	if (nup>0) {
		up = ( double * ) calloc( nup, sizeof(double) );
		liquid_params lp={(2.0 * dim * phi) / M_PI,phi,dim,0.0,up,nup};
		return lp;
		}
	else{
		liquid_params lp={(2.0 * dim * phi) / M_PI,phi,dim,0.0,up,0};
		return lp;
	}
}

void liquid_params_free(liquid_params * lp){
	if (lp->up!=NULL){free(lp->up);}
}

liquid_params 
liquid_params_unit_phi(double dim, int nup){
	double * up;
	int i1;
	if (nup>0) {
		up = ( double * ) calloc( nup, sizeof(double) );
		for (i1=0; i1<nup; i1++) { up[i1]=1.0; }
		liquid_params lp={(2.0 * dim ) / M_PI,1.0,dim,1.0,up,nup};
		return lp;
		}
	else{
		liquid_params lp={(2.0 * dim ) / M_PI,1.0,dim,1.0,up,0};
		return lp;
	}
}

liquid_params 
liquid_params_sum(liquid_params lp1, liquid_params lp2){
	assert( lp1.dim == lp2.dim );
	assert( lp1.nup == lp2.nup );
	liquid_params sum = lp1;
	sum.rho += lp2.rho;
	sum.phi += lp2.phi;
	sum.Tem += lp2.Tem;
	int i1; for(i1=0; i1<sum.nup ; ++i1){sum.up[i1] += lp2.up[i1];}
	return sum;
}

liquid_params 
liquid_params_dif(liquid_params lp1, liquid_params lp2){
	assert( lp1.dim == lp2.dim );
	assert( lp1.nup == lp2.nup );
	liquid_params sum = lp1;
	sum.rho -= lp2.rho;
	sum.phi -= lp2.phi;
	sum.Tem -= lp2.Tem;
	int i1; for(i1=0; i1<sum.nup ; ++i1){sum.up[i1] += lp2.up[i1];}
	return sum;
}

liquid_params 
liquid_params_scale(liquid_params lp, double scale){
	liquid_params lps = lp;
	lps.rho *= scale;
	lps.phi *= scale;
	lps.Tem *= scale;
	int i1; for(i1=0; i1<lp.nup ; ++i1){lps.up[i1] *= scale;}
	return lps;
}

liquid_params 
liquid_params_div(liquid_params lp1, liquid_params lp2){
	assert( lp1.dim == lp2.dim );
	assert( lp1.nup == lp2.nup );
	int i1;
	liquid_params lpdiv = lp1;
	lpdiv.phi *= 1.0 / lp2.phi;
	lpdiv.rho = lpdiv.phi;
	if ( lp2.Tem != 0.0 ) { lpdiv.Tem *= 1.0 / lp2.Tem; } 
	else if ( lp1.Tem == 0.0 ) { lpdiv.Tem = 1.0; }
	else { lpdiv.Tem = 1.0E99; }
	for(i1=0; i1<lp1.nup ; ++i1){
		if ( lp2.up[i1] != 0.0 ) { lpdiv.up[i1] *= 1.0 / lp2.up[i1]; } 
		else if ( lp1.up[i1] == 0.0 ) { lpdiv.up[i1] = 1.0; }
		else { lpdiv.up[i1] = 1.0E99; }
		}
	return lpdiv;
}

double
liquid_params_norm(liquid_params lp){
	double sum;
	int i1;
	sum = ( lp.phi * lp.phi ) + ( lp.Tem * lp.Tem );
	for(i1=0; i1<lp.nup ; ++i1){sum+=lp.up[i1];}
	sum = sqrt(sum);
	return sum;
}

void
s_name_constructor(const char * sys, const char * approx, const char * extension, const int id,
const liquid_params lp, char ** name ){
  	char id_char[3];
  	char phi_char[7];
	char dum_char1[7];
  	char dum_char2[7];
  	int int_phi_2;
  	int int_dum_1, int_dum_2, i1;
	double dum1;
	int_phi_2 =  lp.phi;
  	dum1 = ( lp.phi-int_phi_2 ) * 1e6;
  	int_phi_2 =  round(dum1);
  	sprintf( id_char, "%02d", id );
  	sprintf( phi_char, "%06d", int_phi_2 );
  	strcpy( * name, sys );
  	strcat( * name, "_" );
  	strcat( * name, approx );
  	strcat( * name, "_phi_0_" );
  	strcat( * name, phi_char );
	if ( lp.Tem != 0.0 ){
		int_dum_1 = (int) lp.Tem ;
		dum1 =  ( lp.Tem-int_dum_1 ) * 1e6;
		int_dum_2 = round(dum1);
		strcat( * name, "_T_" );
		sprintf( dum_char1, "%06d", int_dum_1 );
		strcat( * name, dum_char1 );
		strcat( * name, "_" );
		sprintf( dum_char2, "%06d", int_dum_2 );
		strcat( * name, dum_char2 );
		for ( i1=0; i1 < lp.nup; ++i1 ){
			int_dum_1 = lp.up[i1] ;
			dum1 = ( lp.up[i1]-int_dum_1 ) * 1e6;
			int_dum_2 = round(dum1);
			strcat( * name, "_P_" );
			sprintf( dum_char1, "%02d", int_dum_1 );
			strcat( * name, dum_char1 );
			strcat( * name, "_" );
			sprintf( dum_char2, "%02d", int_dum_2 );
			strcat( * name, dum_char2 );
			/*printf("%s\n", *name);*/
		}
	}
  strcat( * name, "_id_" );
  sprintf( id_char, "%02d", id );
  strcat( * name, id_char);
  strcat( * name, ".");
  strcat( * name, extension);
  return;
}

double
s_function_selector_mono_sph( const char * sys, const char * approx, const char * fun, const
double k, const liquid_params lp ){
  	if (strcmp( sys, "HS" ) == 0 ){
		if (strcmp( approx, "PY" ) == 0 ) {
			if (strcmp( fun, "ck" ) == 0 ){ return ck_hs_py( lp.phi, k );}
			else if (strcmp( fun, "is" ) == 0 ){ return is_hs_py( lp.phi, k );}
			else if (strcmp( fun, "sk" ) == 0 ){ return sk_hs_py( lp.phi, k );}
		}
		else if (strcmp( approx, "VW" ) == 0 ) {
			if (strcmp( fun, "ck" ) == 0 ){ return ck_hs_vw( lp.phi, k );}
			else if (strcmp( fun, "is" ) == 0 ){ return is_hs_vw( lp.phi, k );}
			else if (strcmp( fun, "sk" ) == 0 ){ return sk_hs_vw( lp.phi, k );}
		}
		else{ printf("Not a valid structural approximation for ");printf(sys); exit(1);}
	}
	else if ( strcmp( sys, "HSSW" ) == 0 ){
		double lambda = lp.up[0];
		if (strcmp( approx, "SH" ) == 0 ) {
			if (strcmp( fun, "ck" ) == 0 ){ return ck_hssw_vwsh( lp.phi, lp.Tem, lambda, k );}
			else if (strcmp( fun, "is" ) == 0 ){ return is_hssw_vwsh( lp.phi, lp.Tem, lambda, k );}
			else if (strcmp( fun, "sk" ) == 0 ){ return sk_hssw_vwsh( lp.phi, lp.Tem, lambda, k );}
		}
		else{ printf("Not a valid structural approximation for ");printf(sys); exit(1);}
	}
	else if ( strcmp( sys, "HSDBLEXP" ) == 0 ){
		if ( lp.nup == 3 ) {
			double Ta = lp.Tem;
			double Tr = lp.up[0];
			double za = lp.up[1];
			double zr = lp.up[2];
			if (strcmp( approx, "SH" ) == 0 ) {
				if (strcmp( fun, "ck" ) == 0 ){ return ck_dble_exp_vwsh( lp.phi, Ta, Tr, za, zr, k  );}
				else if (strcmp( fun, "is" ) == 0 ){ return is_dble_exp_vwsh( lp.phi, Ta, Tr, za, zr, k  );}
				else if (strcmp( fun, "sk" ) == 0 ){ return sk_dble_exp_vwsh( lp.phi, Ta, Tr, za, zr, k  );}
			}
			else{ printf("Not a valid structural approximation for ");printf(sys); exit(1);}
		}
		else{printf("Incorrect number of parameters for ");printf(sys); exit(1);}
		}
	else if ( strcmp( sys, "HSDBLEXP2" ) == 0 ){
		if ( lp.nup == 3 ) {
			double Ta = lp.Tem;
			double Tr = lp.up[0];
			double za = lp.up[1];
			double zr = lp.up[2];
			if (strcmp( approx, "SH" ) == 0 ) {
				if (strcmp( fun, "ck" ) == 0 ){ return ck_dble_exp_vwsh2( lp.phi, Ta, Tr, za, zr, k  );}
				else if (strcmp( fun, "is" ) == 0 ){ return is_dble_exp_vwsh2( lp.phi, Ta, Tr, za, zr, k  );}
				else if (strcmp( fun, "sk" ) == 0 ){ return sk_dble_exp_vwsh2( lp.phi, Ta, Tr, za, zr, k  );}
			}
			else{ printf("Not a valid structural approximation for ");printf(sys); exit(1);}
		}
		else{printf("Incorrect number of parameters for ");printf(sys); exit(1);}
	}
	else if ( strcmp( sys, "HSDBLEYUK" ) == 0 ){
		if ( lp.nup == 3 ) {
			double Ta = lp.Tem;
			double Tr = lp.up[0];
			double za = lp.up[1];
			double zr = lp.up[2];
			if (strcmp( approx, "SH" ) == 0 ) {
				if (strcmp( fun, "ck" ) == 0 ){ return ck_dble_yukawa_vwsh( lp.phi, Ta, Tr, za, zr, k  );}
				else if (strcmp( fun, "is" ) == 0 ){ return is_dble_yukawa_vwsh( lp.phi, Ta, Tr, za, zr, k  );}
				else if (strcmp( fun, "sk" ) == 0 ){ return sk_dble_yukawa_vwsh( lp.phi, Ta, Tr, za, zr, k  );}
			}
			else{ printf("Not a valid structural approximation for ");printf(sys); exit(1);}
		}
		else{printf("Incorrect number of parameters for ");printf(sys); exit(1);}
		}
	else if ( strcmp( sys, "HSDBLEYUK2" ) == 0 ){
		if ( lp.nup == 3 ) {
			double Ta = lp.Tem;
			double Tr = lp.up[0];
			double za = lp.up[1];
			double zr = lp.up[2];
			if (strcmp( approx, "SH" ) == 0 ) {
				if (strcmp( fun, "ck" ) == 0 ){ return ck_dble_yukawa_vwsh2( lp.phi, Ta, Tr, za, zr, k  );}
				else if (strcmp( fun, "is" ) == 0 ){ return is_dble_yukawa_vwsh2( lp.phi, Ta, Tr, za, zr, k  );}
				else if (strcmp( fun, "sk" ) == 0 ){ return sk_dble_yukawa_vwsh2( lp.phi, Ta, Tr, za, zr, k  );}
			}
			else{ printf("Not a valid structural approximation for ");printf(sys); exit(1);}
		}
		else{printf("Incorrect number of parameters for ");printf(sys); exit(1);}
	}
	else if ( strcmp( sys, "HD" ) == 0 ){
		if (strcmp( approx, "ROTH" ) == 0 ) {
			if (strcmp( fun, "ck" ) == 0 ){ return ck_hd_roth( lp.phi, k );}
			else if (strcmp( fun, "is" ) == 0 ){ return is_hd_roth( lp.phi, k );}
			else if (strcmp( fun, "sk" ) == 0 ){ return sk_hd_roth( lp.phi, k );}
		}
		else{ printf("Not a valid structural approximation for ");printf(sys); exit(1);}
	}
	else{printf("Not a valid structural system \n");exit(1);}
}

void
gsl_vector_s_function_selector_mono_sph(gsl_vector * sk, const char * sys, 
const char * approx, const char * fun, const gsl_vector * k, liquid_params lp ){
	if ( k!=NULL && sk!=NULL ){
		int i1; for (i1=0; i1<k->size; ++i1){ sk->data[i1]=s_function_selector_mono_sph( sys, approx, fun,  k->data[i1], lp ); }
		}
	else{ printf("Warning, check for unitialized structure or wave-vector vectors \n");}
	return;
}

double
radial_distribution_3D(double r, double rho, const structure_grid Sg){
	int i1;
	int size=Sg.k->size;
	double k; double kw; double S; double Lorch; double delta;
	double g;
	gsl_vector * integrand=gsl_vector_alloc(size);
	delta=M_PI / Sg.k->data[size-1];
	for (i1=0; i1< Sg.k->size; ++i1){
		k=Sg.k->data[i1];
		kw=Sg.kw->data[i1];
		S=Sg.S->data[i1];
		Lorch =  sin(delta*k); /* Correction function for truncated FT which is also divided by k, eliminated in the integrand */
		integrand->data[i1] =  ( S - 1.0 ) * sin( k * r ) * kw * Lorch ;
	}
	g = 1.0 + ( gsl_vector_sum(integrand) / ( 2.0 * M_PI * M_PI * rho * r * delta ) );
	
	gsl_vector_free(integrand);
	return g;
}

void
gsl_vector_radial_distribution_3D(gsl_vector * g, gsl_vector * r, double rho, const structure_grid Sg){
	if ( g !=NULL && r!=NULL ){
		int i1; for (i1=0; i1<r->size; ++i1){ g->data[i1]=radial_distribution_3D(r->data[i1], rho, Sg); }
		}
	else{ printf("Warning, check for unitialized structure or wave-vector vectors \n");}
	return;
}

void
s_grid_save_file( const structure_grid sg, const char * folder, const char * prefix, const char * suffix ){
	int i1,ts;
	ts=400*sizeof(char);
	char * complete_file = (char *)malloc(ts);
	strcpy(complete_file,folder);
	strcat(complete_file,prefix);
	strcat(complete_file,suffix);
	FILE * s_File = fopen(complete_file, "w");
	fprintf(s_File, "# || 1 k || 2 f(k) || 3 kw || # \n" );
	for (i1=0; i1<sg.k->size; ++i1){ 
		fprintf(s_File, "%1.9e \t %1.9e \t %1.9e \n", 
		sg.k->data[i1], sg.S->data[i1], sg.kw->data[i1] );
		}
	fclose(s_File);
	free(complete_file);
	return;
}

