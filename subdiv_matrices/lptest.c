/*
	--- simple program that prints out the Loop eigenstructures that were
		 precomputed.

		 Author : Jos Stam (jstam@aw.sgi.com), Alias|wavefront
*/

#include <stdio.h>


#define IX(i,j,n) ((i)+(n)*(j))

typedef
   struct
   {
      double * val;
      double * vecI;
      double ** Phi;
   } EVALSTRUCT;

/* --- this routine reads the eigen-data from the file structure and returns a
       pointer to a data structure containing it */

static EVALSTRUCT ** read_eval ( int * pNmax )
{
	EVALSTRUCT ** ev;
	FILE * f;
	int Nmax, i, N, K;

#ifdef _WIN32
	if ( !(f=fopen("lpdata50NT.dat","rb")) ) return ( NULL );
#else
	if ( !(f=fopen("lpdata50.dat","r")) ) return ( NULL );
#endif

	fread ( &Nmax, sizeof(int), 1, f );

	ev = (EVALSTRUCT **) malloc ( (Nmax-2)*sizeof(EVALSTRUCT *) );

	for ( i=0 ; i<Nmax-2 ; i++ )
	{
		N = i+3;
		K = N+6;

		ev[i] = (EVALSTRUCT *) malloc ( sizeof(EVALSTRUCT) );
		ev[i]->val = (double *) malloc ( K*sizeof(double) );
		ev[i]->vecI = (double *) malloc ( K*K*sizeof(double) );
		ev[i]->Phi = (double **) malloc ( 3*sizeof(double) );
		ev[i]->Phi[0] = (double *) malloc ( K*12*sizeof(double) );
		ev[i]->Phi[1] = (double *) malloc ( K*12*sizeof(double) );
		ev[i]->Phi[2] = (double *) malloc ( K*12*sizeof(double) );

		fread ( ev[i]->val, sizeof(double), K, f );
		fread ( ev[i]->vecI, sizeof(double), K*K, f );
		fread ( ev[i]->Phi[0], sizeof(double), K*12, f );
		fread ( ev[i]->Phi[1], sizeof(double), K*12, f );
		fread ( ev[i]->Phi[2], sizeof(double), K*12, f );
	}

	fclose ( f );

	*pNmax = Nmax;

	return ( ev );
}

/* --- just print out the eigen data */

static void print_eval ( EVALSTRUCT ** ev, int Nmax )
{
	int i, j, k, l, N, K;

	printf ( "Nmax = %d\n\n", Nmax );

	for ( i=0 ; i<Nmax-2 ; i++ )
	{
		N = i+3;
		K = N+6;

		printf ( "N=%d\n\n", N );

		printf ( "Eigenvalues:\n\n" );

		for ( j=0 ; j<K ; j++ ) printf ( "%6.3f ", ev[i]->val[j] );
		printf ( "\n\n" );

		printf ( "Inverse of Eigenvectors:\n\n" );

		for ( j=0 ; j<K ; j++ )
		{
			for ( k=0 ; k<K ; k++ ) printf ( "%6.3f ", ev[i]->vecI[IX(j,k,K)] );
			printf ( "\n" );
		}
		printf ( "\n\n" );

		printf ( "Coefficients of the Eigenbasis functions\n\n" );

		for ( k=0 ; k<3 ; k++ )
		{
			printf ( "k=%d :\n\n", k );

			for ( j=0 ; j<K ; j++ )
			{
				for ( l=0 ; l<12 ; l++ ) printf("%6.3f ",ev[i]->Phi[k][IX(j,l,K)] );
				printf ( "\n" );
			}
			printf ( "\n\n" );
		}

		printf ( "\n\n" );

	}

	return;
}

main ( int argc, char ** argv )
{
	EVALSTRUCT ** ev;
	int Nmax;

	ev = read_eval ( &Nmax );
	if ( !ev ) exit ( 1 );

	print_eval ( ev, Nmax );

	exit ( 0 );
}
