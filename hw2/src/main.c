/* UCSB CS240A, Winter Quarter 2014
 * Main and supporting functions for the Conjugate Gradient Solver on a 5-point stencil
 *
 * NAMES: Tanmoy Sanyal
 * PERMS: 7550049
 * DATE: 
 */
#include "mpi.h"
#include "hw2harness.h"
#include <stdio.h>
#include <stdlib.h>

double* load_vec( char*, int* );
void save_vec( int, double*, int, double, double );

#define ROOT 0
double* dist_vec( int, int, double* );
void assemb_vec( int, double*, double* );
double pddot( int, int, double*, double* );
void pdaxpy( int, double, double, double*, double* );
void pmatvec(double* );


int main( int argc, char* argv[] ) {
	int writeOutX = 0;
	int n, k;
	int maxiterations = 3;
	int niters=0;
 	double norm;
	double *b, *loc_b; // local b vector
	double *x, *loc_x; //local sol vector
	double time;
	double t1, t2;
	
	int comm_rank, p, BLOCK_SIZE;
	int i;

	MPI_Init( &argc, &argv );
 	MPI_Comm_size(MPI_COMM_WORLD, &p);
 	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

	/*************************************** I/0 *************************************************/
	// 1st case runs model problem, 2nd Case allows you to specify your own b vector
	if ( argc == 3 ) {
		k = atoi( argv[1] );
		n = k*k;
		//build global b on the root proc from cs240_getB()
		if (comm_rank == ROOT) {
			b = (double *)malloc(n * sizeof(double) );
			for (i = 0; i < n; i++) b[i] = 1.0;
		}

	} else if  ( !strcmp( argv[1], "-i" ) && argc == 4 ) {
		// load b from file on the root proc
		if (comm_rank == ROOT) b = load_vec( argv[2], &k );
	
	} else {
		if (comm_rank == ROOT) {
		printf( "\nCGSOLVE Usage: \n\t"
			"Model Problem:\tmpirun -np [number_procs] cgsolve [k] [output_1=y_0=n]\n\t"
			"Custom Input:\tmpirun -np [number_procs] cgsolve -i [input_filename] [output_1=y_0=n]\n\n");
		}
	}
	writeOutX = atoi( argv[argc-1] ); // Write X to file if true, do not write if unspecified

	

	/************************************ INITIALIZE LOCAL VARIABLES ******************************/
	BLOCK_SIZE = (int) (n/p);

	// local b vector
	loc_b = dist_vec(comm_rank, BLOCK_SIZE, b);
	
	// local solution vector
	if (comm_rank == ROOT) {
		x = (double *)malloc(n *sizeof(double) );
		for (i = 0; i < n; i++) x[i] = 0.0;
	}
	loc_x = dist_vec(comm_rank, BLOCK_SIZE, x);

	// local residue and direction vector
	double* r = (double *)malloc(BLOCK_SIZE *sizeof(double) );
	double* d = (double *)malloc(BLOCK_SIZE *sizeof(double) );
	for (i = 0; i < BLOCK_SIZE; i++){
		r[i] = loc_b[i];
		d[i] = loc_b[i];
	}

	// iteration variables
    double alpha, beta, r_sq , rnew_sq;
    int isConverged = 0;
    double FRACTOL = 1e-3;
    double FRACTOL_sq = FRACTOL * FRACTOL; 
	


    /**************************************SOLVER ITERATIONS *************************************/ 
	// Start Timer
	t1 = MPI_Wtime();

    	while (niters < maxiterations && !isConverged){
    		/*
    		// compute new step size
    		Ad = pmatvec(d); r_sq = pddot(r, r); alpha = r_sq / pddot(d, Ad);

    		// update local copy of solution vector
    		pdaxpy(1.0, alpha, loc_x, d);
    	
    		// update local copy of residual vector
    		pdaxpy(1.0, -alpha, r_new, Ad);
    	
    		// compute new local copy of search direction vector
    		beta = pddot(r_new, r_new) / r_sq;
    		for (i = 0; i < sizeof(r)/sizeof(r[0]); i++) r[i] = r_new[i];
    		pdaxpy(1.0, beta, r, d);

    		// check norm
    		norm = r_sq / pddot(loc_b, loc_b);
    		if (norm < FRACTOL_sq) isConverged = 1;
    		*/

    		// update iterations
    		niters++;

    		// tests
    		
 		}
 	
 	// End Timer
	t2 = MPI_Wtime();
	
	
	/*************************************** POST PROCESSING ***************************************/
	// assemble local solution vectors
	assemb_vec(BLOCK_SIZE, loc_x, x);

	// Output
	if (comm_rank == ROOT){
		printf( "\nProblem size (k): %d\n",k);
		if(niters>0) printf( "\nNorm of the residual after %d iterations: %lf\n",niters,norm);
		printf( "\nElapsed time during CGSOLVE: %lf\n", t2-t1);
		if ( writeOutX) save_vec( k, x, niters, norm, t2-t1 );

	}
	
	/*Deallocate 
	if(niters > 0){
		free(loc_b); free(b);
		free(loc_x); free(x);
	}*/
	
	MPI_Finalize();
	
	return 0;
}



/**************************************** SUPPORTING FUNCTIONS *************************************/
 // distribution scheme returning local vectors for each processor
 double* dist_vec(int comm_rank, int BLOCK_SIZE, double* glob_vec){
 	double* loc_vec  = (double *)malloc ( BLOCK_SIZE * sizeof(double) );
	MPI_Scatter(glob_vec, BLOCK_SIZE, MPI_DOUBLE,
				loc_vec, BLOCK_SIZE, MPI_DOUBLE, 
				ROOT, MPI_COMM_WORLD);

	return loc_vec;
 }


 // assemble local vectors from each processor
 void assemb_vec(int BLOCK_SIZE, double* loc_vec, double* glob_vec){
 	MPI_Gather(loc_vec, BLOCK_SIZE, MPI_DOUBLE,
			   glob_vec, BLOCK_SIZE, MPI_DOUBLE, 
			   ROOT, MPI_COMM_WORLD);
 }


// parallel ddot
double pddot(int comm_rank, int BLOCK_SIZE, double* v, double* w){
	int i;
	double loc_dotp = 0.0;
	double glob_dotp;
	// compute local dot product
	for (i = 0; i < BLOCK_SIZE; i++)
		loc_dotp += v[i] * w[i];
	// sum the local dot products and send to all procs
	MPI_Allreduce(&loc_dotp, &glob_dotp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return glob_dotp;
}


// parallel daxpy
void pdaxpy(int BLOCK_SIZE, double alpha, double beta, double* v, double* w){
	int i;
	// compute local daxpy
	for (i = 0; i < BLOCK_SIZE; i++) 
		v[i] = alpha * v[i] + beta * w[i];
}


// parallel matvec
void pmatvec(double* w){
	printf("pass");
}


// Load Function
// NOTE: does not distribute data across processors
double* load_vec( char* filename, int* k ) {
	FILE* iFile = fopen(filename, "r");
	int nScan;
	int nTotal = 0;
	int n;
	
	if ( iFile == NULL ) {
		printf("Error reading file.\n");
		exit(0);
	}
	
	nScan = fscanf( iFile, "k=%d\n", k );
	if ( nScan != 1 ) {
		printf("Error reading dimensions.\n");
		exit(0);
	}
	
	n = (*k)*(*k);
	double* vec = (double *)malloc( n * sizeof(double) );
	
	do {
		nScan = fscanf( iFile, "%lf", &vec[nTotal++] );
	} while ( nScan >= 0 );
	
	if ( nTotal != n+1 ) {
		printf("Incorrect number of values scanned n=%d, nTotal=%d.\n",n,nTotal);
		exit(0);
	}
	
	return vec;
}


// Save a vector to a file (include data about k, niters, norm and time_elapsed)
void save_vec( int k, double* x, int niters, double norm, double delta_t ) { 
	FILE* oFile;
	int i;
	oFile = fopen("xApprox.txt","w");
	
	fprintf( oFile, "#k\tniters\tnorm\tdelta_t\n" );
	fprintf( oFile, "%d\t%d\t%lf\t%lf\n", k, niters, norm, delta_t );
	
	for (i = 0; i < k*k; i++) { 
    	fprintf( oFile, "%lf\n", x[i]);
 	} 

	fclose( oFile );
}
