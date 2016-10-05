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
void save_vec( int, int, int, int, double, double, double* );

#define ROOT 0
double* gen_B( int, int, int, int* );
double ddot( int, double*, double* );
void daxpy( int, double, double, double*, double* );
double* matvec( int, int, int, int, double* );


int main( int argc, char* argv[] ) {
	int writeOutX = 0;
	int n, k;
	int maxiterations = 3;
	int niters=0;
 	double norm;
	double* b; 
	double* x;
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
		//generate vector b on each proc
		b = gen_B(comm_rank, n, p, &BLOCK_SIZE);

	} else if  ( !strcmp( argv[1], "-i" ) && argc == 4 ) {
		// load b from file
		b = load_vec( argv[2], &k );
	
	} else {
		if (comm_rank == ROOT) {
		printf( "\nCGSOLVE Usage: \n\t"
			"Model Problem:\tmpirun -np [number_procs] cgsolve [k] [output_1=y_0=n]\n\t"
			"Custom Input:\tmpirun -np [number_procs] cgsolve -i [input_filename] [output_1=y_0=n]\n\n");
		}
	}
	writeOutX = atoi( argv[argc-1] ); // Write X to file if true, do not write if unspecified

	
	/************************************ INITIALIZE VARIABLES ON EACH PROC******************************/
	x = (double *)malloc(BLOCK_SIZE * sizeof(double) );
	double* r = (double *)malloc(BLOCK_SIZE * sizeof(double) );
	double* r_new = (double *)malloc(BLOCK_SIZE * sizeof(double) );
	double* Ad = (double *)malloc(BLOCK_SIZE * sizeof(double) );
	double* d = (double *)malloc(BLOCK_SIZE * sizeof(double) );
	norm = 0.0;
	for (i = 0; i < BLOCK_SIZE; i++){
		x[i] = 0.0;
		r[i] = b[i];
		r_new[i] = b[i];
		d[i] = b[i];
	}
	
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
 		}
 	
 	// End Timer
	t2 = MPI_Wtime();
	double* testprod;
	testprod = matvec(comm_rank, p, k, BLOCK_SIZE, b);
	daxpy(BLOCK_SIZE, 1.0, 1.0, x, testprod);

	/*************************************** POST PROCESSING ***************************************/
	// synchronized write to file 
	if ( writeOutX){
		int ii;
		for (ii = 0; ii < p; ++ii){
			if (comm_rank == ii) save_vec(comm_rank, BLOCK_SIZE, k, niters, t2-t1, norm, x );
			MPI_Barrier(MPI_COMM_WORLD); // prevent other procs from writing before this proc has finished
		}
	}
	
	// print output on stdout from proc 0
	if (comm_rank == ROOT){
		printf( "\nProblem size (k): %d\n",k);
		if(niters>0) printf( "\nNorm of the residual after %d iterations: %lf\n",niters,norm);
		printf( "\nElapsed time during CGSOLVE: %lf\n", t2-t1);
	}
	
	//Deallocate 
	if(niters > 0){
		free(b);
		free(x);
		free(r);
		free(r_new);
		free(d);
		free(Ad);
	}
	
	MPI_Finalize();
	
	return 0;
}



/**************************************** SUPPORTING FUNCTIONS *************************************/
// generate parallel copy of vector b
double* gen_B(int comm_rank, int n, int p, int* BLOCK_SIZE){
	int i;
	*BLOCK_SIZE = (int)(n/p);
	double* vec = (double*)malloc(*BLOCK_SIZE * sizeof(double));
	for (i = 0; i < *BLOCK_SIZE; i++)
		//vec[i] = cs240_getB( i+comm_rank*(*BLOCK_SIZE), n);
		vec[i] = i + comm_rank*(*BLOCK_SIZE);
	return vec;
}


// parallel ddot
double ddot(int BLOCK_SIZE, double* v, double* w){
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
void daxpy(int BLOCK_SIZE, double alpha, double beta, double* v, double* w){
	int i;
	// compute local daxpy
	for (i = 0; i < BLOCK_SIZE; i++) 
		v[i] = alpha * v[i] + beta * w[i];
}


// parallel matvec
double* matvec(int comm_rank, int p, int k, int BLOCK_SIZE, double* w){
	double* v = (double *)malloc(BLOCK_SIZE * sizeof (double) );
	int i, mpierr;
	int sendtag = 0; int recvtag = 0;
	MPI_Status status;

	/*w_top refers to top row of current proc inherited from neigh proc
	w_bottom refers to bottom row of current proc inherited from neigh proc*/
	
	if (comm_rank == ROOT){
		double* sendbuf = (double *)malloc(k * sizeof(double) );
		double* w_bottom = (double *)malloc(k * sizeof(double) );

		// populate send buffers
		for (i = 0; i < k; i++) sendbuf[i] = w[BLOCK_SIZE - k + i];	

		// send bottom row to proc 1 receive bottom row from proc 1
		MPI_Send(sendbuf, k, MPI_DOUBLE, 1, sendtag, MPI_COMM_WORLD);
		MPI_Recv(w_bottom, k, MPI_DOUBLE, 1, recvtag, MPI_COMM_WORLD, &status);

		// compute matvec using 5 point stencil
		for (i = 0; i < BLOCK_SIZE; i++) {
			v[i] = 4*w[i];
			if (i % k == 0) v[i] -= w[i+1]; // left boundary
			if ( (i+1) % k == 0) v[i] -= w[i-1]; // right boundary
			if (i < k) v[i] -= w[i+k]; 
			if (i > k && i < BLOCK_SIZE-k) v[i] -= (w[i-k] + w[i+k]); // middle rows
			if (i >= BLOCK_SIZE-k && i < BLOCK_SIZE) v[i] -= w_bottom[i - BLOCK_SIZE + k]; // bottom row
		}


	} else if (comm_rank == p-1) {
		double* w_top = (double *)malloc(k * sizeof(double) ); 
		double* sendbuf = (double *)malloc(k * sizeof(double) ); 

		// populate send buffers
		for (i = 0; i < k; i++) sendbuf[i] = w[i];	

		// send top row to proc p-2 and receive own top row from proc p-2
		MPI_Recv(w_top, k, MPI_DOUBLE, p-2, recvtag, MPI_COMM_WORLD, &status);
		MPI_Send(sendbuf, k, MPI_DOUBLE, p-2, sendtag, MPI_COMM_WORLD); 

		// compute matvec using 5 point stencil
		for (i = 0; i < BLOCK_SIZE; i++) {
			v[i] = 4*w[i];
			if (i % k == 0) v[i] -= w[i+1]; // left boundary
			if ( (i+1) % k == 0) v[i] -= w[i-1]; // right boundary 
			if (i < k) v[i] -= w_top[i]; // top row
			if (i > k && i < BLOCK_SIZE-k) v[i] -= (w[i-k] + w[i+k]); // middle rows
		}


	} else {
		double* sendbuf_top = (double *)malloc(k * sizeof(double) );
		double* sendbuf_bottom = (double *)malloc(k * sizeof(double) );
		double* w_top = (double *)malloc(k * sizeof(double) );
		double* w_bottom = (double *)malloc(k * sizeof(double) );

		// populate send buffers
		for (i = 0; i < k; i++){
			sendbuf_top[i] = w[i];
			sendbuf_bottom[i] = w[BLOCK_SIZE - k + i];
		}		

		// send top row to prev proc and receive own top row from prev proc
		MPI_Recv(w_top, k, MPI_DOUBLE, comm_rank-1, recvtag, MPI_COMM_WORLD, &status);
		MPI_Send(sendbuf_top, k, MPI_DOUBLE, comm_rank-1, sendtag, MPI_COMM_WORLD); 

		// send bottom row to next proc and receive own bottom row from next proc
		MPI_Send(sendbuf_bottom, k, MPI_DOUBLE, comm_rank+1, sendtag, MPI_COMM_WORLD); 
		MPI_Recv(w_bottom, k, MPI_DOUBLE, comm_rank+1, recvtag, MPI_COMM_WORLD, &status);

		// compute matvec using 5 point stencil
		for (i = 0; i < BLOCK_SIZE; i++) {
			v[i] = 4*w[i];
			if (i % k == 0) v[i] -= w[i+1]; // left boundary
			if ( (i+1) % k == 0) v[i] -= w[i-1]; // right boundary
			if (i < k) v[i] -= w_top[i]; // top row
			if (i > k && i < BLOCK_SIZE-k) v[i] -= (w[i-k] + w[i+k]); // middle rows
			if (i >= BLOCK_SIZE-k && i < BLOCK_SIZE) v[i] -= w_bottom[i - BLOCK_SIZE + k]; // bottom row
		}
	}

	return v;
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
void save_vec( int comm_rank, int BLOCK_SIZE, int k, int niters, double delta_t, double norm, double* x ) { 
	FILE* oFile;
	int i;

	oFile = fopen("xApprox.txt","a");
	
	if (comm_rank == ROOT){
		fprintf( oFile, "#k\tniters\tnorm\tdelta_t\n" );
		fprintf( oFile, "%d\t%d\t%lf\t%lf\n", k, niters, norm, delta_t );
	}
	
	for (i = 0; i < BLOCK_SIZE; i++) { 
    	fprintf( oFile, "%lf\n", x[i]);
 	} 

	fclose( oFile );
}
