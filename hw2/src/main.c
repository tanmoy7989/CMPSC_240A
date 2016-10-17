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
void save_vec( int, int, int, int, int, double, double, double* );

#define ROOT 0
double* gen_B( int, int, int, int* );
double ddot( int, double*, double* );
void daxpy( int, double, double, double*, double* );
double* matvec( int, int, int, int, double* );


int main( int argc, char* argv[] ) {
	int writeOutX = 0; int verifyX = 0;
	int n, k;
	int maxiterations = 1000;
	int niters=0;
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
	if ( argc == 4 ) {
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
			"Model Problem:\tmpirun -np [number_procs] cgsolve [k] [showoutput (0/1)] [verify (0/1)] \n\t"
			"Custom Input:\tmpirun -np [number_procs] cgsolve -i [input_filename] [showoutput (0/1)] [verify (0/1)] \n\n");
		}
		MPI_Finalize();
		return 1;
	}
	
	writeOutX = atoi( argv[argc-2] ); // Write X to file if true
	verifyX = atoi( argv[argc-1] ); // verify output if true (only use for small k values as this needs assembling a global vector)

	
	/************************************ INITIALIZE VARIABLES ON EACH PROC******************************/
	x = (double *)malloc(BLOCK_SIZE * sizeof(double) );
	double* r = (double *)malloc(BLOCK_SIZE * sizeof(double) );
	double* Ad = (double *)malloc(BLOCK_SIZE * sizeof(double) );
	double* d = (double *)malloc(BLOCK_SIZE * sizeof(double) );
	
	double alpha, beta, rtr , rtrold;
	double relnorm = 1.0; double normb = ddot(BLOCK_SIZE, b, b);
	
	for (i = 0; i < BLOCK_SIZE; i++){
		x[i] = 0.0;
		r[i] = b[i];
		d[i] = b[i];
	}
	rtr = ddot(BLOCK_SIZE, r, r);
	
    double FRACTOL = 1e-2;
    double FRACTOL_sq = FRACTOL * FRACTOL; 
	

    /**************************************SOLVER ITERATIONS *************************************/ 
	// Start Timer
	t1 = MPI_Wtime();

    	while (niters < maxiterations && relnorm > FRACTOL_sq){
    		// compute new step size
    		Ad = matvec(comm_rank, p, k, BLOCK_SIZE, d); alpha = rtr / ddot(BLOCK_SIZE, d, Ad);

    		// update local copy of solution vector
    		daxpy(BLOCK_SIZE, 1.0, alpha, x, d);
    	
    		// update local copy of residual vector
    		daxpy(BLOCK_SIZE, 1.0, -alpha, r, Ad);
    	
    		// compute new local copy of search direction vector
    		rtrold = rtr; rtr = ddot(BLOCK_SIZE, r, r);
    		beta =  rtrold/ rtr;
    		daxpy(BLOCK_SIZE, 1.0, beta, d, r);

    		// check norm
    		relnorm = rtr / normb;

    		// update iterations
    		niters++;
 		}
 	
 	// End Timer
	t2 = MPI_Wtime();

	/*************************************** POST PROCESSING ***************************************/
	// synchronized write to file w
	if ( writeOutX ){
		int ii;
		for (ii = 0; ii < p; ++ii){
			if (comm_rank == ii) save_vec(comm_rank, p, BLOCK_SIZE, k, niters, t2-t1, relnorm, x );
			MPI_Barrier(MPI_COMM_WORLD); // prevent other procs from writing before this proc has finished
		}
	}
	
	// print output on stdout from proc 0
	if (comm_rank == ROOT){
		printf( "\nProblem size (k): %d\n",k);
		if(niters>0) printf( "\nNorm of the residual after %d iterations: %lf\n",niters,relnorm);
		printf( "\nElapsed time during CGSOLVE: %lf\n", t2-t1);
	}

	// verify solution using test harness
	double* x_global;
	if ( verifyX ) {
		if (comm_rank == ROOT) x_global = (double *)malloc((BLOCK_SIZE*p) * sizeof(double));
		MPI_Gather(x, BLOCK_SIZE, MPI_DOUBLE, x_global, BLOCK_SIZE, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

		if (comm_rank == ROOT){
			int correct = cs240_verify(x_global, k, t2-t1);
			printf("Correct = %d\n", correct);
		}
	}

	//Deallocate 
	if(niters > 0){
		free(b); free(x); free(r); free(d); free(Ad);
		if (comm_rank == ROOT && verifyX) free(x_global);
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
		vec[i] = cs240_getB( i+comm_rank*(*BLOCK_SIZE), n);
		//vec[i] = i + 3 + comm_rank*(*BLOCK_SIZE);
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
	int i;
	int sendtag = 0; int recvtag = 0;
	MPI_Status status;

	// stencil variables
	double neigh_top, neigh_bottom, neigh_left, neigh_right;

	// single_proc
	if (p == 1){
		for (i = 0; i < BLOCK_SIZE; i++){
			if (i % k == 0) neigh_left = 0.0; else neigh_left = w[i-1];
			if ( (i+1) % k == 0) neigh_right = 0.0; else neigh_right = w[i+1];
			if (i < k) neigh_top = 0.0; else neigh_top = w[i-k];
			if (i >= BLOCK_SIZE-k) neigh_bottom = 0.0; else neigh_bottom = w[i+k];
			v[i] = 4*w[i] - neigh_left - neigh_right - neigh_top - neigh_bottom;
			}
		return v;
	}

	// 1 ghost layer (sort of)
	double *w_top, *w_bottom;

	/*w_top refers to top row of current proc inherited from neigh proc
	w_bottom refers to bottom row of current proc inherited from neigh proc*/

	if (comm_rank == ROOT){
		double* sendbuf = (double *)malloc(k * sizeof(double) );
		w_top = NULL;
		w_bottom = (double *)malloc(k * sizeof(double) );

		// populate send buffers
		for (i = 0; i < k; i++) sendbuf[i] = w[BLOCK_SIZE - k + i];	

		// send bottom row to proc 1 receive bottom row from proc 1
		MPI_Send(sendbuf, k, MPI_DOUBLE, 1, sendtag, MPI_COMM_WORLD);
		MPI_Recv(w_bottom, k, MPI_DOUBLE, 1, recvtag, MPI_COMM_WORLD, &status);


	} else if (comm_rank == p-1) {
		double* sendbuf = (double *)malloc(k * sizeof(double) ); 
		w_top = (double *)malloc(k * sizeof(double) ); 
		w_bottom = NULL;

		// populate send buffers
		for (i = 0; i < k; i++) sendbuf[i] = w[i];	

		// send top row to proc p-2 and receive own top row from proc p-2
		MPI_Recv(w_top, k, MPI_DOUBLE, p-2, recvtag, MPI_COMM_WORLD, &status);
		MPI_Send(sendbuf, k, MPI_DOUBLE, p-2, sendtag, MPI_COMM_WORLD); 


	} else {
		double* sendbuf_top = (double *)malloc(k * sizeof(double) );
		double* sendbuf_bottom = (double *)malloc(k * sizeof(double) );
		w_top = (double *)malloc(k * sizeof(double) );
		w_bottom = (double *)malloc(k * sizeof(double) );

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
	}

	// apply 5 point stencil
	for (i = 0; i < BLOCK_SIZE; i++){
		if (i % k == 0) neigh_left = 0.0; else neigh_left = w[i-1];

		if ( (i+1) % k == 0) neigh_right = 0.0; else neigh_right = w[i+1];

		if (i < k){
			if (comm_rank == ROOT) neigh_top = 0.0; else neigh_top = w_top[i];
			if (BLOCK_SIZE > k) neigh_bottom = w[i+k];
		}

		if (i > k && i < BLOCK_SIZE-k) {
			neigh_top = w[i-k];
			neigh_bottom = w[i+k];
		}
		
		if ( i >= BLOCK_SIZE-k && i < BLOCK_SIZE){
			if (BLOCK_SIZE > k) neigh_top = w[i-k];
			if (comm_rank == p-1) neigh_bottom = 0.0; else neigh_bottom = w_bottom[i-BLOCK_SIZE+k];
		}
		
		v[i] = 4*w[i] - neigh_left - neigh_right - neigh_top - neigh_bottom;
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
void save_vec( int comm_rank, int p, int BLOCK_SIZE, int k, int niters, double delta_t, double norm, double* x ) { 
	FILE* oFile;
	int i;

	oFile = fopen("xApprox.txt","a");
	
	if (comm_rank == ROOT){
		fprintf( oFile, "#k\tp\tniters\tnorm\tdelta_t\n" );
		fprintf( oFile, "%d\t%d\t%d\t%lf\t%lf\n", k, p, niters, norm, delta_t );
	}
	
	for (i = 0; i < BLOCK_SIZE; i++) { 
    	fprintf( oFile, "%lf\n", x[i]);
 	} 

	fclose( oFile );
}
