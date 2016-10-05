#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>


void main( int argc, char* argv[] ) {
	int comm_rank, p, BLOCK_SIZE;
	int i , n;
	MPI_Status status;

	MPI_Init( &argc, &argv );
 	
 	MPI_Comm_size(MPI_COMM_WORLD, &p);
 	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
 	n = 10;
 	
 	if (comm_rank == 0){
 	        double* sendbuf = (double *)malloc(n * sizeof(double));
 	        for (i = 0; i < n; i++) sendbuf[i] = i;
 	        printf("\nProc id: %d\n", comm_rank);
 	        for (i = 0; i < n; i++) printf("\n%lf\n", sendbuf[i]);
 	        MPI_Send(sendbuf, n, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
 	}
 	
 	if (comm_rank == 1){
 	        double* recvbuf = (double *)malloc(n * sizeof(double));
         	MPI_Recv(recvbuf, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
 	        printf("\nProc id: %d\n", comm_rank);
 	        for (i = 0; i < n; i++) printf("\n%lf\n", recvbuf[i]);
 	        
 	}
 	MPI_Finalize();
}
 	
