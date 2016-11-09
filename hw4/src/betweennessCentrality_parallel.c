#include "defs.h"
#include <stdio.h>
#include <omp.h>

#define MAX_THREADS 48

/*
 * Parallel Version
 *
 */
double betweennessCentrality_parallel(graph* G, double* BC) {
  int *S; 	/* stack of vertices in order of distance from s. Also, implicitly, the BFS queue */
  plist* P;  	/* predecessors of vertex v on shortest paths from s */
  double* sig; 	/* No. of shortest paths */
  int* d; 	/* Length of the shortest path between every pair */
  double* del; 	/* dependency of vertices */
  int *in_degree, *numEdges;
  int *pListMem;	
  int* Srcs; 
  int *start, *end;
  int seed = 2387;
  double elapsed_time;
  int i, j, k, p, count, myCount;
  int v, w, vert;
  int numV, num_traversals, n, m, phase_num;

  /* numV: no. of vertices to run BFS from = 2^K4approx */
  //numV = 1<<K4approx;
  n = G->nv;
  m = G->ne;
  numV = 100;

  /* Permute vertices */
  Srcs = (int *) malloc(n*sizeof(int));
  for (i=0; i<n; i++) {
    Srcs[i] = i;
  }

  /* Start timing code from here */
  elapsed_time = get_seconds();

  /* Initialize predecessor lists */
  /* Number of predecessors of a vertex is at most its in-degree. */
  P = (plist *) calloc(n*MAX_THREADS, sizeof(plist));
  in_degree = (int *) calloc(n+1, sizeof(int));
  numEdges = (int *) malloc((n+1)*sizeof(int));
  for (i=0; i<m; i++) {
    v = G->nbr[i];
    in_degree[v]++;
  }
  prefix_sums(in_degree, numEdges, n);
  pListMem = (int *) malloc(m*sizeof(int));
  for (i=0; i<MAX_THREADS; i++) {
  	plist* thisP = &P[i*n];
  	for (j=0; j<n; j++){
    	thisP[j].list = pListMem + numEdges[j];
    	thisP[j].degree = in_degree[j];
    	thisP[j].count = 0;
  	}
  }
  free(in_degree);
  free(numEdges);
	
  /* Allocate shared memory */ 
  S   = (int *) malloc(n*MAX_THREADS * sizeof(int));
  sig = (double *) malloc(n*MAX_THREADS * sizeof(double));
  d   = (int *) malloc(n*MAX_THREADS * sizeof(int));
  del = (double *) calloc(n*MAX_THREADS, sizeof(double));
	
  start = (int *) malloc(n*MAX_THREADS * sizeof(int));
  end = (int *) malloc(n*MAX_THREADS * sizeof(int));

  num_traversals = 0;

  for (i=0; i<n*MAX_THREADS; i++) {
    d[i] = -1;
  }

  /* keep track of how many threads were actually used */
  int nthreads; // note: this is shared because it'll be needed to update BC after the main loop at the end

  /* reducer array to hold BC updates */
  double* BC_reducer = (double *) calloc(n*MAX_THREADS, sizeof(double));

  /* zero out given BC array */
  for(i=0; i<n; i++) BC[i] = 0.0;

  /***********************************/
  /*** MAIN LOOP *********************/
  /***********************************/
  #pragma omp parallel for
  for (p=0; p<n ; p++) {

  		// get thread id
  		int tid = omp_get_thread_num();

  		// find actual number of threads used
  		nthreads = omp_get_num_threads();

  		/* declare thread private copies of stuff */
  		int *myS; 	/* stack of vertices in order of distance from s. Also, implicitly, the BFS queue */
  		plist* myP;  	/* predecessors of vertex v on shortest paths from s */
  		double* mysig; 	/* No. of shortest paths */
  		int* myd; 	/* Length of the shortest path between every pair */
  		double* mydel; 	/* dependency of vertices */
  		int *mystart, *myend;
  		int i, j, k, count, phase_num, myCount;
  		int v, w, vert;
  		double* myBC;

		i = Srcs[p];
		if (G->firstnbr[i+1] - G->firstnbr[i] == 0) continue;
		else {
				if (num_traversals == numV + 1) continue;
				else {
					#pragma omp atomic
						num_traversals++;
				}
		}

		/* point thread private copy to appropriate location on corresponding shared array */
		myS = &S[tid*n];
		myP = &P[tid*n];
		mysig = &sig[tid*n];
		myd = &d[tid*n];
		mydel = &del[tid*n];
		mystart = &start[tid*n];
		myend = &end[tid*n];
		myBC = &BC_reducer[tid*n];
		
		mysig[i] = 1;
		myd[i] = 0;
		myS[0] = i;
		mystart[0] = 0;
		myend[0] = 1;
		
		count = 1;
		phase_num = 0;

		while (myend[phase_num] - mystart[phase_num] > 0) {
				myCount = 0;
				// BFS to destination, calculate distances
				for ( vert = mystart[phase_num]; vert < myend[phase_num]; vert++ ) {
					v = myS[vert];
					for ( j=G->firstnbr[v]; j<G->firstnbr[v+1]; j++ ) {
						w = G->nbr[j];
						if (v != w) {

							//w found for the first time? 
							if (myd[w] == -1) {
								myS[myend[phase_num] + myCount] = w;
								myCount++;
								myd[w] = myd[v] + 1; 
								mysig[w] = mysig[v]; 
								myP[w].list[myP[w].count++] = v;
							} else if (myd[w] == myd[v] + 1) {
								mysig[w] += mysig[v]; 
								myP[w].list[myP[w].count++] = v;
							}
						
						}
					}
	 			}
			
				// Merge all local stacks for next iteration
				phase_num++ ; 
				mystart[phase_num] = myend[phase_num-1];
				myend[phase_num] = mystart[phase_num] + myCount;
				count = myend[phase_num];
			
		}

		while (phase_num > 0) {
			for (j=mystart[phase_num]; j<myend[phase_num]; j++) {
				w = myS[j];
				for (k = 0; k < myP[w].count; k++) {
					v = myP[w].list[k];
					mydel[v] = mydel[v] + mysig[v]*(1+mydel[w])/mysig[w];
				}
				
				//#pragma omp atomic
				myBC[w] += mydel[w];
			}

			phase_num--;
		}
		
		for (j=0; j<count; j++) {
			w = myS[j];
			myd[w] = -1;
			mydel[w] = 0;
			myP[w].count = 0;
		}
		
  }
  /***********************************/
  /*** END OF MAIN LOOP **************/
  /***********************************/

  /********* END PARALLEL SECTION **********/

  // accumulate thread local BC reductions into output array BC
  for(i=0; i<nthreads; i++){
  	for(j=0; j<n; j++){
  		BC[j] += BC_reducer[i*n+j];
  	}
  }

  free(BC_reducer);

  free(S);
  free(pListMem);
  free(P);
  free(sig);
  free(d);
  free(del);
  free(start);
  free(end);
  elapsed_time = get_seconds() - elapsed_time;
  free(Srcs);

  return elapsed_time;
}
