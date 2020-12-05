#include<stdio.h>
#include <mpi.h>

#define N 700

int n;

int main(int argc, char *argv[]) {
	int i, j, k, id, p, pid, size, st, cell, left, right, rid;

	double time;

	MPI_Status status;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &p);
	
	right = id;
	left = (p+id-1)%p;
	
	n = N;
	if (N%p != 0) {
		n += (p-N%p);		
	}

	int a[n][n];
	int b[n][n];
	int c[n][n];
	int bt[n][n];


	if (id == 0) {
		// populate 1-D array
		// ------------------------------------------------
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				if (i >= N || j >= N) {
					// padding
					a[i][j] = 0;
					b[i][j] = 0;					
				} else {
					// assign value
					a[i][j] = 2;
					b[i][j] = 3;
				}	
			}
		}
		// ------------------------------------------------

		// sequential matrix multiplication
		// ------------------------------------------------
		time = MPI_Wtime();

		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				cell = 0;
				for (k = 0; k < n; ++k) {
					cell += (a[j][k] * b[k][i]);				
				}
				c[j][i] = cell;
			}	
		}
		time  = MPI_Wtime() - time;
		printf("t_s: %.4f\n", time);
		// printf("%d 1 %.4f, ", n, time);
		// ------------------------------------------------

		// transpose b 
		// ------------------------------------------------
		for (i = 0; i < n; i++) {
	        for (j = 0; j < n; j++) {
	            bt[i][j] = b[j][i];
			}
        }
		// ------------------------------------------------

		// ------------------------------------------------
		time = MPI_Wtime();
		size = (n)/p;

		// broadcast populated array to other nodes
		for (pid=1; pid < p; pid++) {
    		MPI_Send(&a[pid*size][0], size*n, MPI_INT,pid,1, MPI_COMM_WORLD);
      		MPI_Send(&bt[pid*size][0], size*n, MPI_INT, pid, 1, MPI_COMM_WORLD); 
	    }

	    // ring communication to share complete bt among all processors
	    for (rid = 1; rid < p; ++rid) {  	  			
  			MPI_Send(&bt[right*size][0], size*n, MPI_INT, 1, 3, MPI_COMM_WORLD); 				
			right = (n+right-1)%p;
		}


	    // sequential matrix multiplication
		for (i = 0; i < n; ++i) {
			for (j = 0; j < size; ++j) {
				cell = 0;
				for (k = 0; k < n; ++k) {
					cell += (a[j][k] * bt[j][k]);				
				}
				c[j][i] = cell;
			}	
		}

		// wait until the results from all nodes
		for (pid=1; pid < p; pid++) {
    		MPI_Recv(&st, 1, MPI_INT, pid, 2, MPI_COMM_WORLD, &status);
			MPI_Recv(&size, 1, MPI_INT, pid, 2, MPI_COMM_WORLD, &status);
			MPI_Recv(&c[st][0], size*n, MPI_INT, pid, 2, MPI_COMM_WORLD, &status); 
	    }

	    time = MPI_Wtime() - time;
	    printf("t_p: %.4f\n", time);
	    // printf("%d %d %.4f\n", n, p, time);

		// printf("Result \n");				
		// // print result
		// for (i = 0; i < N; ++i)
		// {
		// 	for (j = 0; j < N; ++j)
		// 	{
		// 		printf("%d ", c[i][j]);
		// 	}
		// 	printf("\n");	
		// }
		// ------------------------------------------------

	} else {
		// ------------------------------------------------
		size = n/p;		
		st = id*size;
	    MPI_Recv(&a, size*n, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
	    MPI_Recv(&bt[st][0], size*n, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);

	    // ring communication to share complete bt among all processors
	    // continue with matrix multiplication after comlete bt is reveived
  		for (rid = 1; rid < p; ++rid) {  
  			if (id == p-1) {
				MPI_Recv(&bt[left*size][0], size*n, MPI_INT, id-1, 3, MPI_COMM_WORLD, &status);
			}else if (id%2 == 0) {
	  			MPI_Send(&bt[right*size][0], size*n, MPI_INT, id+1, 3, MPI_COMM_WORLD); 
		    	MPI_Recv(&bt[left*size][0], size*n, MPI_INT, id-1, 3, MPI_COMM_WORLD, &status);				
			} else {
		    	MPI_Recv(&bt[left*size][0], size*n, MPI_INT, id-1, 3, MPI_COMM_WORLD, &status);
	  			MPI_Send(&bt[right*size][0], size*n, MPI_INT, id+1, 3, MPI_COMM_WORLD); 		    	
			}
			right = (n+right-1)%p;
			left = (n+left-1)%p;
		}

		// if (id == 5)
		// {			
		// 	printf("collected 1\n");				
		// 	// print result
		// 	for (i = 0; i < n; ++i)
		// 	{
		// 		for (j = 0; j < n; ++j)
		// 		{
		// 			printf("- %d ", bt[i][j]);
		// 		}
		// 		printf("\n");	
		// 	}
		// }		
	    
		// sequential matrix multiplication
		for (i = 0; i < n; ++i) {
			for (j = 0; j < size; ++j) {
				cell = 0;
				for (k = 0; k < n; ++k) {
					cell += (a[j][k] * bt[j][k]);				
				}
				c[j][i] = cell;
			}	
		}

		// send results to root
		MPI_Send(&st, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	    MPI_Send(&size, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	    MPI_Send(&c, size*n, MPI_INT, 0, 2, MPI_COMM_WORLD);
		// ------------------------------------------------

	}

	// cleanup 
	MPI_Finalize();

	return 0;
}