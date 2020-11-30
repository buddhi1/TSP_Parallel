#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>

using namespace std;

#define vr 10

int *permutationArray; 
long lastRow = 0;
long numPerm = 1; 

void printArray(int v[], int s) {
	for (int i = 0; i < s; ++i)
	{
		cout << v[i] << " ";
	}
	cout << endl;
}

// populate graph
void setGraph(int graph[vr][vr]) {
	int num;
	for (int i = 0; i < vr; ++i)
	{
		for (int j = i; j < vr; ++j)
		{
			if (i == j)
			{
				graph[i][j] = 0;
			} else {
				num = rand() % 10 + 1;
				graph[i][j] = num; 
				graph[j][i] = num; 
			}
		}
	}
}

void setPermutationRow(int a[],int n) {
	for (int j = 0; j < n; j++) {
    	*(permutationArray + lastRow*(vr-1) + j) = a[j]; 	
	}
	lastRow++;	
}

// Generating permutation using Heap Algorithm
void heapPermutation(int a[vr-1], int size, int n)
{
    // if size becomes 1 then prints the obtained
    // permutation
    if (size == 1) {
    	// printArray(a, n);
        setPermutationRow(a, n);
        return;
    }
 
    for (int i = 0; i < size; i++) {
        heapPermutation(a, size - 1, n);
 
        // if size is odd, swap 0th i.e (first) and 
        // (size-1)th i.e (last) element
        if (size % 2 == 1)
            swap(a[0], a[size - 1]);
 
        // If size is even, swap ith and 
        // (size-1)th i.e (last) element
        else
            swap(a[i], a[size - 1]);
    }
}

void getAllPermutattions(int a[vr-1], int n) {
	do {
		setPermutationRow(a, n);
	}while (next_permutation(a, a+n));
}

void generatePermutationArray(int p) {
	
	int ver[vr-1], j = 0;

	lastRow = 0;
	numPerm = 1; 

	for (int i = 0; i < vr; i++)
	{
		if (i != p)
		{
			ver[j++] = i;
		}

		if (i >= 2)
		{
			numPerm *= i;
		}

	}
	permutationArray = (int *)malloc(numPerm * (vr-1) * sizeof(int)); 
	// cout << "num: " << numPerm << endl;
	// heapPermutation(ver ,vr-1, vr-1);
	getAllPermutattions(ver, vr-1);
	// cout << "Permutation table done!" << endl;
}

int getMin(int x, int y) {
	if (x < y)
	{
		return x;
	}
	return y;
}

int TSP_parallel(int graph[][vr], int origin, int permutationCount) // implement traveling Salesman Problem. 
{
	int m_p = INT_MAX; // store minimum weight of a graph 

	omp_set_num_threads(2);
	#pragma omp parallel shared(graph)
	{
		// #pragma omp single
		// cout << "threads " << omp_get_num_threads() <<endl;

		#pragma omp for reduction(min: m_p) 
		for (long permId = 0; permId < permutationCount; ++permId)
		{
			int cur_pth = 0;
			int k = origin;
			int i;

			for (i = 0; i < vr-1; i++) {
				if (i == 0)
				{				
					// cur_pth += graph[0][ver[i]];
					cur_pth += graph[0][*(permutationArray + permId*(vr-1) + i)];
					if (i == (vr-2))
					{					
						cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][origin];
					}
				} else {
					// cur_pth += graph[ver[i-1]][ver[i]];
					cur_pth += graph[*(permutationArray + permId*(vr-1) + (i-1))][*(permutationArray + permId*(vr-1) + i)];
					if (i == (vr-2))
					{					
						cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][origin];
					}
				}
			}

			m_p = getMin(m_p, cur_pth); // to update the value of minimum weight
		}
	}

	return m_p;
}

void TSP_Sequential(int graph[][vr], int origin) // implement traveling Salesman Problem. 
{
	int m_p = INT_MAX; // store minimum weight of a graph 

	for (long permId = 0; permId < numPerm; ++permId)
	{
		int cur_pth = 0;
		int k = origin;
		int i;

		for (i = 0; i < vr-1; i++) {
			if (i == 0)
			{				
				// cur_pth += graph[0][ver[i]];
				cur_pth += graph[0][*(permutationArray + permId*(vr-1) + i)];
				if (i == (vr-2))
				{					
					cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][origin];
				}
			} else {
				// cur_pth += graph[ver[i-1]][ver[i]];
				cur_pth += graph[*(permutationArray + permId*(vr-1) + (i-1))][*(permutationArray + permId*(vr-1) + i)];
				if (i == (vr-2))
				{					
					cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][origin];
				}
			}
		}
		m_p = getMin(m_p, cur_pth); // to update the value of minimum weight
	}
	// cout << "Min Path(Parallel) " << m_p endl;
}

// organize block
void organizeBlockSize(int pid, int numPerm, int procCount, int &blockSize, unsigned long long int &idx, int colSize) {
	// assigns more work to some threads if n%p != 0
    // processor id 0 to id n%p-1 will do one extra row element 
    if (pid < numPerm%procCount) {
        blockSize++;
    }
    // starting point changes if n%p != 0
    // idx is the starting index of each subarray
    if (numPerm%procCount-pid > 0 && pid != 0)
    {
    	idx += pid*colSize;
    } else if(pid != 0) {
    	idx += numPerm%procCount*colSize;
    }
    // cout << blockSize << "*****inside " << idx << endl;

}


int main(int argc, char *argv[]) 
{

	int graph[vr][vr];
	setGraph(graph);
	int origin = 0, id, procCount, blockSize, minPath, tmpMin;
	unsigned long long int idx = 0;

	double time, time2;

	MPI_Status status;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &procCount);

	if (id == 0)
	{
		// for (int i = 0; i < vr; ++i)
		// {
		// 	for (int j = 0; j < vr; ++j)
		// 	{
		// 		cout << graph[i][j] << " ";
		// 	}
		// 	cout << endl;
		// }

		generatePermutationArray(origin);
		// for (int i = 0; i <  numPerm; i++) {
		// 	for (int j = 0; j < vr-1; j++) {
		// 		printf("%d ", *(permutationArray + i*(vr-1) + j)); 
		// 	}
		// 	cout << endl;
		// }

		time = MPI_Wtime();
		TSP_Sequential(graph, origin);	
		time  = MPI_Wtime() - time;
		// cout << "t_s: " << time << endl;
		// cout << endl;

		time2 = MPI_Wtime();		
		blockSize = numPerm / procCount;
		// broadcast blocks to each processor
		for (int pid = 1; pid < procCount; ++pid)
		{
			idx = (unsigned long long int) pid*blockSize*(vr-1);

			organizeBlockSize(pid, numPerm, procCount, blockSize, idx, (vr-1));
		    
		    // cout << idx << " -- " << blockSize << " " << pid << endl;

    		MPI_Send(&numPerm, 1, MPI_INT,pid,1, MPI_COMM_WORLD);			
    		MPI_Send((permutationArray+idx), (vr-1)*blockSize, MPI_INT,pid,1, MPI_COMM_WORLD);			
    		MPI_Send(graph, vr*vr, MPI_INT,pid,1, MPI_COMM_WORLD);			
			blockSize = numPerm / procCount;			
		}
		// cout << "Broadcasting over" << endl;

		idx = 0;
		organizeBlockSize(0, numPerm, procCount, blockSize, idx, (vr-1));
	    // cout << 0 << " -- " << blockSize << " " << 0 << endl;
	    tmpMin = TSP_parallel(graph, origin, blockSize);

	    // cout << "Master computation done!" << endl;

	} else {
	    MPI_Recv(&numPerm, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	
		blockSize = numPerm / procCount;
		idx = 0;
		organizeBlockSize(id, numPerm, procCount, blockSize, idx, (vr-1));
		permutationArray = (int *)malloc(blockSize * (vr-1) * sizeof(int)); 


	    MPI_Recv(permutationArray, (vr-1)*blockSize, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	
	    MPI_Recv(graph, vr*vr, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	
	    tmpMin = TSP_parallel(graph, origin, blockSize);
	}
	// MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&tmpMin, &minPath, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

    if (id == 0)
    {
    	// cout << "Parallel Result is: " << minPath << endl;
		time2  = MPI_Wtime() - time2;
		// cout << "t_p: " << time2 << endl;
		cout << vr << " " << procCount << " " << time << " " << time2 << endl;
    }

	// cleanup 
	MPI_Finalize();

	return 0;
}