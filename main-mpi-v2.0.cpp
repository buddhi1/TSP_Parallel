#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>

using namespace std;

#define vr 11

vector<int> vertices(vr-1);
long numPerm = 1; 

int numThreads = 1;

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

void setNumPermutations() {
	numPerm = 1; 

	for (int i = 1; i < vr; i++)
	{
		numPerm *= i;
	}
	// cout << "Permutations: " << numPerm << endl;
}

void convertNumToFactoradic(unsigned long long int num, int output[]) {

	int x;
	for (int i = 1, j = vr-2; i < vr; ++i, j--)
	{
		if (num != 0)
		{
			x = num % i;
			num /= i;
			output[j] = x;
		} else {
			output[j] = 0;			
		}
		// cout << "---==" << output[j] << endl;
	}
}

void getPermutation(vector<int> fact, int permutation[], unsigned long long int num) {
	int factoradic[vr-1], idx;

	if (num ==0)
	{
		for (int i = 0; i < vr-1; ++i)
		{
			permutation[i] = vertices[i];			
		}
	} else {
		convertNumToFactoradic(num, factoradic);
		for (int i = 0; i < vr-1; ++i)
		{
			idx = factoradic[i];
			permutation[i] = fact[idx];
			fact.erase(fact.begin()+idx);
			// cout << factoradic[i] << " *" << fact[factoradic[i]] << endl;
		}
	}
}

int getMin(int x, int y) {
	if (x < y)
	{
		return x;
	}
	return y;
}

int TSP_parallel(int graph[][vr], int origin, int permutationCount, int st) 
{
	int m_p = INT_MAX; // store minimum weight of a graph 

	omp_set_num_threads(numThreads);
	#pragma omp parallel shared(graph)
	{
		// #pragma omp single
		// cout << "threads " << omp_get_num_threads() <<endl;

		#pragma omp for reduction(min: m_p) 
		for (long permId = 0; permId < permutationCount; ++permId)
		{
			int cur_pth = 0;
			int k = origin;
			int i, permutation[vr-1];

			getPermutation(vertices, permutation, st+permId);	
			// for (int i = 0; i < vr-1; ++i)
			// {
			// 	cout << "+++" << permutation[i] << " " ;		
			// }
			// cout  << endl;


			for (i = 0; i < vr-1; i++) {
				if (i == 0)
				{				
					cur_pth += graph[0][permutation[i]];
					// cur_pth += graph[0][*(permutationArray + permId*(vr-1) + i)];
					if (i == (vr-2))
					{					
						cur_pth += graph[permutation[vr-2]][origin];
					}
				} else {
					cur_pth += graph[permutation[i-1]][permutation[i]];
					// cur_pth += graph[*(permutationArray + permId*(vr-1) + (i-1))][*(permutationArray + permId*(vr-1) + i)];
					if (i == (vr-2))
					{					
						cur_pth += graph[permutation[vr-2]][origin];
					}
				}
			}
			// cout << "&&& " << permutationCount << " " << cur_pth << " " << st+permId << endl;

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
		int i, permutation[vr-1];

		getPermutation(vertices, permutation, permId);	
		// cout << "-- " << endl;
		// for (int i = 0; i < vr-1; ++i)
		// {
		// 	cout << "+++" << permutation[i] << " " ;		
		// }
		// cout  << endl;	

		for (i = 0; i < vr-1; i++) {
			if (i == 0)
			{				
				cur_pth += graph[0][permutation[i]];
				if (i == (vr-2))
				{					
					cur_pth += graph[permutation[(vr-2)]][origin];
				}
			} else {
				cur_pth += graph[permutation[i-1]][permutation[i]];
				if (i == (vr-2))
				{					
					cur_pth += graph[permutation[(vr-2)]][origin];
				}
			}
		}

		m_p = getMin(m_p, cur_pth); // to update the value of minimum weight
	}
	// cout << "Min path(Sequential) " << m_p << endl;
}

// organize block
void organizeBlockSize(int pid, int numPerm, int procCount, int &blockSize, unsigned long long int &idx) {
	// assigns more work to some threads if n%p != 0
    // processor id 0 to id n%p-1 will do one extra row element 
    if (pid < numPerm%procCount) {
        blockSize++;
    }
    // starting point changes if n%p != 0
    // idx is the starting index of each subarray
    if (numPerm%procCount-pid > 0 && pid != 0)
    {
    	idx += pid;
    } else if(pid != 0) {
    	idx += numPerm%procCount;
    }
    // cout << blockSize << "*****inside " << idx << endl;

}


int main(int argc, char *argv[]) 
{

	int graph[vr][vr], origin = 0, id, procCount, blockSize;
	int tmpMin, minPath;
	unsigned long long int idx, tt;

	double time, time2;
	
	MPI_Status status;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
	MPI_Comm_size (MPI_COMM_WORLD, &procCount);

	numThreads =  atoi(argv[1]);

	
	for (int i = 1; i < vr; ++i)
	{
		vertices[i-1] = i;
	}

	if (id == 0)
	{
		setGraph(graph);

		// for (int i = 0; i < vr; ++i)
		// {
		// 	for (int j = 0; j < vr; ++j)
		// 	{
		// 		cout << graph[i][j] << " ";
		// 	}
		// 	cout << endl;
		// }

		setNumPermutations();

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
			idx = (unsigned long long int) pid*blockSize;

			organizeBlockSize(pid, numPerm, procCount, blockSize, idx);
		    
		    // cout << idx << " -- " << blockSize << " " << pid << endl;

    		MPI_Send(&numPerm, 1, MPI_INT, pid,1, MPI_COMM_WORLD);						
    		MPI_Send(&idx, 1, MPI_UNSIGNED_LONG_LONG, pid,1, MPI_COMM_WORLD);						
    		MPI_Send(graph, vr*vr, MPI_INT, pid,1, MPI_COMM_WORLD);			
			blockSize = numPerm / procCount;			
		}
		// cout << "Broadcasting over" << endl;

		idx = 0;
		organizeBlockSize(0, numPerm, procCount, blockSize, idx);
		// cout << "eee " << id << " " << idx << " " << blockSize << endl;

	    // cout << 0 << " -- " << blockSize << " " << 0 << endl;
	    tmpMin = TSP_parallel(graph, origin, blockSize, idx);

	    // cout << "Master computation done!" << endl;

	} else {
	    MPI_Recv(&numPerm, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	
		blockSize = numPerm / procCount;
	    MPI_Recv(&idx, 1, MPI_UNSIGNED_LONG_LONG, 0, 1, MPI_COMM_WORLD, &status);	
		
		organizeBlockSize(id, numPerm, procCount, blockSize, tt); 

		// cout << "eee " << id << " " << idx << " " << blockSize << endl;
	    MPI_Recv(graph, vr*vr, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);	
	    tmpMin = TSP_parallel(graph, origin, blockSize, idx);
	}
	// MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&tmpMin, &minPath, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

    if (id == 0)
    {
    	// cout << "Parallel Result is: " << minPath << endl;
		time2  = MPI_Wtime() - time2;
		// cout << "t_p: " << time2 << endl;
		cout << vr << " " << procCount << " " << numThreads << " " << time << " " << time2 << endl;
    }

	// cleanup 
	MPI_Finalize();

	return 0;
}