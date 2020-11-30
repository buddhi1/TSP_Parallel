#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

#define vr 13

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

void generatePermutationArray(int origin) {
	
	int ver[vr-1], j = 0;

	lastRow = 0;

	for (int i = 0; i < vr; i++)
	{
		if (i != origin)
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

void TSP_Parallel(int graph[][vr], int origin, int numThreads) // implement traveling Salesman Problem. 
{
	int m_p = INT_MAX; // store minimum weight of a graph 
	omp_set_num_threads(numThreads);

	#pragma omp parallel for reduction(min: m_p) shared(graph) 
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
		m_p = min(m_p, cur_pth); 
		// cout << "Min Path(Parallel) " << m_p << endl;
	}
}

void TSP_Sequential(int graph[][vr], int origin) // implement traveling Salesman Problem. 
{
	int m_p = INT_MAX; // store minimum weight of a graph 

	// #pragma omp parallel for shared(graph)
	for (long permId = 0; permId < numPerm; ++permId)
	{
		int cur_pth = 0;
		int k = origin;
		int i;

		// #pragma omp parallel for reduction(+: cur_pth)
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

		m_p = min(m_p, cur_pth); 
	}
	// cout << "Min Path(Sequential) " << m_p << endl; 
}


int main() 
{
	int graph[vr][vr];
	setGraph(graph);
	int origin = 0;

	double t_s, t_p, t_p2;

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

	t_p = omp_get_wtime();	
	TSP_Sequential(graph, origin);	
	t_p = omp_get_wtime() - t_p;
	// cout << "t_s=" << t_p << endl;

	for (int numThreads = 2 ; numThreads < 16; ++numThreads)
	{	
		t_p2 = omp_get_wtime();	
		TSP_Parallel(graph, origin, numThreads);	
		t_p2 = omp_get_wtime() - t_p2;
		// cout << "t_p=" << t_p << endl;
		cout << vr << " " << numThreads << " " << t_p << " " << t_p2 << endl;
	}

	return 0;
}