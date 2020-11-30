#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

#define vr 13

vector<int> vertices(vr-1);
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

void TSP_Parallel(int graph[][vr], int origin, int numThreads) // implement traveling Salesman Problem. 
{
	int m_p = INT_MAX; // store minimum weight of a graph 

	omp_set_num_threads(numThreads);
	#pragma omp parallel shared(graph)
	{
		// #pragma omp single
		// cout << "threads " << omp_get_num_threads() <<endl;

		#pragma omp for reduction(min: m_p) 
		for (long permId = 0; permId < numPerm; ++permId)
		{
			int cur_pth = 0;
			int k = origin;
			int i, permutation[vr-1];

			getPermutation(vertices, permutation, permId);	
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

	// cout << "Min Path(Parallel) " << m_p << endl;
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
	// cout << "Min Path(Parallel) " << m_p << endl;
}


int main() 
{
	int graph[vr][vr];
	int origin = 0;

	double t_s, t_p, t_p2;

	setGraph(graph);
	// for (int i = 0; i < vr; ++i)
	// {
	// 	for (int j = 0; j < vr; ++j)
	// 	{
	// 		cout << graph[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }

	for (int i = 1; i < vr; ++i)
	{
		vertices[i-1] = i;
	}

	setNumPermutations();

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

	for (int numThreads = 2; numThreads < 16; ++numThreads)
	{
		t_p2 = omp_get_wtime();	
		TSP_Parallel(graph, origin, numThreads);	
		t_p2 = omp_get_wtime() - t_p2;
		// cout << "t_p=" << t_p2 << endl;
		cout << vr << " " << numThreads << " " << t_p << " " << t_p2 << endl;
	}

	return 0;
}