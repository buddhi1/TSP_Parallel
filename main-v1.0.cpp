#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

#define vr 5

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

int TSP(int grph[][vr], int p) // implement traveling Salesman Problem. 
{
	vector<int> ver; 
	int cd = 0;

	for (int i = 0; i < vr; i++)
	{
		if (i != p)
		{
			ver.push_back(i);
		}
	}
	int m_p = INT_MAX; // store minimum weight of a graph

	do 
	{
		cd++;
		int cur_pth = 0;
		int k = p;

		// for (int i = 0; i < ver.size(); i++) {
		for (int i = 0; i < vr-1; i++) {
			// cout << ver[i] << " ";
			cur_pth += grph[k][ver[i]];
			k = ver[i];
		}
		// cout << endl;
		cur_pth += grph[k][p];
		// cout << "cur_path: " << cur_pth << endl;

		m_p = min(m_p, cur_pth); // to update the value of minimum weight
	}
	while (next_permutation(ver.begin(), ver.end()));

	cout << "TSP Sequential: " << cd << endl;
	return m_p;
}

int TSP3(int graph[][vr], int p) // implement traveling Salesman Problem. 
{
	vector<int> ver; 
	int cd = 0, i;

	for (int i = 0; i < vr; i++)
	{
		if (i != p)
		{
			ver.push_back(i);
		}
	}
	int m_p = INT_MAX; // store minimum weight of a graph

	omp_set_num_threads(8);
	#pragma omp parallel shared(graph, ver)
	{
		do 
		{
			cd++;
			int cur_pth = 0;
			int k = p;

			// #pragma omp parallel for reduction(+: cur_pth)
			for (i = 0; i < vr-1; i++) {
				if (i == 0)
				{				
					cur_pth += graph[0][ver[i]];
				} else {
					cur_pth += graph[ver[i-1]][ver[i]];
				}
			}
			cur_pth += graph[ver[ver.size()-1]][p];

			m_p = min(m_p, cur_pth); // to update the value of minimum weight
		}
		while (next_permutation(ver.begin(), ver.end()));
	}
	cout << "TSP33333: " << cd << endl;
	return m_p;
}

void generatePermutationArray(int p) {
	
	int ver[vr-1], j = 0;

	lastRow = 0;

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
	cout << "num: " << numPerm << endl;
	// heapPermutation(ver ,vr-1, vr-1);
	getAllPermutattions(ver, vr-1);
	cout << "Permutation table done!" << endl;
}

int TSP2(int graph[][vr], int p) // implement traveling Salesman Problem. 
{
	int m_p = INT_MAX; // store minimum weight of a graph 

	#pragma omp parallel for shared(graph)
	for (long permId = 0; permId < numPerm; ++permId)
	{
		int cur_pth = 0;
		int k = p;
		int i;

		// #pragma omp parallel
		{

			// #pragma omp for reduction(+: cur_pth)
			for (i = 0; i < vr-1; i++) {
				if (i == 0)
				{				
					// cur_pth += graph[0][ver[i]];
					cur_pth += graph[0][*(permutationArray + permId*(vr-1) + i)];
					if (i == (vr-2))
					{					
						cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][p];
					}
				} else {
					// cur_pth += graph[ver[i-1]][ver[i]];
					cur_pth += graph[*(permutationArray + permId*(vr-1) + (i-1))][*(permutationArray + permId*(vr-1) + i)];
					if (i == (vr-2))
					{					
						cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][p];
					}
				}
			}
		}
		// cur_pth += graph[ver[ver.size()-1]][p];
		// cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][p];
		// cout << "cur_path2: " << cur_pth << endl;

		m_p = min(m_p, cur_pth); // to update the value of minimum weight

	}
	// delete(permutationArray);

	return m_p;
}

int TSP22(int graph[][vr], int p) // implement traveling Salesman Problem. 
{
	int m_p = INT_MAX; // store minimum weight of a graph 

	// #pragma omp parallel for shared(graph)
	for (long permId = 0; permId < numPerm; ++permId)
	{
		int cur_pth = 0;
		int k = p;
		int i;

		// #pragma omp parallel for reduction(+: cur_pth)
		for (i = 0; i < vr-1; i++) {
			if (i == 0)
			{				
				// cur_pth += graph[0][ver[i]];
				cur_pth += graph[0][*(permutationArray + permId*(vr-1) + i)];
				if (i == (vr-2))
				{					
					cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][p];
				}
			} else {
				// cur_pth += graph[ver[i-1]][ver[i]];
				cur_pth += graph[*(permutationArray + permId*(vr-1) + (i-1))][*(permutationArray + permId*(vr-1) + i)];
				if (i == (vr-2))
				{					
					cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][p];
				}
			}
		}
		// cur_pth += graph[ver[ver.size()-1]][p];
		// cur_pth += graph[*(permutationArray + permId*(vr-1) + (vr-2))][p];
		// cout << "cur_path2: " << cur_pth << endl;

		m_p = min(m_p, cur_pth); // to update the value of minimum weight

	}
	return m_p;
}


int main() 
{
	// int graph[][vr] = { { 0, 5, 10, 15 }, //values of a graph in a form of matrix
	// 					{ 5, 0, 20, 30 },
	// 					{ 10, 20, 0, 35 },
	// 					{ 15, 30, 35, 0 }
	// 					};
	int graph[vr][vr];
	setGraph(graph);
	int p = 0;

	double t_s, t_p, t_p2;

	for (int i = 0; i < vr; ++i)
	{
		for (int j = 0; j < vr; ++j)
		{
			cout << graph[i][j] << " ";
		}
		cout << endl;
	}

	generatePermutationArray(p);

	for (int i = 0; i <  numPerm; i++) {
		for (int j = 0; j < vr-1; j++) {
			printf("%d ", *(permutationArray + i*(vr-1) + j)); 
		}
		cout << endl;
	}

	t_s = omp_get_wtime();
	cout<< "Sequential result is: "<< TSP(graph, p) << endl;
	t_s = omp_get_wtime() - t_s;
	cout << "t_s=" << t_s << endl;

	t_p = omp_get_wtime();	
	cout<< "Parallel result is: "<< TSP2(graph, p) << endl;	
	t_p = omp_get_wtime() - t_p;
	cout << "t_p=" << t_p << endl;

	// t_p = omp_get_wtime();	
	// cout<< "Parallel result is: "<< TSP3(graph, p) << endl;	
	// t_p = omp_get_wtime() - t_p;
	// cout << "t_p=" << t_p << endl;

	return 0;
}