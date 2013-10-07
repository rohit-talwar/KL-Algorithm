/*

Input format : 
No_of_Vertices No_of_Edges
#Line 1 Vertex_a Vertex_b Weight_of_Edge
... Till number of edges


Author - Rohit Talwar
*/
#include<cstdio>
#include<iostream>
#include<queue>
#include<vector>
#include<list>
#include<string>
#include<climits>
#include<stack>
#include<utility>
#include<deque>
#include<bitset>
#include<set>
#include<map>
#include<cstring>
#include<cstdlib>
#include<algorithm>

#define pii pair<int, int>
#define PB push_back

using namespace std;
typedef long long LL;

int n;    					// no of nodes or vertices
vector< pii > adj[100] ;	 		// adjacency list
int grp[100];					// Stores the group name of the ith vertex
vector<int> gval;				// the max g values obtained after each iteration
vector<int> swapA, swapB;			// the vertices from A and B to be swapped

vector< pii > D[2] ;		// there are two grps -- need to change the initial partition to accomodate many grps

bool comp(pii a, pii b){
	return a.second > b.second ;
}

void printPar(){
	for(int i=1; i<=n ; ++i)
		printf("vertex %d is in %d group\n", i, grp[i]);
}

void printAdjList(void){
	for(int i=1; i<=n ; ++i){
		printf("For vertex %d\n", i);
		for(vector< pii >::iterator it = adj[i].begin() ; it!=adj[i].end() ; ++it){
			printf("vertex = %d Weight = %d\n", it->first, it->second);
		}
	}
}

void initialPart(){
	for(int i=0; i<=n ; ++i)
		grp[i] = 0;
	int i=0;
	while(i < (n>>1)){
		int num = rand()%n + 1;
		if(grp[num] != 1){
			grp[num] = 1;
			++i;
		}
	}
	
/*	grp[0] = 0;
	grp[1] = 1;
	grp[2] = 0;
	grp[3] = 1;
	grp[4] = 0;
*/	for(int i=1 ; i<=n; ++i)
		printf("Vertex %d Group %d\n", i, grp[i]);
}

void Cal_D_Val(){  			// extensible to many grps also
	D[0].clear();
	D[1].clear();
	for(int i=1; i<=n ; ++i){
		int base = grp[i] ;
		int external = 0, internal = 0;
		for(vector< pii >::iterator it = adj[i].begin() ; it!=adj[i].end() ; ++it){
			if(grp[it->first] != base){
				external += it->second ;
			} else{
				internal += it->second ;
			}
		}
		D[base].PB(pii(i, external-internal) );
	}
}

void printD(){
	for(int i=0; i<=1 ; ++i){
		printf("For i = %d\n", i);
		for(vector< pii >::iterator it = D[i].begin() ; it!=D[i].end() ; ++it){
			printf("vertex = %d D Val = %d\n", it->first, it->second);
		}
	}
}

void sort_D(){
	sort(D[0].begin(), D[0].end(), comp);
	sort(D[1].begin(), D[1].end(), comp);
	printf("INSIDE SORT_D \n");
	printD();
}

void updateD(int a, int b){
	int a_edge[100]={0}, b_edge[100]={0} ;

	for(vector< pii >::iterator it = adj[a].begin() ; it!=adj[a].end() ; ++it)
		a_edge[it->first] = it->second;
	
	for(vector< pii >::iterator it = adj[b].begin() ; it!=adj[b].end() ; ++it)
		b_edge[it->first] = it->second;
	
	for(vector< pii >::iterator it = D[0].begin() ; it!=D[0].end() ; ++it)
		it->second = it->second + 2*a_edge[it->first] - 2*b_edge[it->first] ;
		
	for(vector< pii >::iterator it = D[1].begin() ; it!=D[1].end() ; ++it)
		it->second = it->second + 2*b_edge[it->first] - 2*a_edge[it->first] ;

	printf("INSIDE UPDATE_D \n");
	printD();
	
}

int nextSwap(){
	int maxg = -9999, a=-1, b=-1;
	for(vector< pii >::iterator ita = D[0].begin() ; ita!=D[0].end() ; ++ita){
		for(vector< pii >::iterator itb = D[1].begin() ; itb!=D[1].end() ; ++itb){
			if(maxg > (ita->second + itb->second)){
				break;
			}else{
				int Cab=0;
				for(vector< pii >::iterator it = adj[ita->first].begin(); it!=adj[ita->first].end() ; ++it){
					if(it->first == itb->first){
						Cab = it->second;
					printf("Cross edge found between %d and %d of wt %d\n",ita->first, it->first, it->second);
						break;
					}
				}
				int newg = ita->second + itb->second - 2*Cab ;
				cout<<"NEWG "<<newg<<endl;
				if(newg > maxg ){
					maxg = newg;
					a = ita->first;
					b = itb->first;
				}
			}
		}
	}
	if(a==-1 && b==-1){
		printf(" No further iterations possible\n Size of D[0] = %d\n", D[0].size() );
		return 0;
	}
	// end of one iteration
	swapA.PB(a);
	swapB.PB(b);
	gval.PB(maxg);
	cout<<"a "<<a<<" b "<<b<<endl;
	// remove a and b from D
	for(vector< pii >::iterator it = D[0].begin() ; it!=D[0].end() ; ++it){
		if(it->first == a){
			D[0].erase(it);
			break;
		}
	}
	for(vector< pii >::iterator it = D[1].begin() ; it!=D[1].end() ; ++it){
		if(it->first == b){
			D[1].erase(it);
			break;
		}
	}
	updateD(a, b);
	return 1;
}

int maximiseSwap(){
	int tempSum = 0, maxSum = -999999, k=-1;
	for(int i=0; i<gval.size() ; ++i){
		tempSum += gval[i];
		if(tempSum > maxSum){
			k = i+1;
			maxSum = tempSum;
		}
	}
	printf("maxSum %d\n", maxSum);
	if(maxSum <= 0){				// no more iterations required
		return 0; 
	}else{
		// Swap the elements till k!
		for(int i=0 ; i<k ; ++i){
			grp[swapA[i]] = 1;
			grp[swapB[i]] = 0;
		}
		// printing partition
		printPar();
	}
	swapA.clear();
	swapB.clear();
	gval.clear();
	return 1;
}

int main(){
	
	scanf("%d", &n);
	int edges;
	scanf("%d", &edges);
	for(int i=0 ; i<edges ; ++i){
		int a,b,w;
		scanf("%d %d %d", &a, &b, &w);
		adj[a].PB(pii(b, w));
		adj[b].PB(pii(a, w));
	}
	//printAdjList();
	initialPart();
	while(1){
		Cal_D_Val();
		sort_D();
		while(nextSwap()) ;
		if(maximiseSwap() ==0)
			break;
	}
	printf("FINAL ANSWER\n");
	printPar();
	return 0;
}
