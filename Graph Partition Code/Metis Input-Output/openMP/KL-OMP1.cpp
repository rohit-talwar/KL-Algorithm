/*
	Input format : Refer to the METIS manual version 4.0 : Chapter 4.5  
	
	The graphs given in the input folder adhere to the above format

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
#include<omp.h>

#define pii pair<int, int>
#define PB push_back
#define MAXSIZE 100000
using namespace std;
typedef long long LL;

#define CHUNK 100
int n, m;    					// no of nodes or vertices and edges
vector< pii > adj[MAXSIZE] ;	 		// adjacency list
int grp[MAXSIZE];					// Stores the group name of the ith vertex
vector<int> gval;				// the max g values obtained after each iteration
vector<int> swapA, swapB;			// the vertices from A and B to be swapped

vector< pii > D[2] ;			// there are two grps -- need to change the initial partition to accomodate many grps
vector<int> vert_wt;				// the weights of the vertices

bool comp(pii a, pii b){
	return a.second > b.second ;
}

void printPar(){
	for(int i=1; i<=n ; ++i){
	//	printf("vertex %d is in %d group\n", i, grp[i]);
		printf("%d\n",grp[i]);
	}
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
/*	
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
*/
	for(int i=1 ; i<n/2+1 ; ++i)
		grp[i] = 0 ;
	for(int i=n/2+1 ; i<=n ; ++i)
		grp[i] = 1 ;
//	for(int i=1 ; i<=n; ++i)
//		printf("Vertex %d Group %d\n", i, grp[i]);
}

void Cal_D_Val(){           // extensible to many grps also
	D[0].clear();
	D[1].clear();
	int external = 0, internal = 0;
	int i,j ;
#pragma omp parallel for private(external, internal, i,j) schedule(dynamic, CHUNK)
	for(i=1; i<=n ; ++i){
		external = 0;
		internal = 0 ;
		for(j=0; j<adj[i].size() ; ++j){
			if(grp[adj[i][j].first] != grp[i]){
				external += adj[i][j].second ;
			} else{
				internal += adj[i][j].second ;
			}
		}
#pragma omp critical 
		{
		D[grp[i]].PB(pii(i, external-internal) );   // can replace this with an array which can store the latest
		}	
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
#pragma omp parallel 
	{
#pragma omp sections nowait
		{

#pragma omp section
	sort(D[0].begin(), D[0].end(), comp);

#pragma omp section
	sort(D[1].begin(), D[1].end(), comp);
//	printf("INSIDE SORT_D \n");
//	printD();
		}
	}
}

void updateD(int a, int b){
	int a_edge[MAXSIZE]={0}, b_edge[MAXSIZE]={0} ;
	int i,j,k,l;
#pragma omp parallel shared(a_edge, b_edge)
	{
#pragma omp for private(i) schedule(dynamic,CHUNK) nowait				// can do a section paralllelization but
	for(i=0 ; i < adj[a].size() ; ++i)				// wont be beneficial as number of cores is more
		a_edge[adj[a][i].first] = adj[a][i].second;


#pragma omp for private(j) schedule(dynamic,CHUNK) nowait				// can do a section paralllelization but
	for(j=0 ; j<adj[b].size() ; ++j)
		b_edge[adj[b][j].first] = adj[b][j].second;

//#pragma omp barrier

#pragma omp for private(k) schedule(dynamic,CHUNK) nowait				
	for(k = 0 ; k<D[0].size() ; ++k)
		D[0][k].second = D[0][k].second + 2*a_edge[D[0][k].first] - 2*b_edge[D[0][k].first] ;
		
#pragma omp for private(k) schedule(dynamic,CHUNK) nowait				
	for(l=0; l<D[1].size() ; ++l)
		D[1][l].second = D[1][l].second + 2*a_edge[D[1][l].first] - 2*b_edge[D[1][l].first] ;

//	printf("INSIDE UPDATE_D \n");
//	printD();
	}
}

int nextSwap(){
	int maxg = -9999999, a=-1, b=-1;
	int i,j,k,l, Cab, newg ;
#pragma omp parallel for shared(maxg,a,b) private(i,j,k,Cab, newg) schedule(dynamic, CHUNK) collapse(2) 
	for(i = 0 ; i<D[0].size() ; ++i) {
		for(j=0 ; j< D[1].size() ;  ++j){
			if(maxg > (D[0][i].second + D[1][j].second)){
				continue ;
			}else{
				Cab=0;
				for(k=0 ; k<adj[D[0][i].first].size(); ++k ){
					if(adj[D[0][i].first][k].first == D[1][j].first){
						Cab = D[1][j].second;
				//printf("Cross edge found between %d and %d of wt %d\n",ita->first, it->first, it->second);
					//	cout<<endl;
						break;
					}
				}
				newg = D[0][i].second + D[1][j].second - 2*Cab ;
				//cout<<"NEWG "<<newg<<endl;
#pragma omp critical
				{
				if(newg > maxg ){
					
					maxg = newg;
					a = D[0][i].first ;
					b = D[1][j].first ; 
					
				}
				}
			}
		}
	}
	if(a==-1 && b==-1){
		//printf(" No further iterations possible\n Size of D[0] = %d\n", D[0].size() );
		return 0;
	}
	// end of one iteration
	swapA.PB(a);
	swapB.PB(b);
	gval.PB(maxg);
	//cout<<"a "<<a<<" b "<<b<<endl;
	// remove a and b from D
	
#pragma omp parallel for private(i) schedule(guided) 	
	for(i =0 ; i< D[0].size() ; ++i){
		if(D[0][i].first == a){
			D[0].erase(D[0].begin() + i);
		}
	}

#pragma omp parallel for private(i) schedule(guided) 	
	for(i = 0;  i<D[1].size() ; ++i){
		if(D[1][i].first == b){
			D[1].erase(D[1].begin() + i);
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
//	printf("maxSum %d\n", maxSum);
	if(maxSum <= 0){				// no more iterations required
		return 0; 
	}else{
		// Swap the elements till k!
		for(int i=0 ; i<k ; ++i){
			grp[swapA[i]] = 1;
			grp[swapB[i]] = 0;
		}
		// printing partition
		//printPar();
	}
	swapA.clear();
	swapB.clear();
	gval.clear();
	return 1;
}

char line[1000000000];						// the read line

int getNextNum(int pos, int &num){				// gets the number in integer form & stores into num	
	char a[15];						// Returns the pointer position pos where 
	int i = 0;						// to read for getting the next number
	while(line[pos] != '\n' && line[pos] != EOF && line[pos] != ' ' && line[pos] != '\0'){
		a[i++] = line[pos++];
	}
	a[i] = '\0';
	if(i!=0){
		num = atoi(a);
		++pos;
	}
	return pos;
}

void takeInp(char inpFile[100]){

	FILE *inp;
	inp = fopen(inpFile, "r");
	fscanf(inp," %[^\n]",line);
	int pos = 0, mode=0, len;	
	pos = getNextNum(pos, n);
	pos = getNextNum(pos, m);
	pos = getNextNum(pos, mode);
//	cout<<" Vertices "<<n<<" Edegs "<<m<<" Mode "<<mode<<endl;
	if(mode==0){					// Only vertices will be in input: vertices & edges have no weights
		int a = 0;				// the vertice in concern
		while(fscanf(inp," %[^\n]",line)!=EOF){
			++a;				// Starting vertex is = 1 in METIS
			pos = 0;
			len = strlen(line);
			while(pos < len){
				int nbr; 		// Neighbouring vertex
				pos = getNextNum(pos, nbr);
				adj[a].PB(pii(nbr , 1));
				adj[nbr].PB(pii(a, 1));
			}
		}
	}else if(mode==1){				// Edges also have weights!
		int a = 0;				// the vertice in concern
		while(fscanf(inp," %[^\n]",line)!=EOF){
			++a;				// Starting vertex is = 1 in METIS
			pos = 0;
			len = strlen(line);
			while(pos < len){
				int nbr, wt; 		// Neighbouring vertex
				pos = getNextNum(pos, nbr);
				pos = getNextNum(pos, wt);
				adj[a].PB(pii(nbr , wt));
				adj[nbr].PB(pii(a, wt));
			}
		}
	}else if(mode==10){				// vertice has a weight
		int a = 0;				// the vertice in concern
		while(fscanf(inp," %[^\n]",line)!=EOF){
			int nbr, wt; 		// Neighbouring vertex
			++a;				// Starting vertex is = 1 in METIS
			pos = 0;
			len = strlen(line);
			pos = getNextNum(pos, wt);

			while(pos < len){
				pos = getNextNum(pos, nbr);
				adj[a].PB(pii(nbr , 1));
				adj[nbr].PB(pii(a, wt));
			}
		}
	}else if(mode==11){
	}

	fclose(inp);
}
int crossingedge ;

int calc_Cutset(){
	crossingedge = 0;
	/*
	for(int i=1 ; i<=n ; ++i){
		int base1 = grp[i] ;
		for(int j=i+1 ; j<=n ; ++j){
			int base2 = grp[j] ;
			if(base1 != base2){
				crossingedge += 1 ;
			}
		}
	}
	*/
	for(int i=1 ; i<=n ; ++i){
		for(int j=0; j<adj[i].size() ; ++j){
			if(grp[adj[i][j].first] != grp[i]){
				crossingedge += adj[i][j].second ;
			}
		}
	}
	crossingedge /= 2 ;
}

int main(){
	crossingedge = 0;
	char inpFile[100];
	scanf("%s",inpFile);
	takeInp(inpFile);
	//printAdjList();
	initialPart();
	calc_Cutset() ;
	printf("crossing edge %d\n", crossingedge) ;
	while(1){
		Cal_D_Val();
		sort_D();
		while(nextSwap()) ;
		if(maximiseSwap() ==0)
			break;
	}
	printf("FINAL ANSWER\n");
	calc_Cutset() ;
	printf("crossing edge %d\n", crossingedge) ;
	printPar();

	return 0;
}
