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

#define pii pair<int, int>
#define PB push_back
#define MAXSIZE 100000
using namespace std;
typedef long long LL;

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
//		printf("vertex %d is in %d group\n", i, grp[i]);
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
//	for(int i=1 ; i<=n; ++i)
//		printf("Vertex %d Group %d\n", i, grp[i]);
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
//	printf("INSIDE SORT_D \n");
//	printD();
}

void updateD(int a, int b){
	int a_edge[MAXSIZE]={0}, b_edge[MAXSIZE]={0} ;

	for(vector< pii >::iterator it = adj[a].begin() ; it!=adj[a].end() ; ++it)
		a_edge[it->first] = it->second;
	
	for(vector< pii >::iterator it = adj[b].begin() ; it!=adj[b].end() ; ++it)
		b_edge[it->first] = it->second;
	
	for(vector< pii >::iterator it = D[0].begin() ; it!=D[0].end() ; ++it)
		it->second = it->second + 2*a_edge[it->first] - 2*b_edge[it->first] ;
		
	for(vector< pii >::iterator it = D[1].begin() ; it!=D[1].end() ; ++it)
		it->second = it->second + 2*b_edge[it->first] - 2*a_edge[it->first] ;

//	printf("INSIDE UPDATE_D \n");
//	printD();
	
}

int nextSwap(){
	int maxg = -9999999, a=-1, b=-1;
	for(vector< pii >::iterator ita = D[0].begin() ; ita!=D[0].end() ; ++ita){
		for(vector< pii >::iterator itb = D[1].begin() ; itb!=D[1].end() ; ++itb){
			if(maxg > (ita->second + itb->second)){
				break;
			}else{
				int Cab=0;
				for(vector< pii >::iterator it = adj[ita->first].begin(); it!=adj[ita->first].end() ; ++it){
					if(it->first == itb->first){
						Cab = it->second;
				//printf("Cross edge found between %d and %d of wt %d\n",ita->first, it->first, it->second);
					//	cout<<endl;
						break;
					}
				}
				int newg = ita->second + itb->second - 2*Cab ;
				//cout<<"NEWG "<<newg<<endl;
				if(newg > maxg ){
					maxg = newg;
					a = ita->first;
					b = itb->first;
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
	return (maxg>0 ? 1 : 0 ) ;
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
	//printf("maxSum %d\n", maxSum);
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
	
	char inpFile[100];
	scanf("%s",inpFile);
	takeInp(inpFile);
//	printAdjList();
	initialPart();
	int ct=0 , ct1=0;
	while(ct<10){
		int ct1 = 0; 
		Cal_D_Val();
		sort_D();
		while(ct1<50){
			if(!nextSwap()){
				++ct1;
			}else{
				ct1 = 0;
			}
		}
		if(maximiseSwap() ==0)
			break;
	}
//	printf("FINAL ANSWER\n");
	calc_Cutset() ;
	printf("%d\n",crossingedge) ; 
	printPar();

	return 0;
}
