#include <bits/stdc++.h>
using namespace std; 

const int INF = 1e5; 
const int nax = 6;

vector<int> adj[nax];
vector<pair<int,int>> grph[nax]; 
vector<bool> vis(nax,false);

void addEdge(int a, int b, int wt) {
	grph[a].push_back(make_pair(b,wt));
	// grph[b].push_back(make_pair(a,wt));
}

void addEdge(int a, int b) {
	adj[a].push_back(b); 
	// adj[b].push_back(a);
}
//---------------------------------------------------------------------------------------------------
vector<int> dist(nax,0), parent(nax,0); // to be used in BFS

void printBfs(int src) {
	dist[src] = 0; 
	parent[src] = -1; 

	queue<int> q; 
	q.push(src); 
	vis[src] = true; 

	while(!q.empty()) {
		int dest = q.front(); 
		q.pop(); 

		cout << dest << " " ; 
		for(int node : adj[dest]) {
			if(!vis[node]) {
				vis[node] = true; 
				q.push(node);
				dist[node] = dist[dest] + 1; 
				parent[node] = dest;
			}
		}
	}
}
//----------------------------------------------------------------------------------------------------
vector<int> time_in(nax,0), time_out(nax,0); 
int dfs_timer = -1; 

void printDfs(int src) {
	time_in[src] = ++dfs_timer;
	cout << src << " "; 

	for(int node : adj[src]) {
		if(!vis[node]) {
			vis[node] = true;
			printDfs(node);
		}
	}
	time_out[src] = ++dfs_timer;
}
//----------------------------------------------------------------------------------------------------
bool cycle_Detection_Undirected_Graph_BFS(int src) {
	queue<pair<int, int>> q; 
	vis[src] = true;
	q.push({src, -1});

	while(!q.empty()) {
		int node = q.front().first; 
		int parent = q.front().second; // to check the visited vertex is not the absolute ancestor
		q.pop(); 

		for(int it : adj[node]) {
			if(!vis[it]) {
				vis[it] = true; 
				q.push({it, node});
			} 
			else if(parent != it)
				return true; 
		}
	}
	return false; 
}
//----------------------------------------------------------------------------------------------------
bool cycle_Detection_Undirected_Graph_DFS(int src, int parent) { 
	vis[src] = true; 

	for(int node : adj[src]) {
		if(!vis[node]) { // impt : if deep down any recursive call returns true, then func becomes true; 
			if(cycle_Detection_Undirected_Graph_DFS(node, src))
				return true; 
			cycle_Detection_Undirected_Graph_DFS(node, src);
		} 
		else if(node != parent) return true; 
	}	
	return false; 
}
//----------------------------------------------------------------------------------------------------
// BIRPARTITE CHECK .
vector<int> color(nax, -1);  // color array will also be used as visited array. 

bool bipartite_BFS(int src) {	
	queue<int> q; 
	q.push(src); 
	color[src] = 1; 

	while(!q.empty()) {
		int node = q.front(); 
		q.pop(); 

		for(int it : adj[node]) {
			if(color[it] == -1) {
				color[it] = 1 - color[node];
				q.push(it);
			}
			else if(color[it] == color[node]) 
				return false; 
		}
	}
	return true; 	
}

bool bipartite_DFS(int src) {
	if(color[src] == -1) 	// color only when it is not colored before. 
		color[src] = 1; 

	for(int node : adj[src]) {
		if(color[node] == -1) {
			color[node] = 1 - color[src]; 
			if(!bipartite_DFS(node))
				return false; 
		}
		else if(color[node] == color[src])
			return false;
	}
	return true; 
}

bool checkBipartite() {
	for(int i = 1; i < nax; i++) {
		if(color[i] == -1) {
			// if(!bipartite_BFS(i)) { 	// any component fails to be bipartite. 
			if(!bipartite_DFS(i)) {
				return false; 
			}
		}
	}
	return true; 
}
//----------------------------------------------------------------------------------------------------
// CYCLE DETECTION DIRECTED GRAPH
vector<bool> currDfs(nax, false); 

bool cycle_Detection_Directed_Graph_DFS(int src) {
	vis[src] = true; 
	currDfs[src] = true; 

	for(int node : adj[src]){ 
		if(!vis[node]){ 
			if(cycle_Detection_Directed_Graph_DFS(node)) 
				return true; 
		} else if(currDfs[node]) {
			return true; 
		}
	}
	currDfs[src] = false; 
	return false;
}
//----------------------------------------------------------------------------------------------------
// TOPOLOGICAL SORT
void topologicalSort(int src, stack<int>& st) {
	vis[src] = true; 
	for(int node : adj[src]) {
		if(!vis[node]) {
			topologicalSort(node, st);
		}
	}
	st.push(src);
}
// USING BFS - KAHN'S ALGORITHM :- 

// start topo sort by pushing the nodes with in-degree 0, because they have a dependency of 0 (no vertex 
// before them). then use them to reduce the in-degree of their outgoing vertices by 1, means, reducing their 
// dependency by this current vertex. as soon as some vertex reaches in-degree 0, push it in the queue and 
// obviously this will be the next node in the topologicalSort.

void topologicalSort_KAHN() {
	vector<int> indegree(nax, 0); 

	for(int i = 0; i < nax; i++) {
		for(int node : adj[i]) { 
			indegree[node]++;
		}
	}

	queue<int> q; 
	for(int i = 0; i < nax; i++) {
		if(indegree[i] == 0) {
			q.push(i);
		}
	}

	while(!q.empty()) {
		int node = q.front(); 
		q.pop(); 
		cout << node << " "; 

		for(int to : adj[node]) {
			indegree[to]--;

			if(indegree[to] == 0) {
				q.push(to);
			}
		}
	}
}
//----------------------------------------------------------------------------------------------------
// CYCLE DETECTION IN DIRECTED GRAPHS USING BFS :

// since, a topologicalSort is only possible for DAG, so it can be concluded that if the topologicalSort
// is not possible to be generated for the given graph, then it is cyclic (in case of directed only). 
// so, if it is cyclic then the total count of the vertices produced in the topologicalSort will 
// not be equal to the vertices in the graph.
// so instead of cout << node <<" "; do :-> cnt++; at the end return (cnt == nax ? no : yes); 

//----------------------------------------------------------------------------------------------------
//  SHORTEST PATH IN UNDIRECTED GRAPH (USING BFS)

vector<int> shortest_Dist_Undirected(int src) {
	vector<int> dist(nax, INF); 

	dist[src] = 0; 
	queue<int> q; 
	q.push(src); 

	while(!q.empty()) {
		int node = q.front(); 
		q.pop(); 

		for(int to : adj[node]) {
			if(dist[to] > dist[node] + 1) {
				dist[to] = 1 + dist[node]; 
				q.push(to);
			}
		}
	}
	return dist;
}
//----------------------------------------------------------------------------------------------------
// SHORTEST PATH IN DIRECTED ACYCLIC GRAPH USING TOPOLOGICAL SORT. 
// TIME : 2 * O(V + E)

class shortest_Path_DAG {
	void findTopoSort(int src, stack<int>& st) {
		vis[src] = true; 
		for(auto it : grph[src]) {
			if(!vis[it.first]) {
				findTopoSort(it.first, st); 
			}
		}
		st.push(src); 
	}
public:
	vector<int> shortest_Paths(int src) {
		stack<int> st; 
		for(int i = 0; i < nax; i++) {
			if(!vis[i]) {
				findTopoSort(i, st);
			}
		}

		vector<int> dist(nax, INF); 
		dist[src] = 0; 

		while(!st.empty()) {
			int node = st.top(); 
			st.pop(); 
			
			if(dist[node] != INF) {	
			// this node is reached before, only then its neighbour vertices can be reached
				for(auto it : grph[node]) {
					if(dist[it.first] > dist[node] + it.second) {
						dist[it.first] = dist[node] + it.second;
					}
				}
			}
		}
		return dist;
	}
};
//----------------------------------------------------------------------------------------------------
// DIJKSTRA'S ALGORITHM 
// TIME : O((N+M)*log(N)) or O(N*logN)
vector<int> dijkstra(int src) {
	vector<int> distTo(nax, INF); 
	distTo[src] = 0; 

	priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> pq; 

	pq.push(make_pair(0,src)); 

	while(!pq.empty()) {
		int dist = pq.top().first; // dist = distTo[prev] 
		int prev = pq.top().second; 
		pq.pop(); 

		for(auto it : grph[prev]){ 
			int next = it.first; 
			int nextDist = it.second; // wt for this edge. 

			if(distTo[next] > distTo[prev] + nextDist) {	//	if(distTo[next] > dist + nextDist) {
				distTo[next] = distTo[prev] + nextDist; 	//		distTo[next] = dist + nextDist;
				pq.push(make_pair(distTo[next], next)); 	//		pq.push(make_pair(distTo[next], next));
			}												//	}
		}
	}
	return distTo; 
}
//----------------------------------------------------------------------------------------------------
// DISJOINT SET UNION 
int root[nax]; 
int siz[nax];

void dsu_make() {
	for(int i = 0; i < nax; i++) {
		root[i] = i; 
		siz[i] = 1; 
	}
}
int dsu_find(int x) {
	if(root[x] == x) return x; 
	return root[x] = dsu_find(root[x]);
}
void dsu_unite(int a, int b) {
	a = dsu_find(a); 
	b = dsu_find(b); 
	if(a != b) {
		if(siz[a] < siz[b]) 
			swap(a,b); 
		root[b] = a; 
		siz[a] += siz[b]; 
	}
}
bool dsu_same(int a, int b) {
	return dsu_find(a) == dsu_find(b); 
}

// below version of DSU is not tested
class DSU {
	vector<int> root, siz;
public:
	DSU(int _n) {
		root.resize(n);
		siz.resize(n);
		for(int i = 0; i < _n; i++) {
			root[i] = i; 
			siz[i] = 1; 
		}
	}
	int find(int a){
		while(root[a] != a) {
			root[a] = root[root[a]];
		}
		return a; 
	}
	void unite(int a, int b) {
		a = find(a);
		b = find(b);
		if(a != b) {
			if(siz[a] < siz[b]) {
				swap(a, b);
			}
			root[b] = a; 
			siz[a] += siz[b];
		}
	}
};

//----------------------------------------------------------------------------------------------------
// MINIMUM SPANNING TREE : KRUSKAL'S ALGORITHM 
int kruskal() {
	vector< pair< int, pair< int,int > > > egs; 	// AVOID USING THIS HERE. USE THIS IN MAIN()
	for(int i = 1; i < nax; i++) {
		for(int j = 0; j < grph[i].size(); j++) {
			egs.push_back(make_pair(grph[i][j].second, make_pair(i, grph[i][j].first)));
		}
	}
	sort(egs.begin(), egs.end());  vector<pair<int,int>> mst_egs; 

	dsu_make(); 
 	
	int mst_wt = 0; 
	for(auto it : egs) {
		int wt = it.first; 
		int u = it.second.first; 
		int v = it.second.second; 	

		if(!dsu_same(u,v)) {
			mst_wt += wt; 
			mst_egs.push_back(make_pair(u,v)); 
			dsu_unite(u,v); 
		}
	} 
	for(auto it : mst_egs) { 
		cout << it.first << " " << it.second << endl; 
	}
	return mst_wt;
}
//----------------------------------------------------------------------------------------------------
// BRIDGES IN OFFLINE USING DFS : O(N+M)
void bridgeDFS(int node, int parent, int timer, vector<int>& tin, vector<int>& low) {
	vis[node] = true; 
	tin[node] = low[node] = timer++; 

	for(int to : adj[node]) {
		if(to == parent) continue; // so that we don't back to the node we came from 

		else if(!vis[to]) {		// visit this and after its dfs finish, update the low[node]. low[node] is 
			bridgeDFS(to, node, timer, tin, low); // updated freqntly to mark the presence of some another path
			low[node] = min(low[node], low[to]); // from its neighbours to reach (to). and this value will be 
			// equal to the lowest tin[x], where x is the descendant of (node) in the dfs.  

			if(low[to] > tin[node]) {	// this means there is no back edge and all the ancestors of node 
				cout << to << " " << node << '\n';  // including itself have to be visited first, no backedge
			}			// from this (to) or its descendants to (node) or its ancestors. 
		}
		else { // this means for such visited (to), the tin[to] == low[to] == timer++ as this was the ancestor 
			low[node] = min(low[node], tin[to]); // of (node) in the dfs. 
		}
	}
}
void bridges() {
	vector<int> tin(nax,-1); 
	vector<int> low(nax,-1); 
	int timer = 0; 
	for(int i = 1; i < nax; i++) {
		if(!vis[i]) {
			bridgeDFS(i, -1, timer, tin, low); 
		}
	}
}
//----------------------------------------------------------------------------------------------------
// ARTICULATION POINTS 

void articulationDFS(int node, int parent, int timer, vector<int>& tin, vector<int>& low) {
	vis[node] = true; 
	tin[node] = low[node] = timer++; 

	for(int to : adj[node]) {
		if(to == parent) continue; 

		if(!vis[to]) {
			articulationDFS(to, node, timer, tin, low); 
			low[node] = min(low[node], low[to]); 

			if(low[to] >= tin[node] && parent != -1) {
				cout << node << " " << '\n'; 
			}
		} else {
			low[node] = min(low[node], tin[to]); 
		}
	}
}

void articulationPoints() {
	vector<int> tin(nax,-1); 
	vector<int> low(nax,-1); 
	int timer = 0; 
	for(int i = 1; i < nax; i++) {
		if(!vis[i]) {
			articulationDFS(i, -1, timer, tin, low); 
		}
	}
}
//----------------------------------------------------------------------------------------------------
// KOSARAJU'S ALGORITHM (STRONGLY CONNECTED COMPONENTS) 
class kosaraju{
	vector<int> transpose[nax];
	void topoSort(int src, stack<int>& st) {
		vis[src] = true; 
		for(int node : adj[src]) {
			if(!vis[node]) {
				topoSort(node, st); 
			}
		}
		st.push(src); 
	}
	void revDfs(int src) {
		vis[src] = true; 
		cout << src << " " ; 
		for(int to : transpose[src]) {
			if(!vis[to]) {
				revDfs(to); 
			}
		}
	}
public:
	void scc() {
		stack<int> st; 
		for(int i = 0; i < nax; i++) {
			if(!vis[i]) {
				topoSort(i, st); 
			}
		}

		for(int i = 1; i < nax; i++) {
			vis[i] = false; 
			for(int node : adj[i]) {
				transpose[node].push_back(i); 
			}
		}

		while(!st.empty()) {
			int node = st.top(); 
			st.pop(); 

			if(!vis[node]) {
				cout << "SCC: "; 
				revDfs(node); 
				cout << endl; 
			}
		}
	}
};
//----------------------------------------------------------------------------------------------------
// BELLMAN FORD ALGORITHM 
struct node {
	int u, v, wt; 
	node(int u, int v, int wt) {
		this->u = u; 
		this->v = v; 
		this->wt = wt; 
	}
}; 
vector<node> egs; 
void getEgs(int u, int v, int wt) {
	egs.push_back(node(u,v,wt)); 
}
vector<int> BellmanFord(int src) {
	// total nodes are stll nax. 
	vector<int> distTo(nax, INF); 
	distTo[src] = 0; 
// max. path that can be traversed in a graph is of length (n-1). so , relax every edge for (n-1) times. 
	for(int i = 0; i < nax - 1; i++) {
		for(auto it : egs) {
			if(distTo[it.v] > distTo[it.u] + it.wt) {
				distTo[it.v] = distTo[it.u] + it.wt; 
			}
		}
	}
	// do it one more time to find if the dist reduces then its a negative cycle. 
	for(auto it : egs) {
		if(distTo[it.v] > distTo[it.u] + it.wt) {
			cout << "Negative Cycle" << endl; 
			break; // found atleast one negative cycle. 
		}
	}
	return distTo;
}
//----------------------------------------------------------------------------------------------------

// DRIVER FUNCTION =============================================================================================
int32_t main() {
	ios::sync_with_stdio(0); cin.tie(0); cout.tie(0);
	
	// SAMPLE DIRECTED GRAPH

	// addEdge(1,2); 
	// addEdge(2,3); 
	// addEdge(3,4); 
	// addEdge(4,5);
	// addEdge(3,6);
	// addEdge(6,5);
	// addEdge(7,2);
	// addEdge(7,8);
	// addEdge(8,9);
	// addEdge(9,7);

	// addEdge(5,0); 
	// addEdge(5,2);
	// addEdge(4,0);
	// addEdge(4,1);
	// addEdge(2,3);
	// addEdge(3,1);

	// addEdge(1,3); 
	// addEdge(3,2); 
	// addEdge(2,1); 
	// addEdge(3,5); 
	// addEdge(5,4); 
	// addEdge(4,6); 
	// addEdge(6,5); 
	

	// ----------------------------------------------------------------------------------
	// SAMPLE UNDIRECTED GRAPH

	// addEdge(1,2);
	// addEdge(2,4);
	// addEdge(3,5);
	// addEdge(5,10);
	// addEdge(5,6);
	// addEdge(6,7);
	// addEdge(10,9);
	// addEdge(9,8);
	// addEdge(8,11);
	// addEdge(7,8);

	// addEdge(0,1);
	// addEdge(1,2);
	// addEdge(2,6);
	// addEdge(6,7);
	// addEdge(7,8);
	// addEdge(6,8);
	// addEdge(0,3);
	// addEdge(1,3);
	// addEdge(3,4);
	// addEdge(4,5);
	// addEdge(5,6);

	// addEdge(1,2); 
	// addEdge(2,3); 
	// addEdge(1,4);
	// addEdge(3,4);
	// addEdge(4,5);
	// addEdge(5,6);
	// addEdge(6,9);
	// addEdge(9,8);
	// addEdge(8,7);
	// addEdge(7,6);
	// addEdge(8,10);
	// addEdge(10,11);
	// addEdge(11,12);
	// addEdge(10,12);

	// -------------------------------------------------------------------------------------
	// 	SAMPLE WEIGHTED DAG
	// addEdge(1,2,2);
	// addEdge(2,5,5);
	// addEdge(1,4,1);
	// addEdge(4,3,3);
	// addEdge(3,2,4);
	// addEdge(3,5,1);

	// addEdge(5,4,9); 
	// addEdge(5,1,4);
	// addEdge(1,4,1);
	// addEdge(1,2,2);
	// addEdge(4,3,5);
	// addEdge(4,2,3);
	// addEdge(2,3,3);
	// addEdge(2,6,7);
	// addEdge(3,6,8);
	
	getEgs(0,1,5); 
	getEgs(1,5,-3); 
	getEgs(5,3,1); 
	getEgs(3,2,6); 
	getEgs(1,2,-2); 
	getEgs(2,4,3); 
	getEgs(3,4,-2); 


	// -------------------------------------------------------------------------------------

	// for(int i = 1; i < nax; i++) {
	// 	if(!vis[i])
	// 		printBfs(i);
	// } cout << endl; 

	// for(int i = 1; i < nax; i++) {
	// 	cout << "parent[" << i << "]: " << parent[i] << endl; 
	// 	cout << "dist[" << i << "]: " << dist[i] << endl; 
	// }

	// for(int i = 1; i < nax; i++) {
	// 	if(!vis[i]) {
	// 		printDfs(i);
	// 	}
	// }
	// cout << endl; 
	// for(int i = 1; i < nax; i++) {
	// 	cout << "time_in[" << i << "]: " << time_in[i] << endl; 
	// 	cout << "time_out[" << i << "]: " << time_out[i] << endl; 
	// }

	// for(int i = 1; i < nax; i++) {
	// 	if(!vis[i]) {
	// 		if(cycle_Detection_Undirected_Graph_DFS(i,-1)) {
	// 			cout << "CYCLE PRESENT" << endl; 
	// 			break;
	// 		}
	// 	}
	// }	

	// cout << checkBipartite() << endl; 

	// for(int i = 1; i < nax; i++) {
	// 	if(!vis[i]) {
	// 		if(cycle_Detection_Directed_Graph_DFS(i)){
	// 			cout << "CYCLE" << endl; 
	// 			break;
	// 		}
	// 	}
	// }

	// stack<int> st;
	// for(int i = 0; i < nax; i++) {
	// 	if(!vis[i]) 
	// 		topologicalSort(i, st);
	// }

	// while(!st.empty()) {
	// 	cout << st.top() << " " ; 
	// 	st.pop(); 
	// }

	// topologicalSort_KAHN();	
	
	// vector<int> dist = shortest_Dist_Undirected(0); 
	// for(int i : dist) cout << i << " " ; cout << endl; 

	// vector<int> dist = shortest_Path_DAG().shortest_Paths(0); 
	// for(int i : dist) cout << i << " " ; cout << endl;
	
	// vector<int> distTo = dijkstra(1); 
	// for(int i = 1; i < nax; i++) cout << distTo[i] << " " ; cout << endl; 
	// cout << kruskal(); 

	// bridges(); 

	// articulationPoints(); 
	// kosaraju().scc();

	vector<int> distTo = BellmanFord(0);
	for(int i : distTo) cout << i << " " ; cout << endl; 

}
