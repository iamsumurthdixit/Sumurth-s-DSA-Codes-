#include<bits/stdc++.h>
using namespace std;

struct node{
    int data;
    node* left; 
    node* right;
};
//===============================================================================================
struct node* add(int data){
	node* root = new node();
	root->data = data;
	root->left = NULL;
	root->right = NULL;
	return root;
}
//===============================================================================================
int height(node* root){
    if(!root) return 0;
    int l = height(root->left);
    int r = height(root->right);
    if(l > r) return 1+l;
    else return 1+r;
}
// void reverse_level_order(node* root){
//     int h = height(root);
//     for(int i=h;i>=1;i--)
//         print_given_level(root,i);
// }
// void print_given_level(node* root, int level){
//     if(!root) return;
//     if(level == 1){
//         cout<<root->data<<" ";
//     }
//     print_given_level(root->left, level-1);
//     print_given_level(root->right, level-1);
// }
//===============================================================================================
void reverse_level_order(node* root){
	if(!root) return;
	stack<node*> S;
	queue<node*> Q;
	Q.push(root);
	while(!Q.empty()){
		node* temp = Q.front();
		Q.pop();
		S.push(temp);

		if(temp->right)
			Q.push(temp->right);
		if(temp->left) 
			Q.push(temp->left);
	}
	while(!S.empty()){
		node* temp = S.top();
		cout<<temp->data<<" ";
		S.pop();
	}
}
//===============================================================================================
void spiral_level_order(node* root){
	if(!root) return ;
	stack<node*> s1;
	stack<node*> s2;

	s1.push(root);

	while(!s1.empty() || !s2.empty()){
		while(!s1.empty()){
			node* temp = s1.top();
			s1.pop();
			cout<<temp->data<<" ";
			if(temp->right)
				s2.push(temp->right);
			if(temp->left)
				s2.push(temp->left);
		}
		while(!s2.empty()){
			node* temp = s2.top();
			s2.pop();
			cout<<temp->data<<" ";
			if(temp->left)
				s1.push(temp->left);
			if(temp->right)
				s1.push(temp->right);
		}
	}
}
//===============================================================================================
// static int cnt = 0;
// int cnt_nodes(node* root){
// 	if(!root) return 0;
// 	cnt_nodes(root->left);
// 	cnt_nodes(root->right);
// 	cnt++;
// 	return cnt;
// }
//===============================================================================================

// 	ANTI CLOCK BOUNDRY TRAVERSAL OF A TREE

void print_boundry_left(node* root){
	//  Done in TOP DOWN manner
	if(!root) return;

	if(root->left){
		cout<<root->data<<" ";
		print_boundry_left(root->left);
	}
	else if(root->right){
		cout<<root->data<<" ";
		print_boundry_left(root->right);
	}
}
void print_boundry_right(node* root){
	// Done in BOTTOM UP manner coz going up
	if(!root) return ;
	
	if(root->right){
		print_boundry_right(root->right);
		cout<<root->data<<" ";
	} 
	else if(root->left){
		print_boundry_right(root-left);
		cout<<root->data<<" ";
	}
}
void print_leaves(node* root){
	// INORDER like 
	if(!root) return;
	print_leaves(root->left);
	if(!root->left && !root->right)
		cout<<root->data<<" ";
	print_leaves(root->right);
}
void print_boundry(node* root){
	// driver function 
	if(!root) return ;

	cout<<root->data<<" ";
	print_boundry_left(root->left);

	print_leaves(root->left);
	print_leaves(root->right);

	print_boundry_right(root->right);
}
//===============================================================================================
bool sum_parent(node* root){
	if(!root) return true;
	if(!root->left && !root->right)	
		return true;
	if(root->left && root->right){
		if(root->data == root->left->data + root->right->data)
			return (sum_parent(root->left) && sum_parent(root->right) );
		else return false;
	}
	if(root->left && !root->right){
		if(root->data == root->left->data)
			return sum_parent(root->left);
		else
			return false;
	}
	if(root->right && !root->left){
		if(root->data == root->right->data)
			return sum_parent(root->right);
		else 
			return false;
	}
	return true;
}
//===============================================================================================
int sum(node* root){
	if(!root) return 0;
	return root->data + sum(root->left) + sum(root->right);	
}
bool isSumTree(node* node){
	int ls , rs;
	if(node == NULL || node->left == NULL && node->right == NULL)
		return true;
	ls = sum(root->left);
	rs = sum(root->right);

	if((node->data == ls + rs) && isSumTree(node->left) && isSumTree(node->right))
		return true;
	return false;			
}
//===============================================================================================
bool sum_covered_uncovered(node* root){
	if(!root) return true;
	int sum_co = 0; 
	int sum_un = 0;

}
//===============================================================================================
//  Case 1: searching data 
int get_level_utility(node* root, int data, int level){
	if(!root) return 0;
	if(root->data == data) 
		return level;
	int downlevel = get_level_utility(root->left,data,level+1);
	if(downlevel != 0)
		return downlevel;
	downlevel = get_level_utility(root->right,data,level+1);
		return downlevel;	
}
int get_level(node* root, int data){
	int level = 1;
	return get_level_utility(node* root,int data,level);
}
// Case 2: searching node
int level(node* root, node* ptr, int lev){
	if(!root) return 0;
	if(root == ptr) return lev;
	int l = level(root->left,ptr,lev+1);
	if(l != 0)
		return l;
	return level(root->right,ptr,lev+1);
}
//===============================================================================================
bool check_siblings(node* root, node* a, node* b){
	if(!root) return false;
	return ((root->left == a && root->right == b)||
		(root->right == a && root->left === b)|| 
		check_siblings(root->left,a,b)||
		check_siblings(root->right,a,b));
}
bool check_cousins(node* root, node* a, node* b){
	if(level(root,a,1) == level(root,b,1) && !check_siblings(root,a,b))
			return true;
		else 
			return false;
}
//===============================================================================================
bool check_leaves_same_level_utility(node* root, int level, int* first_leaf_level){
	if(!root) return true;
	
	if(root->left == NULL && root->right == NULL){
		if(*first_leaf_level == 0){
			*first_leaf_level = level;
			return true;
		}
		return (level == *first_leaf_level); // only checking for leaf nodes 
	}
	return (check_leaves_same_level_utility(root->left,level+1,*first_leaf_level) 
			&& check_leaves_same_level_utility(root->right,level+1,*first_leaf_level));
}
bool check_leaves_same_level(node* root){
	if(!root) return true;
	int level = 0;
	int first_leaf_level = 0;
	return check_leaves_same_level_utility(root,level,&first_leaf_level);
}
//===============================================================================================
int total_sum(node* root){
	if(root) 
		return root->data + sum(root->left) + sum(root->right);
	else 
		return 0;
}
int uncovered_sum_left(node* root){
	if(!root) return 0;
	if(!root->left && !root->right)
		return root->data;
	if(root->left) 
		return root->data + uncovered_sum_left(root->left);
	else if(root->right)
		return root->data + uncovered_sum_left(root->right);
}
int uncovered_sum_right(node* root){
	if(!root) return 0;
	if(!root->left && !root->right) 
		return root->data;
	if(root->right) 
		return root->data + uncovered_sum_right(root->right);
	else if(root->left) 
		return root->data + uncovered_sum_right(root->left);
}
int total_uncovered_sum(node* root){
	return root->data + uncovered_sum_left(root->left) + uncovered_sum_right(root->right);
}
//===============================================================================================
void traversal(node* root, vector<node*> &v){
	if(!root) return;
	if(root->left == NULL && root->right == NULL)
		v.push_back(root);
	traversal(root->left,v);
	traversal(root->right,v);
}
bool same_leaf_traversal(node* root1, node* root2){
	vector<node*> v1, v2;
	traversal(root1,v1);
	traversal(root2,v2);
	 // check to kr hi loge ab 
}
//===============================================================================================
int depth(node* root){
	int d = 0;
	while(node  != NULL){
		node = node->left;
		++d;
	}
	return d;
}
bool isPerfect(node* root, int d, int lev){
	if(root == NULL) return true;
	if(root->left == NULL && root->right == NULL){
		return (d == lev + 1);
	}
	if(root->left == NULL || root->right == NULL)
		return false;
	return isPerfect(root->left,d,lev+1) &&
		isPerfect(root->right,d,lev+1);
}
bool isPerfectTree(node* root){
	int d = depth(root);
	return isPerfect(root,d,0);
}
//===============================================================================================
bool mirror_tree(node* root1, node* root2){
	if(!root1 && !root2)
		return true;
	if(!root1 && root2 || root1 && !root2)
		return false;
	if(root1->data == root2->data){
		return mirror_tree(root1->left,root2->right) &&
			mirror_tree(root1->right,root2->left);
	}
	else return false;
}
//===============================================================================================
void find(queue<node*> &q, node* root, int data){
	if(!root) 
		return;
	if(root->data == data)
		q.push(root);
	find(q,root->left,data);
	find(q,root->right,data);
}
bool check_subtree_main(node* root1, node* root2){
	if(root1 == NULL && root2 == NULL)
		return true;
	if(root1 == NULL || root2 == NULL)
		return false;
	if(root1->data != root2->data)
		return false;
	return check_subtree_main(root1->left,root2->left) 
		&& check_subtree_main(root1->right,root2->right);
}
bool check_subtree(node* root1, node* root2){	
	if(!root1 && !root2)
		return true;
	if(!root1 || !root2)
		return false;
	queue<node*> q;
	find(q,root2,root1->data);
	if(q.empty()){
		return false;
	}
	bool yes = false;
	while(!q.empty() && !yes){
		node* temp = q.front();
		q.pop();
		yes = check_subtree_main(root1,temp);
	}
	return yes;
}
//===============================================================================================
// Optimized version of above code to find subtree
void random_code(node* root, int n, int a[], int indx){
	if(indx == n || !root)
		return false;
	if(root->left == NULL && root->right == NULL){
		return ((root->data == a[indx]) && (indx == n-1));
	}
	if(indx < n-1){
		return (root->data == a[indx] && (random_code(root->left,n,a,indx+1)
			|| random_code(root->right,n,a,indx+1)));
	}
}
//===============================================================================================
int height(node* a){
	if(!a) return 0;
	return 1 + max(height(a->left),height(a->right));
}
int diameter(node* root){  // O(N^2) : N = nodes
	if(!root) return 0;

	int l_height = height(root->left);
	int r_height = height(root->right);

	int l_diameter = diameter(root->left);
	int r_diameter = diameter(root->right);

	return max(l_height + r_height + 1, max(l_diameter,r_diameter));	
}
//===============================================================================================
int Optim_diameter(node* root, int* height){ // O(N)
	int l_h = 0, r_h = 0;
	int l_dia = 0, r_dia =0 ;

	if(root == NULL){
		*height = 0;
		return 0;
	}
	l_dia = Optim_diameter(root->left, &l_h);
	r_dia = Optim_diameter(root->right, &r_h);
	
	*height = 1 + max(l_h, r_h);
	return max(l_h + r_h + 1, max(l_dia, r_dia));
}
//===============================================================================================
int Height(node* root, int &ans){
	if(!root) return 0;
	int lh = Height(root->left, ans);
	int rh = Height(root->right, ans);

	ans = max(ans, 1 + lh + rh);
	return 1 + max(lh,rh);
}
int diameter_best(node* root){
	if(!root) return 0;
	int ans = -1;
	int ht_tree = Height(root,ans);
	return ans;
}
//===============================================================================================
void print_path(int path[] , int pathlen){
	for(int i=0 ; i<pathlen; i++)
		cout<<path[i]<<" ";
	cout<<endl;
}
void print_utility(node* root, int path[] ,int pathlen){
	if(!root) return ;
	path[pathlen] = root->data;
	pathlen++;

	if(root->left == NULL && root->right == NULL){
		print_path(path,pathlen);
	} else {
		print_utility(root->left, path, pathlen);
		print_utility(root->right, path, pathlen);
	}
}
void root_to_leaf_paths(node* root){
	if(!root) return ;
	int path[100]; 
	int pathlen = 0;
	print_utility(root,path,pathlen);
}
//===============================================================================================
void print_middle_level_utility(node* slw, node* fst){
	if(slw==NULL || fst==NULL) 
		return ;
	if(fst->left == NULL && fst->right == NULL) {	// fast pointer reached end print slow pointer
		cout << slw->data << " ";
		return ;
	}
	if(fst->left->left) {
		print_middle_level_utility(slw->left, fst->left->left);
		print_middle_level_utility(slw->right, fst->left->left);
	} else {
		print_middle_level_utility(slw->left, fst->left);
		print_middle_level_utility(slw->right, fst->left);
	}

}
void print_middle_level(node* root){
	print_middle_level_utility(root,root);
}
//===============================================================================================
bool isLeaf(node* x) {
	if(!x) return false;
	if(x->left == NULL && x->right == NULL) 
		return true;
	return false;
}
int left_leaf_sum(node* root, int & sum) { 
	if(!root)
		return 0;
	else {
		if(isLeaf(node->left))
			sum += node->left->data;
		else 
			left_leaf_sum(root->left, sum);
		left_leaf_sum(root->right, sum);
	}
	return sum;
}
//===============================================================================================
// printing path to a node from root considering all noded keys to be distinct
bool hasPath(node* root, vector<int> & v, int x) {
	if(!root)
		return false;

	v.push_back(root->data) ;	// just pushing the root data

	if(root->data == x) 		// found x so let it be present and return true;
		return true;
	if(hasPath(root->left,v,x) || hasPath(root->right,v,x))
		return true;			// x is present in left or right subtree 
								// so this node forms a path , let it be present
	v.pop_back();				// else this node is not in the path to x remove it and return false
	return false;			
}
void print_path_to_node(node* root, int x) {
	 vector<int> v;
	 if(hasPath(root,v,x)) {
	 	for(auto i : v) 
	 		cout << i << " -> ";
	 } else {
	 	cout << "No path" << endl;
	 }
}
//===============================================================================================
			// LCA method 1 : time O(N) ; N = numebr of nodes 
// shortened code from above code , returns v and true is present else false
bool find_path(node* root, vector<int> & v, int x) {
	if(!root) return false;
	v.push_back(root->data);
	if(root->data == x) 
		return true;
	if(find_path(root->left, v, x) || find_path(root->right, v, x))
		return true;
	v.pop_back(root->data);
	return false;
}
int findLCA(node* root, int a, int b) {
	vector<int> path1, path2;
	if(!find_path(root, path1, a) || find_path(root, path2, b)) 
		return -1; // either element not present in the tree
	int i;
	for(i = 0; i < path1.size() && i < path2.size(); i++) {
		if(path1[i] != path2[i]) 
			break;
	}
	return path1[i-1]; // last common element in the path or LCA 
}
//===============================================================================================
		// LCA method 2 : time O(N) ; N = numebr of nodes
// finding LCA in single traversal without extra space assuming both the keys are presnt 
// return node* to the LCA
struct node* findLCA2(node* root, int a, int b) {
	if(!root) return NULL;
	
	if(root->data == a || root->data == b) 
		return root;		// if any node matches root data then this is LCA
	
	// checking in left and right subtrees
	node* left_LCA = findLCA2(root->left, a , b);		
	node* right_LCA = findLCA2(root->right, a, b);

	// if one node is in left subtree and other is in right subtree and none of 
	// the two calls are NULL then this root is the LCA
	if(left_LCA && right_LCA)
		return root;

	// if both are NOT in left subtree , check in right subtree and vice versa
	return left_LCA == NULL ? right_LCA : left_LCA ;
}
//===============================================================================================
int find_Level(node* root, int k, int & level) {
	if(!root) return -1;
	if(root->data == k)
		return level;
	int lh = find_Level(root->left, k , 1 + level); // finding in left subtree
	if(lh == -1)			// if not in left subtree
		return find_Level(root->right, k , 1 + level);
	return lh; // returns this when lh is not -1
}
int find_distance_btwn_2_nodes(node* root, int a, int b) { // sum of distances from lca to both nodes 
	node* lca = findLCA2(root, a, b);		// assuming both keys are present 
	int d1 = find_Level(lca, a, 0);
	int d2 = find_Level(lca, b, 0);
	return d1 + d2;
}
//===============================================================================================
bool printAncestors(node* root, int target) { // this function can be used to print the path to a node without 
	if(!root) return false;						// storing the intermediate nodes in a vector
	if(root->data == target) {
		cout << target <<" ";
		return true;
	}
	if(printAncestors(root->left, target) || printAncestors(root->right, target)) {
		cout << root->key << " ";
		return true;
	}
	return false;
}
bool print_common_path_to_2_nodes(node* root, int a, int b) { // common ancestors
	node* lca = findLCA2(root, a, b);
	if(lca == NULL) return false;
	printAncestors(root, lca->data); // important function
}
//===============================================================================================

//===============================================================================================

//===============================================================================================

//===============================================================================================

//===============================================================================================

//===============================================================================================

//===============================================================================================

//===============================================================================================

//===============================================================================================

//============ DRIVER FUNCTION =====================================================================================================================
int main(){
	
	node* root = add(12);
	root->left = add(5);
	root->right = add(7);
	root->left->left = add(3);
	root->right->right = add(1);
	// cout<<sum_parent(root);
	cout<<check_leaves_same_level(root);
	// assert(check_leaves_same_level(root) == 1);

	return 0;
}
