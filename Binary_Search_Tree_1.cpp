#include<bits/stdc++.h>
using namespace std;

struct node{
	int data;
	node* left;
	node* right;
};
//---------------------------------------------------------------------------------------------------
struct node* add(int x) {
	node* New = new node();
	New->data = x;
	New->left = NULL;
	New->right = NULL;
	return New;
}
//---------------------------------------------------------------------------------------------------
void print_inorder(node* root) {
	if(!root) return ;
	print_inorder(root->left);
	cout << root->data << " ";
	print_inorder(root->right);
}
//---------------------------------------------------------------------------------------------------
struct node* insert(node* Node, int key) {		// O(H) : Height of BST
	if(!Node) 	
		return add(key);  // base case : no root present OR reached the end of the BST
	
	if(key < Node->data) 
		node->left = insert(Node->left, key);
	else if(key > Node->data) 
		node->right = insert(Node->right, key);
	return Node;
}

struct node* search(node* root, int key) {		// O(H) : Height of BST
	if(!root) return NULL;
	if(root->data == key) 
		return root;
	if(root->data > key)
		return search(root->left, key);
	else 
		return search(root->right, key);
}

//---------------------------------------------------------------------------------------------------
void greater_sum_bst(node* root, int & sum) { // done in reverse inorder O(N) 
	if(!root) return ;
	greater_sum_bst(root->right, sum);

	sum += root->data;				// adding the root data
	root-data = sum - root->data;		

	greater_sum_bst(root->left, sum);
}
void greater_sum_bst_caller(node* root) {
	int sum = 0;
	greater_sum_bst(root, sum);
}
//---------------------------------------------------------------------------------------------------
// SORTED ARRAY TO BALANCED BST
/* 
	mid of the sorted array is made root , recursively done for the left half of the array to create
	left sub tree and right half of the array to create right sub tree 
*/
struct node* sorted_array_to_bst(int a[], int strt, int end) {  // O(N) : N = number of nodes  // end = n - 1
	if(strt > end) return NULL;
	int mid = (strt + end) / 2;

	struct node* root = add(a[mid]);
	root->left = sorted_array_to_bst(int a[], 0, mid - 1);
	root->right = sorted_array_to_bst(int a[], mid + 1, end);
	return root;
}		// https://www.geeksforgeeks.org/sorted-array-to-balanced-bst/
//---------------------------------------------------------------------------------------------------
// IMPORTANT : convert binary tree to BST such that it's structure remains same 
// 1 . copy inorder traversal order of the binary tree and store in a array
// 2 . sort the array 
// 3 . again do inorder of the binary tree and update the values from the array
// Time = O(NlogN) 
void inorder(node* root, vector<int> & keys) {
	if(!root) return ;
	inorder(root->left, keys);
	keys.push_back(root->data);
	inorder(root->right, keys);
}
void assign_keys(node* root, vector<int> & keys, int &i) {
	if(!root) return ;
	assign_keys(root->left, keys, i);
	root->data = keys[i];
	i++;
	assign_keys(root->right, keys, i);
}
void binary_tree_to_BST(node* root) {
	vector<int> keys ;
	inorder(root, keys);
	sort(keys.begin(),keys.end());
	int i = 0;
	assign_keys(root, keys, i);
}
//---------------------------------------------------------------------------------------------------
// CREATE BALANCED BST OF MINIMUM HEIGHT FROM ANY BST ; Time : O(N) : N = number of nodes 
struct node* create_balanced_bst(vector<int> & keys, int strt, int end) {
	if(strt > end) 
		return NULL;
	int mid = (strt + end) / 2;
	node* root = add(keys[mid]);
	root->left = create_balanced_bst(keys, 0, mid - 1);
	root->right = create_balanced_bst(keys, mid + 1, end);
	return root;
}
void normal_bst_to_balanced_bst(node* root) {
	vector<int> keys;
	inorder(root, keys);
	node* root = create_balanced_bst(keys, 0, keys.size() - 1);
}
//---------------------------------------------------------------------------------------------------
int minimum_value_in_bst(node* Root) {
	node* root = Root;
	while(root->left)
		root = root->left;
	return root->data;
}
//---------------------------------------------------------------------------------------------------
			// ARRAYS CONCEPT : FIND NGE (NEXT GREATER ELEMENT ON RIGHT IN O(N))
					void nge(vector<int> & v) {
						int n = v.size();

						stack<int> s;
						s.push(v[0]);

						for(int i = 1; i < n; i++) {
							if(s.empty()){
								s.push(v[i]);
								continue;	
							}
							while(s.empty() == false && s.top() < v[i]) {	//	O(1)
								cout << s.top() << " --> " << v[i] << endl;
								s.pop();
							}
							s.push(v[i]);
						}
						while(s.empty() == false) {
							cout << s.top() << " --> " << -1 << endl;
							s.pop();
						}
					}
//---------------------------------------------------------------------------------------------------
bool can_preorder_represent_bst(vector<int> & pre) {
	int n = pre.size();
	stack<int> s;
	int root = INT_MIN;

	for(int i = 0; i < n; i++) {
		if(pre[i] < root) 
			return false;
		while(!s.empty() && s.top() < pre[i]) {
			root = s.top();
			s.pop();
		}
		s.push(pre[i]);
	}
	return true;
}
//---------------------------------------------------------------------------------------------------
// LCA in a BST 
struct node* LCA(node* root, int a, int b) { 	// time = O(H) : Height of BST , space = O(1)
	if(!root)
		return NULL;
	if(root->data > a && root->data > b) // lca lies in the left subtree 
		return LCA(root->left, a , b);

	if(root->data < a && root->data < b) // lca lies in the right subtree
		return LCA(root->right, a, b);

	return root;	// if one is greater than root data and other smaller, then this node is the lca
}
//---------------------------------------------------------------------------------------------------
// INODER SUCCESSOR AND PREDECESSOR IN BST USING ITERATIVE APPROACH 
void Pred_Succ(node* root, node* & pre, node* & suc, int key) {
	if(!root) return ;

	while(root != NULL) {			// searching for given key in BST 
		if(root->data == key) {
		
			if(root->left) { 		// pre is largest key in left subtree
				pre = root->left;
				while(pre->right)
					pre = pre->right;
			}
			if(root->right) {		// suc is smallest key in right subtree
				suc = root->right;
				while(suc->left)
					suc = suc->left;
			}
			return;
		} 
		else if(root->data < key) {		// key present in right subtree and 
			pre = root;					// pre is root and suc is right child of root 
			suc = root->right;
		}
		else {						// key is present in left subtree and suc is root 
			suc = root;				// pre is left child of root 
			pre = root->left;
		}
	}
}
//---------------------------------------------------------------------------------------------------
void Pred_Succ_2(node* root, node* & pre, node* & suc, int key) {
	if(!root) return ;

	if(root->data == key) {
		if(root->left) {
			pre = root->left;
			while(pre->right) 
				pre = pre->right;
		}
		if(root->right) {
			suc = root->right;
			while(suc->left)
				suc = suc->left;
		}
		return;
	}
	if(root->data > key) {
		suc = root;
		Pred_Succ_2(root->left, pre, suc, key);
	}
	if(root->data < key) {
		pre = root;
		Pred_Succ_2(root->right, pre, suc, key);
	}
}
//---------------------------------------------------------------------------------------------------
void closest_to_given_key_util(node* root, int key, int & diff, int & ans) {
	if(!root) return ;
	if(root->data == key) {
		ans = key;
		return;
	}
	if(diff > abs(root->data - key)){
		diff = abs(root->data - key);
		ans = root->data;
	}
	if(key < root->data) 
		return closest_to_given_key_util(root->left, key, diff, ans);
	else 
		return closest_to_given_key_util(root->right, key, diff, ans);
}
int closest_to_given_key(node* root, int key) {
	int diff = 1e5 + 7; int ans = -1;
	closest_to_given_key_util(root, key, diff, ans);
	return ans;
}
//---------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------------

//========== DRIVER FUNCTION ==================================================================================================
int main() {
	
	return 0;
}