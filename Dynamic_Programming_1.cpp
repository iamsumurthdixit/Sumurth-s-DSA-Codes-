#include<bits/stdc++.h>
using namespace std;
#define deb(x) cerr<<"["<<#x<<": "<<x<<"]"<<endl;
#define Fi first
#define Se second
// #define deb2(x) cerr<<"["<<#x.Fi<<": "<<x.Fi<<" : "<<#x.Se<<": "<<x.Se<<"]"<<endl;
#define ll unsigned long long
#define pii pair<int, int> 
#define pll pair<ll,ll>

const int INF = 1e6 + 6;
const ll nax = 1e5 + 5;
const ll MOD = 1e9 + 7;

//=======================================================================================================
ll ugly(ll n) {
	vector<ll> v(n+1,0);
	ll i2 = 0, i3 = 0, i5 = 0;
	v[0] = 1;
	for(int i = 1; i <= n; i++) {
		ll nxt_by_2 = v[i2] * 2;
		ll nxt_by_3 = v[i3] * 3;
		ll nxt_by_5 = v[i5] * 5;
		
		ll nxt_ugly = min(nxt_by_2, min(nxt_by_5, nxt_by_3));
		v[i] = nxt_ugly;

		if(nxt_ugly == nxt_by_2) 
			++i2;
		if(nxt_ugly == nxt_by_3) 
			++i3;
		if(nxt_ugly == nxt_by_5) 
			++i5;
		// deb(v[i])
	}
	return v[n-1];	// first ugly number is 1 
}
//=======================================================================================================
int tiles(int n) {
	vector<int> dp(n+1,0);
	dp[1] = 1;
	dp[2] = 2;
	for(int i = 3; i <= n; i++) {
		dp[i] = dp[i-1] + dp[i-2];
		// deb(dp[i])
	}
	return dp[n];
}
//=======================================================================================================
int gold_mine(vector<vector<int>> & mine, int n, int m) {
	vector<vector<int>> gold(n, vector<int>(m));
	// int gold[6][6];

	for(int i = 0; i < n; i++){
		gold[i][0] = mine[i][0];
	}
	
	for(int col = 1; col < m; col++){
		for(int row = 0; row < n; row++){
			int r = gold[row][col-1];
			int rup = (row - 1) > 0 ?  gold[row-1][col-1] : 0;
			int rdn = (row + 1) < n ? gold[row+1][col-1] : 0;

			gold[row][col] = mine[row][col] + max(r, max(rup,rdn));
		}
	}

	int mx = -1;
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < m; j++) {
			if(gold[i][j] > mx)
				mx = gold[i][j];
			// cout << gold[i][j] << " ";
		} 
		// cout << endl;
	}
	return mx;
}
//=======================================================================================================
ll permutation(int n, int k) {
	// int fact[1000];
	// fact[0] = 1;
	// for(int i = 1; i <= n; i++) {
	// 	fact[i] = fact[i-1] * i;
	// }
	// return fact[n] / fact[n-k];
	ll P = 1;
	for(int i = 0; i < k; i++)
		P *= (n-i);
	return P;
}
//=======================================================================================================
// ll coin_change_1(int n, vector<int> & coins){
// 	vector<ll> dp(n+1,0);
// 	dp[0] = 1;

// 	for(auto c : coins) {
// 		for(int i = 1; i <= n; i++) {
// 			if(i - c >= 0)
// 				dp[i] += dp[i-c];
// 		}
// 	}
// 	return dp[n];
// }
//=======================================================================================================
ll friend_pairing(int n) {
	ll dp[n+1];
	memset(dp,0,sizeof(dp));
	dp[1] = 1;
	dp[2] = 2;
	for(int i = 3; i <= n; i++){
		dp[i] = dp[i-1] + (n-1) * dp[i-2];
	}
	return dp[n];
}
//=======================================================================================================
bool subset_exits(int n, int set[], int sum) {
	if(sum == 0 || n == 0) return true;
	
	bool dp[n+1][sum+1];

	for(int i = 0; i <= n; i++)
		dp[i][0] = true;			// if sum is 0 , then null set can be chosen
	for(int i = 1; i <= sum; i++)
		dp[0][i] = false;			// if sum != 0 and n = 0 then it is fasle

	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= sum; j++){
			if(j < set[i-1])
				dp[i][j] = dp[i-1][j];
			if(j >= set[i-1])
				dp[i][j] = dp[i-1][j] || dp[i-1][j-set[i-1]];
		}
	}
	return dp[n][sum];
}	
//=======================================================================================================
// int rod_cutting(int n, vector<int> & price) {
// 	int dp[n+1];
// 	dp[0] = 0;
// 	for(int i = 1; i <= n; i++){
// 		int max_val = -1;
// 		for(int j = 0; j < i; j++) 
// 			max_val = max(max_val, price[j] + dp[i-j-1]);
// 		dp[i] = max_val;
// 	}
// 	return dp[n];
// }
//=======================================================================================================
/*
		0/1 KNAPSACK
*/
//=======================================================================================================
int knapsack_01_recursion(int wt[], int val[], int W, int n) {
	if(n == 0 || W == 0) 
		return 0;
	if(wt[n-1] > W)
		return knapsack_01_recursion(wt,val,W,n-1);
	return max(val[n-1] + knapsack_01_recursion(wt,val,W-wt[n-1],n-1), 
		knapsack_01_recursion(wt,val,W,n-1));
}
//=======================================================================================================
int knapsack_01(int wt[], int val[], int W, int n) {
	int dp[n+1][W+1];
	memset(dp,0,sizeof(dp));
	
	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= W; j++){
			
			if(wt[i-1] > j)
				dp[i][j] = dp[i-1][j];
			
			else 
				dp[i][j] = max(val[i-1] + dp[i-1][j-wt[i-1]], dp[i-1][j]);	
		}
	}
	return dp[n][W];
}
//=======================================================================================================
// bool subset_sum(int n, int a[], int x, int sum) {
// 	if(x == sum) 
// 		return true;
// 	if(n == 0 && sum != x) 
// 		return false;
// 	if(a[n-1] > sum) 
// 		return subset_sum(n-1,a,x,sum);
// 	return subset_sum(n-1, a, x + a[n-1], sum) || subset_sum(n-1, a, x, sum);
// }
bool subset_sum(int n, int a[], int sum) {
	if(sum == 0) 
		return true;
	if(n == 0 && sum != 0) 
		return false;
	if(a[n-1] > sum) 
		return subset_sum(n-1,a,sum);
	return subset_sum(n-1, a, sum - a[n-1]) || subset_sum(n-1, a, sum);
}
//=======================================================================================================
bool subset_sum_dp(int n, int a[], int sum) {
	bool dp[n+1][sum+1];

	for(int i = 0; i < sum + 1; i++)
		dp[0][i] = true;
	for(int i = 1; i < n + 1; i++)
		dp[i][0] = false;

	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= sum; j++){

			if(a[i-1] > j) 
				dp[i][j] = dp[i-1][j];
			else
				dp[i][j] = dp[i-1][j] || dp[i-1][j-a[i-1]];
		}
	}
	return dp[n][sum];
}
//=======================================================================================================
bool equal_sum_partition(int n, int a[]) {
	int sum = 0;
	for(int i = 0; i < n; i++) 
		sum += a[i];
	if(sum & 1) return false;
	int s = sum / 2;
	
	bool dp[n+1][s+1];

	for(int i = 0; i < s + 1; i++)
		dp[0][i] = true;
	for(int i = 1; i < n + 1; i++)
		dp[i][0] = false;

	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= s; j++) {
			if(j < a[i-1])
				dp[i][j] = dp[i-1][j];
			else 
				dp[i][j] = dp[i-1][j] || dp[i-1][j-a[i-1]];
		}
	}
	return dp[n][s];
}
//=======================================================================================================
int count_subsets_of_sum(int n, int a[], int sum) {
	int dp[n+1][sum+1];
	// memset(dp,0,sizeof(dp));
	for(int i = 0; i <= n; i++)
		dp[i][0] = 1;
	for(int i = 1; i <= sum; i++)
		dp[0][i] = 0;

	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= sum; j++) {
			if(a[i-1] > j)
				dp[i][j] = dp[i-1][j];
			else
				dp[i][j] = dp[i-1][j] + dp[i-1][j-a[i-1]];
		}
	}
	return dp[n][sum];
}
//=======================================================================================================
int minimum_subset_sum(int n, int a[]) {
	int sum = 0;
	for(int i = 0; i < n; i++) 
		sum += a[i];
	int total_sum = sum;
	sum = (sum + 1) / 2;
	bool dp[n+1][sum+1];
	for(int i = 0; i < n + 1; i++)
		dp[i][0] = true;
	for(int i = 1; i < sum + 1; i++)
		dp[0][i] = false;

	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= sum; j++){
			if(a[i-1] > j)
				dp[i][j] = dp[i-1][j];
			else
				dp[i][j] = dp[i-1][j] || dp[i-1][j-a[i-1]];
		}
	}
	int min_diff = INF; 
	for(int j = 0; j <= sum; j++){
		if(dp[n][j]) {
			if(abs(total_sum - 2 * j) < min_diff) {
				min_diff = total_sum - 2 * j;
			}
		}
	}
	return min_diff;
}
//=======================================================================================================
int count_subsets_with_diff(int n, int a[], int diff) {
	int total_sum = 0; 
	for(int i = 0; i < n; i++)
		total_sum += a[i];						// s1 + s2 = T; s1 - s2 = diff; s1 = diff + s2;
	int sum = (diff + total_sum) / 2;			// s1 = diff + T - s1; s1 = (T + diff) / 2

	int dp[n+1][sum+1];
	for(int i = 0; i <= n; i++)
		dp[i][0] = 1;
	for(int i = 1; i <= sum; i++)
		dp[0][i] = 0;

	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= sum; j++){
			if(a[i-1] > j) 
				dp[i][j] = dp[i-1][j];
			else
				dp[i][j] = dp[i-1][j] + dp[i-1][j-a[i-1]];
		}
	}
	return dp[n][sum];
}
//=======================================================================================================
int target_sum(int n, int a[], int target) {
	// entire code is -> int count_subsets_with_diff(int n, int a[], int diff) 
}
//=======================================================================================================
/*
		UNBOUNDED KNAPSACK
*/
//=======================================================================================================
int unbounded_knapsack_recursion(int wt[], int val[], int W, int n){
	if(n == 0) {
		if(W == 0) 
			return 1;
		else 
			return 0;
	}
	if(W == 0)
		return 1;
	if(wt[n-1] > W) 
		return unbounded_knapsack_recursion(wt, val, W, n-1);
	else 
		return max(val[n-1] + unbounded_knapsack_recursion(wt, val, W - wt[n-1], n), 
			unbounded_knapsack_recursion(wt, val, W, n-1));
}
//=======================================================================================================
int unbounded_knapsack(int wt[], int val[], int W, int n) {
	int dp[n+1][W+1];
	for(int i = 0; i < n + 1; i++) 
		dp[i][0] = 1;
	for(int i = 1; i < W + 1; i++)
		dp[0][i] = 0;

	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= W; j++) {
			if(wt[i-1] > j) // if we are not taking the element it is processed , we will not take it again
				dp[i][j] = dp[i-1][j];
			else  			// here , if we are taking the element once, then we need to call it again for same i
				dp[i][j] = max(val[i-1] + dp[i][j-wt[i-1]], dp[i-1][j]);
		}
	}
	return dp[n][W];
}
//=======================================================================================================
int rod_cutting(int L, int price[]) {
	int len[L];
	for(int i = 0; i < L; i++)
		len[i] = i + 1;

	int dp[L+1][L+1];	// first dimension for length and second for profit
	for(int i = 0; i < L + 1; i++) {
		dp[i][0] = 0;	// length is 0;
		dp[0][i] = 0;	// profit is 0;
	}

	for(int i = 1; i <= L; i++) {
		for(int j = 1; j <= L; j++) {
			if(len[i-1] > j) 
				dp[i][j] = dp[i-1][j];
			else 
				dp[i][j] = max(price[i-1] + dp[i][j-len[i-1]], dp[i-1][j]);
		}
	}
	return dp[L][L];
}
//=======================================================================================================
int coin_change_1(int n, int cns[], int sum) {		// total unique ways to get sum 
	int dp[n+1][sum+1];
	for(int i = 0; i <= n; i++)
		dp[i][0] = 1;
	for(int i = 1; i <= sum; i++) 
		dp[0][i] = 0;

	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= sum; j++) {
			if(cns[i-1] > j)
				dp[i][j] = dp[i-1][j];
			else
				dp[i][j] = dp[i][j-cns[i-1]] + dp[i-1][j];
		}
	}
	return dp[n][sum];
}
//=======================================================================================================
int coin_change_2(int n, int cns[], int sum) {	// minimum number of coins needed to sum 
	int dp[n+1][sum+1];
	for(int i = 0; i <= n; i++)
		dp[i][0] = 0;
	for(int i = 0; i <= sum; i++)
		dp[0][i] = INF;

	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= sum; j++){
			if(cns[i-1] > j) 
				dp[i][j] = dp[i-1][j];
			else
				dp[i][j] = min(1 + dp[i][j-cns[i-1]], dp[i-1][j]);
		}
	}
	return dp[n][sum];
}
//=======================================================================================================
/*
	LONGEST COMMON SUBSEQUENCE

	Note : for LCS , two strings are I/O .. if 1 is given, then second is some function of first string 
	output is generally integer(can be int arrays like in LIS)
*/
//=======================================================================================================
int longest_common_subsequence_recursion(int n, string a, int m, string b) {
	if(n == 0 || m == 0)
		return 0;
	if(a[n-1] == b[n-1]) 
		return 1 + longest_common_subsequence_recursion(n-1, a, m-1, b);
	else 
		return max(longest_common_subsequence_recursion(n-1, a, m, b), 
				longest_common_subsequence_recursion(n, a, m-1, b));
}
//=======================================================================================================
int longest_common_subsequence(int n, string a, int m, string b) {
	int dp[n+1][m+1];
	for(int i = 0; i < n+1; i++)
		dp[i][0] = 0;
	for(int i = 0; i < m+1; i++)
		dp[0][i] = 0;

	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= m; j++){
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else 
				dp[i][j] = max(dp[i-1][j] , dp[i][j-1]);
		}
	}
	return dp[n][m];
}
//=======================================================================================================
int longest_common_substring(int n, string a, int m, string b) {
	if(n == 0 || m == 0) 
		return 0;
	int dp[n+1][m+1];
	for(int i = 0; i < n + 1; i++)
		dp[i][0] = 0;
	for(int j = 0; j < m + 1; j++)
		dp[0][j] = 0;

	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= m; j++) {
			if(a[i-1] == b[j-1])	// if it is equal then add 1 and continue 
				dp[i][j] = 1 + dp[i-1][j-1];
			else					// otherwise, flow is broken, start with 0 
				dp[i][j] = 0;
		}
	}
	return dp[n][m];
}
//=======================================================================================================
string LCS(int n, string a, int m, string b) {
	int dp[n+1][m+1];
	for(int i = 0; i < n + 1; i++) 
		dp[i][0] = 0;
	for(int i = 0; i < m + 1; i++)
		dp[0][i] = 0;

	for(int i = 1; i < n + 1; i++) {
		for(int j = 1; j < m + 1; j++) {
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else 
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	int i = n ; int j = m; 
	string lcs = "";
	while(i > 0 && j > 0) {
		if(a[i-1] == b[j-1]) {
			lcs += a[i-1];
			i--;
			j--;
		} else {		// move in the direction of maximum (where change happened first)
			if(dp[i-1][j] > dp[i][j-1])
				i--;
			else 
				j--;
		}
	}
	reverse(lcs.begin(), lcs.end());
	return lcs;
}
//=======================================================================================================
int shortest_common_supersequence(int n, string a, int m, string b) {
	int dp[n+1][m+1];
	for(int i = 0; i < n + 1; i++) 
		dp[i][0] = 0;
	for(int i = 0; i < m + 1; i++)
		dp[0][i] = 0;

	for(int i = 1; i < n + 1; i++) {
		for(int j = 1; j < m + 1; j++) {
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else 
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}					// LCS is common , so add it once in supersequence(subtract once in total)
	return ((n + m) - dp[n][m]);
}
//=======================================================================================================
pair<int,int> minimum_insrtn_deltn_a_to_b(int n, string a, int m, string b) {
	int dp[n+1][m+1];
	for(int i = 0; i < n + 1; i++) 
		dp[i][0] = 0;
	for(int i = 0; i < m + 1; i++)
		dp[0][i] = 0;

	for(int i = 1; i < n + 1; i++) {
		for(int j = 1; j < m + 1; j++) {
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else 
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	int deltn = n - dp[n][m];
	int insrtn = m - dp[n][m];
	return make_pair(insrtn,deltn);
}
//=======================================================================================================
/*
	LONGEST PALINDROMIC SUBSEQUENCE

	LPS of string a = LCS(a, reverse(a))
*/
//=======================================================================================================
int longest_palindromic_subsequence(int n, string a) {
	string b = a;
	reverse(b.begin(), b.end());

	int dp[n+1][n+1];
	for(int i = 0; i < n + 1; i++) {
		dp[i][0] = 0;
		dp[0][i] = 0;
	}

	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= n; j++){
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	return dp[n][n];
}
//=======================================================================================================
int longest_increasing_subsequence(int n, int a[]) {
	int b[n];						// LCS : second array = sort(a)
	for(int i = 0; i < n; i++) 
		b[i] = a[i];
	
	sort(b, b+n);
	int dp[n+1][n+1];
	for(int i = 0; i <= n; i++) {
		dp[i][0] = 0;
		dp[0][i] = 0;
	}
	for(int i = 1; i < n + 1; i++){
		for(int j = 1; j < n + 1; j++) {
			if(a[i-1] == b[j-1]) 
				dp[i][j] = 1 + dp[i-1][j-1];
			else
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	return dp[n][n];
}
//=======================================================================================================
int longest_palindromic_substring(int n, string a) {
	// string b = a;
	// reverse(b.begin(), b.end());

	// int dp[n+1][n+1];
	// for(int i = 0; i < n+1 ;i++) {
	// 	dp[i][0] = 0;
	// 	dp[0][i] = 0;
	// }
	// for(int i = 1; i <= n; i++) {
	// 	for(int j = 1; j <= n; j++) {
	// 		if(a[i-1] == b[j-1]) 
	// 			dp[i][j] = 1 + dp[i-1][j-1];
	// 		else
	// 			dp[i][j] = 0;
	// 	}
	// }
	// return dp[n][n];
}
//=======================================================================================================
int minimum_deltn_to_palndrm(int n, string a) {
	string b = a;
	reverse(b.begin(), b.end());
	int dp[n+1][n+1];
	for(int i = 0; i < n + 1; i++) {
		dp[i][0] = 0;
		dp[0][i] = 0;
	}

	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	return n - dp[n][n];
}
//=======================================================================================================
string print_shortest_common_supersequence(int n, string a, int m, string b){
	int dp[n+1][m+1];
	for(int i = 0; i <= n; i++) 
		dp[i][0] = 0;
	for(int j = 0; j <= m; j++) 
		dp[0][j] = 0;

	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else
				dp[i][j] = max(dp[i-1][j] , dp[i][j-1]);
		}
	}
	int i = n; int j = m;
	string scs = "";
	while(i > 0 && j > 0) {
		if(a[i-1] == b[j-1]) {
			scs += a[i-1];
			i--;
			j--;
		} else {
			if(dp[i-1][j] > dp[i][j-1]) {
				scs += a[i-1];
				i--;
			} else {
				scs += b[j-1];
				j--;
			}
		}
	}
	while(i > 0) {
		scs += a[i-1];
		i--;
	}
	while(j > 0) {
		scs += b[j-1];
		j--;
	}
	reverse(scs.begin(), scs.end());
	assert(scs.length() == shortest_common_supersequence(n,a,m,b));
	return scs;
}
//=======================================================================================================
int longest_repeating_subsequence(int n, string a){
	string b = a;
	int dp[n+1][n+1];
	for(int i = 0; i <= n; i++) {
		dp[i][0] = 0;
		dp[0][i] = 0;
	}
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			if(a[i-1] == b[j-1] && i != j) // i != j to ensure there are more than 1 occurences of a char. 
				dp[i][j] = 1 + dp[i-1][j-1];
			else
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	return dp[n][n];
}
//=======================================================================================================
bool sequence_pattern_searching(int n, string a, int m, string b) {
	int dp[n+1][m+1];
	for(int i = 0; i < n + 1; i++) 
		dp[i][0] = 0;
	for(int j = 0; j < m + 1; j++)
		dp[0][j] = 0;
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= m; j++) {
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	return (dp[n][m] == n);
}
//=======================================================================================================
int minimum_insrtn_to_palndrm(int n, string a) {
	string b = a;
	reverse(b.begin(), b.end());
	int dp[n+1][n+1];
	for(int i = 0; i < n + 1; i++ ){
		dp[i][0] = 0;
		dp[0][i] = 0;
	}
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else
				dp[i][j] = max(dp[i][j-1], dp[i-1][j]);
		}
	}
	return n - dp[n][n];
}
//=======================================================================================================
/*
	MATRIX CHAIN MULTIPLICATION
*/

// 		FORMAT AND IDENTIFICATION
// given a array or string

// int solve(int a[], int i, int j) {
// 	if(i > j) 
// 		return 0;
// 	for(int k = i; k < j; k++) {
// 		int temp_ans = solve(a,i,k) + solve(a,k+1,j);

// 		function for ans from temp_ans;; 
// 		ans = function(temp_ans);
// 	}
// 	return ans;
// }

// for n size array , we have n-1 matrices, if i >= 1; dimension of i-th matrix is a[i-1] X a[i];
//=======================================================================================================
int solve(int a[], int i, int j){ 	// first call , i == 1 and j == n-1
	if(i >= j) //i == j means only 1 element, for even 1 matrix we need two dimensions 
		return 0;
	int ans = INF;
	for(int k = i; k <= j-1; k++) {
		int temp_ans = solve(a,i,k) + solve(a,k+1,j) + (a[i-1] * a[k] * a[j]);
		ans = min(ans, temp_ans);
	}
	return ans;
}
//=======================================================================================================

//=======================================================================================================

//=======================================================================================================

//========= DRIVER FUNCTION ===================================================================================================
int main() {
	// cout << ugly(4);
	// cout << tiles(4);

	// int n = 4; int m = 4;
	// vector<vector<int>> mine(n, vector<int>(m));

	// for(int i = 0; i < n; i++) 
	// 	for(int j = 0; j < m; j++) 
	// 		cin >> mine[i][j];


	// for(int i = 0; i < n; i++){
	// 	for(int j = 0; j < m; j++)
	// 		cout << mine[i][j] << " ";
	// 	cout << endl;
	// }

	// cout << gold_mine(mine,n,m);
	// cout << permutation(43,12);

	// int N = 4; vector<int> coins = {1,2,3};
	// cout << coin_change_1(N,coins);
	// cout << friend_pairing(4);

	// int val[] = {60,100,120};
	// int wt[] = {10,20,30};
	// int W = 50;
	// int n = 3;
	// deb(knapsack_01(wt,val,W,n))
	// deb(knapsack_01_recursion(wt,val,W,n))
	// int a[] = {1,1,2,3};
	// int n = 4;
	// int diff = 1;
	// int sum = 13;
	// cout << subset_sum(n,a,sum);
	// cout << subset_sum_dp(n,a,sum);
	// cout << subset_exits(n,a,sum);
	// cout << count_subsets_of_sum(n,a,sum);
	// cout << minimum_subset_sum(n,a);
	// cout << count_subsets_with_diff(n,a,diff);

	// int n = 8;
	// int a[] = {1,5,8,9,10,17,17,20};
	// cout << rod_cutting(n,a);
	
	// int n = 3;
	// int cns[] = {1,7,5};
	// int sum = 11;
	// cout << coin_change_2(n,cns,sum);

	// string a = "dixit";
	// int n = a.length();
	// string b = "adgxt";
	// int m = b.length();

	// deb(longest_common_subsequence_recursion(n,a,m,b))
	// deb(longest_common_subsequence(n,a,m,b))	
	// deb(longest_common_substring(n,a,m,b))
	// deb(LCS(n, a, m, b))
	// deb(shortest_common_supersequence(n,a,m,b))
	// deb(print_shortest_common_supersequence(n,a,m,b))
	// pair<int,int> x = minimum_insrtn_deltn_a_to_b(n,a,m,b);
	// deb(x.Fi)
	// deb(x.Se)
	
	// string a = "abgcba";
	// int n = a.length();
	// deb(longest_palindromic_subsequence(n,a))
	// deb(longest_palindromic_substring(n,a))
	// deb(minimum_deltn_to_palndrm(n,a))
	// deb(minimum_insrtn_to_palndrm(n,a))

	// int a[] = {6,2,5,1,7,4,8,3};
	// int n = sizeof(a) / sizeof(int);

	// deb(longest_increasing_subsequence(n,a))
	// string a = "AABEBCDD";
	// int n = a.length();
	// deb(longest_repeating_subsequence(n,a))

	// string a = "AXY";
	// int n = a.length();
	// string b = "ADXCPY";
	// int m = b.length();
	// deb(sequence_pattern_searching(n,a,m,b))




	
	return 0;
}
