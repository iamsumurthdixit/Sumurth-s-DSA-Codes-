#include<bits/stdc++.h>
using namespace std;
#define deb(x) cerr<<"["<<#x<<": "<<x<<"]"<<endl;
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
ll coin_change_1(int n, vector<int> & coins){
	vector<ll> dp(n+1,0);
	dp[0] = 1;

	for(auto c : coins) {
		for(int i = 1; i <= n; i++) {
			if(i - c >= 0)
				dp[i] += dp[i-c];
		}
	}
	return dp[n];
}
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
int rod_cutting(int n, vector<int> & price) {
	int dp[n+1];
	dp[0] = 0;
	for(int i = 1; i <= n; i++){
		int max_val = -1;
		for(int j = 0; j < i; j++) 
			max_val = max(max_val, price[j] + dp[i-j-1]);
		dp[i] = max_val;
	}
	return dp[n];
}
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
	int a[] = {3,5,6,11};
	int n = 4;
	// int sum = 13;
	// cout << subset_sum(n,a,sum);
	// cout << subset_sum_dp(n,a,sum);
	// cout << subset_exits(n,a,sum);
	// cout << count_subsets_of_sum(n,a,sum);
	// cout << minimum_subset_sum(n,a);
	

	return 0;
}
