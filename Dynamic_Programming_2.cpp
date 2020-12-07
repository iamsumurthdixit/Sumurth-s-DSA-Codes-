#include<bits/stdc++.h>
using namespace std;
#define ll unsigned long long
#define deb(x) cerr<<"["<<#x<<": "<<x<<"]"<<endl;
#define mmst(x,val) memset(x,val,sizeof(x));

const ll INF = 1e7;
const ll MOD = 1e9 + 7;
const ll nax = 1e5 + 5;

//====================================================================================
int coin_chnge_1(int s[], int n, int sum) {
	int dp[n+1][sum+1];
	for(int i = 0; i < n + 1; i++)
		dp[i][0] = 1;
	for(int i = 1; i < sum + 1; i++)
		dp[0][i] = 0;
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= sum; j++) {
			if(s[i-1] > j)
				dp[i][j] = dp[i-1][j];
			else
				dp[i][j] = dp[i][j-s[i-1]] + dp[i-1][j];
		}
	}
	return dp[n][sum];
}
//====================================================================================
int lrgst_subset(int a[], int n) {
	int dp[n];
	mmst(dp,0);
	for(int i = 1; i < n ;i++) {
		for(int j = 0; j < i; j++) {
			if(a[i] % a[j] == 0)
				dp[i] = max(dp[i], 1 + dp[j]);
		}
	}
	return *max_element(dp, dp+n);
}
//====================================================================================
string longest_common_subsequence(int n, string a, int m, string b) {
	int dp[n+1][m+1];
	for(int i= 0 ; i< n + 1; i++)
		dp[i][0] = 0;
	for(int i = 0; i < m + 1; i++)
		dp[0][i] = 0;
	for(int i= 1; i <= n; i++) {
		for(int j = 1; j <= m; j++) {
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	int i = n; int j = m;
	string lcs = "";
	while(i > 0 && j > 0) {
		if(a[i-1] == b[j-1]) {
			lcs += a[i-1];
			i--;
			j--;
		}
		else if(a[i-1] > b[j-1]) 
			i--;
		else 
			j--;
	}
	reverse(lcs.begin(), lcs.end());
	return lcs;
	// this printing also works for longest repeated subsequence 
}
//====================================================================================
int longest_increasing_subsequence_1(int n, int a[]){
	int dp[n];
	for(int i = 0; i < n; i++) 
		dp[i] = 1;
	for(int i = 1; i < n; i++) {
		for(int j = 0; j < i; j++) {
			if(a[i] > a[j]) {
				dp[i] = max(dp[i] , 1 + dp[j]);
			}
		}
	}
	return *max_element(dp , dp + n);
}
//====================================================================================
int longest_increasing_subsequence_optimized(int n, int a[]) {
	// always choose the minimum element and put in the vector
	vector<int> v;
	v.push_back(a[0]);
	for(int i = 1; i < n; i++) {
		if(v.back() < a[i])
			v.push_back(a[i]);
		else {
			// int id = lower_bound(v.begin(), v.end(), a[i]) - v.begin();
			// v[id] = a[i];
			int l = 0; 
			int r = v.size();
			while(l < r) {
				int m = (l + r) / 2;
				if(v[m] > a[i]) 
					r = m ;
				else
					l = m  + 1;
				// deb(m)
			}
			v[l-1] = a[i];
		}
	}
	return v.size();
}
//====================================================================================
int longest_common_subsequence_3_strings(string a, string b, string c) {
	int n = a.length();
	int m = b.length();
	int k = c.length();
	int dp[n+1][m+1][k+1];
	for(int i = 0; i <= n; i++) {
		for(int j = 0; j <= m; j++) {
			for(int l = 0; l <= k; l++) {
				if(i == 0 || j == 0 || l == 0) 
					dp[i][j][l] = 0;
				else if(a[i-1] == b[j-1] && a[i-1] == c[l-1])
					dp[i][j][l] = 1 + dp[i-1][j-1][l-1];
				else
					dp[i][j][l] = max(dp[i-1][j][l] , max(dp[i][j-1][l], dp[i][j][l-1]));
			}
		}
	}
	return dp[n][m][k];
}
//====================================================================================
int longest_bitonic_subsequence(int n, int a[]){
	int lis[n];
	int lds[n];
	for(int i = 0; i < n; i++) {
		lis[i] = 1;
		lds[i] = 1;
	}
	for(int i = 1; i < n; i++){
		for(int j = 0; j < i; j++){
			if(a[i] > a[j]) 
				lis[i] = max(lis[i] , 1 + lis[j]);
		}
	}
	for(int i = n-2; i >= 0; i--) {
		for(int j = n-1; j > i; j--) {
			if(a[i] > a[j]) {
				lds[i] = max(lds[i], 1 + lds[j]);
			}
		}
	}
	int mx = lis[0] + lds[0] - 1;
	for(int i = 1; i < n; i++){ 
		mx = max(mx, lis[i] + lds[i] - 1);
	}
	return mx;
}
//====================================================================================
int max_lis_sum(int n, int a[]) {
	int dps[n];

	for(int i = 0; i < n; i++)
		dps[i] = a[i];

	for(int i = 1; i < n; i++){
		for(int j = 0; j < i; j++){
			if(a[i] > a[j] && dps[i] < dps[j] + a[i])
				dps[i] = dps[j] + a[i];
		}
	}
	int max_sum = -1;
	for(int i = 0; i < n; i++) {
		max_sum = max(max_sum, dps[i]);
	}
	return max_sum;
}
//====================================================================================
int maximum_bitonic_sum(int n, int a[]) {
	int lis_sum[n];
	int lds_sum[n];

	for(int i = 0; i < n; i++) { // at just one element the sum will be equal to the num itself
		lis_sum[i] = a[i];
		lds_sum[i] = a[i];
	}

	for(int i = 1; i < n; i++) {
		for(int j = 0; j < i; j++){
			if(a[i] > a[j] && lis_sum[i] < lis_sum[j] + a[i])
				lis_sum[i] = lis_sum[j] + a[i];
		}
	}
	for(int i = n-2; i >= 0; i--) {
		for(int j = n-1; j > i; j--) {
			if(a[i] > a[j] && lds_sum[i] < lds_sum[j] + a[i])
				lds_sum[i] = lds_sum[j] + a[i];
		}
	}
	int max_sum = lis_sum[0] + lds_sum[0] - a[0];
	for(int i = 1; i < n; i++) 
		max_sum = max(max_sum, lis_sum[i] + lds_sum[i] - a[i]);
	return max_sum;
}
//====================================================================================
int maximum_product_increasing_subsequence(int n, int a[]) {
	int dps[n];

	for(int i = 0; i < n; i++)
		dps[i] = a[i];

	for(int i = 1; i < n; i++){
		for(int j = 0; j < i; j++){
			if(a[i] > a[j] && dps[i] < dps[j] * a[i])
				dps[i] = dps[j] * a[i];
		}
	}
	int max_sum = -1;
	for(int i = 0; i < n; i++) {
		max_sum = max(max_sum, dps[i]);
	}
	return max_sum;
}
//====================================================================================
int subsequences_product_less_than_K_recursion(int n, int a[], int k) {
	if(n == 0 || k == 0)
		return 0;
	if(a[n-1] > k) 
		return subsequences_product_less_than_K_recursion(n-1,a,k);
	return 1 + subsequences_product_less_than_K_recursion(n-1,a,k/a[n-1]) + 
				subsequences_product_less_than_K_recursion(n-1,a,k);
}
//====================================================================================
int subsequences_product_less_than_K(int n, int a[], int k) {
	int dp[n+1][k+1];
	for(int i = 0; i < n + 1; i++)
		dp[i][0] = 0;
	for(int i = 0; i < k + 1; i++)
		dp[0][i] = 0;

	for(int i = 1; i <= n; i++){
		for(int j = 1; j <= k; j++) {
			if(a[i-1] > j)
				dp[i][j] = dp[i-1][j];
			else
				dp[i][j] = 1 + dp[i-1][j/a[i-1]] + dp[i-1][j];
		}
	}
	return dp[n][k];
}
//====================================================================================
int longest_subsequence_difference_adjacents_is_1(int n, int a[]) {
	int dp[n];
	for(int i = 0; i < n; i++) 
		dp[i] = 1;
	for(int i = 1; i < n; i++) {
		for(int j = 0; j < i; j++ ){
			if(abs(a[i] - a[j]) == 1)
				dp[i] = max(dp[i], 1 + dp[j]);
		}
	}
	int mx = dp[0];
	for(int i = 1; i < n; i++) {
		mx = max(mx, dp[i]);
	}
	return mx;
}
//====================================================================================
int maximum_len_subseq_diff_bt_adj_0_or_1(int n, int a[]) {
	int dp[n];
	for(int i = 0; i < n; i++) 
		dp[i] = 1;
	for(int i = 1; i < n; i++) {
		for(int j = 0; j < i; j++) {
			if(abs(a[i] - a[j]) == 0 || abs(a[j] - a[i]) == 1) 
				dp[i] = max(dp[i] , 1 + dp[j]);
		}
	}
	int mx = dp[0];
	for(int i = 1; i < n; i++)
		mx  = max(mx, dp[i]);
	return mx;
}
//====================================================================================
int maximum_len_chain_pairs(vector<pair<int, int>> & vp) {
	int n = vp.size();
	int dp[n];
	for(int i = 0; i < n; i++) 
		dp[i] = 1;
	for(int i = 1; i < n; i++ ){
		for(int j = 0; j < i; j++) {
			if(vp[j].second < vp[i].first)
				dp[i] = max(dp[i], 1 + dp[j]);
		}
	}
	int mx = dp[0];
	for(int i = 1; i < n; i++) 
		mx = max(mx, dp[i]);
	return mx;
}
//====================================================================================
double maximum_average_value_in_grid(int n, int cost[][3]) {
	int dp[n][n];
	dp[0][0] = cost[0][0];
	for(int i = 1; i < n; i++) {
		dp[i][0] = dp[i-1][0] + cost[i][0];
		dp[0][i] = dp[0][i-1] + cost[0][i];
	}
	for(int i = 1; i < n; i++) {
		for(int j = 1; j < n; j++) {
			dp[i][j] = max(dp[i-1][j] , dp[i][j-1]) + cost[i][j];
		}
	}
	// maximum sum path has cost dp[n-1][n-1];
	return (double) dp[n-1][n-1] / (2 * n - 1) ;
}
//====================================================================================
int minimum_cost_path_matrix(int m, int n, int cost[][3]) {
	int dp[m][n];
	dp[0][0] = cost[0][0];
	for(int i = 1; i < m; i++) 
		dp[i][0] = dp[i-1][0] + cost[i][0];
	for(int j = 1; j < n; j++) 
		dp[0][j] = dp[0][j-1] + cost[0][j];

	for(int i = 1; i < m; i++) {
		for(int j = 1; j < n; j++){
			dp[i][j] = min(dp[i-1][j-1], min(dp[i-1][j], dp[i][j-1])) + cost[i][j];
		}
	}
	return dp[m-1][n-1];
}
//====================================================================================
// KADANE'S ALGORITHM 
// LONGEST SUM CONTIGOUS SUBARRAY
int kadanes(int n, int a[]) {	// finds longest positive sum subarray
	int max_so_far = 0;
	int max_ending_here = 0;

	for(int i = 0; i < n; i++) {
		max_ending_here += a[i];

		if(max_ending_here < 0) 
			max_ending_here = 0;

		if(max_so_far < max_ending_here) 
			max_so_far = max_ending_here;
	}
	return max_so_far;
}
//====================================================================================
int kadanes_only_negative(int n, int a[]) {
	int max_so_far = a[0];
	int max_ending_here = a[0];
	for(int i = 1; i < n; i++) {

	}
}
//====================================================================================
int maxsum_subarray_len(int n, int a[]) {
	int max_so_far = 0;
	int max_ending_here = 0;
	int start = 0; int end = 0; int restart = 0;

	for(int i = 0; i < n; i++){
		max_ending_here += a[i];

		if(max_so_far < max_ending_here) {
			start = restart;
			end = i;
		}
		if(max_ending_here < 0) {
			max_ending_here = 0;
			restart = i + 1;
		}
	}	
	return (end - start);
}
//====================================================================================
void print_kadanes(int n, int a[]) {
	// fst tym restart = 0 so that we have a positive sum subarrray 
	// everytime we get a negtve sum subarray we need to start over from the 
	// next index so -> restart = i + 1 and end is updated automatically
	int max_so_far = 0;
	int max_ending_here = 0;
	int start = 0; int end = 0; int restart = 0;
	for(int i = 0; i < n; i++) {
		max_ending_here += a[i];

		if(max_ending_here > max_so_far) { 
			max_so_far = max_ending_here;
			start = restart;
			end = i;
		}
		if(max_ending_here < 0) {
			max_ending_here = 0;
			restart = i + 1;
		}
	}
	for(int i = start; i <= end; i++) 
		cout << a[i] << " ";
	cout << endl;
}
//====================================================================================
int maxsum_pairs_with_diff_k(int n, int a[], int k ) {
	sort(a, a+n);
	int dp[n];
	dp[0] = 0;
	for(int i = 1; i < n; i++) {

		dp[i] = dp[i-1]; // no pairing with any element at first

		if(a[i] - a[i-1] < k) { // checking if pairing possible
			if(i >= 2)	
				dp[i] = max(dp[i], dp[i-2] + a[i] + a[i-1]);
			else
				dp[i] = max(dp[i], a[i] + a[i-1]);
		}
	}
	return dp[n-1];
}
//====================================================================================
// void max_size_square_submatirx_with_all_1(int r, int c, int a[][]) {
// 	int dp[r][c];
// 	for(int i = 0; i < r; i++) 
// 		dp[i][0] = a[i][0];
// 	for(int i = 0; i < c; i++)
// 		dp[0][i] = a[0][i];

// 	for(int i = 1; i < r; i++) {
// 		for(int j = 1; j < c; j++) {
// 			if(a[i][j] == 0)
// 				dp[i][j] = 0;
// 			else
// 				dp[i][j] = 1 + min(dp[i-1][j-1], min(dp[i-1][j], dp[i][j-1]));
// 		}
// 	}
// 	int max_size = dp[0][0]; int start = 0; int end = 0;
// 	for(int i = 0; i < r; i++) 
// 		for(int j = 0; j < c; j++) 
// 			if(dp[i][j] > max_size) {
// 				max_size = dp[i][j];
// 				start = i;
// 				end = j;
// 			}
// 	cout << "max size is : " << max_size << endl;
// 	for(int i = start; i > max_size - start; i--) {
// 		for(int j  = end; j > max_size - end; j--) 
// 			cout << a[i][j] << " ";
// 		cout << endl;
// 	}
// }
// //====================================================================================
// void max_size_square_submatirx_with_all_1_easy(int r, int c, int a[][C]){
// 	int dp[r][c];
// 	for(int i = 0; i < r; i++) 
// 		dp[i][0] = a[i][0];
// 	for(int i = 0; i < c; i++)
// 		dp[0][i] = a[0][i];

// 	int max_size = 0;
// 	for(int i = 1; i < r; i++) {
// 		for(int j = 1; j < c; j++) {
// 			if(a[i][j] == 0)
// 				dp[i][j] = 0;
// 			else {
// 				dp[i][j] = 1 + min(dp[i-1][j-1], min(dp[i-1][j], dp[i][j-1]));
// 				max_size = max(max_size, dp[i][j]);
// 			}
// 		}
// 	}
// 	cout << "max size is : " << max_size << endl;
// 	for(int i = 0; i < max_size; i++) {
// 		for(int j = 0; j < max_size; j++)
// 			cout << 1 << " ";
// 		cout << endl;
// 	}
// }
//====================================================================================
int break_sum(int n) {
	int dp[n+1];
	for(int i = 0; i < n + 1; i++) 
		dp[i] = i;
	for(int i = 1; i < n + 1; i++)
		dp[i] = max(dp[i], dp[i/2] + dp[i/3] + dp[i/4]);
	return dp[n];
}
//====================================================================================
int max_weight_path(vector<vector<int>> & a) {
	int n = a.size();
	int dp[n][n];
	memset(dp,0,sizeof(dp));

	for(int i = 0; i < n; i++) 
		dp[i][0] = a[i][0];
	
	for(int i = 1; i < n; i++) // can only be traversed down 
		dp[i][0] = dp[i-1][0] + a[i][0];
		
	for(int i = 1; i < n; i++){ 
		for(int j = 1; j < n; j++) {
			dp[i][j] = a[i][j] + max(dp[i-1][j-1], dp[i-1][j]);
		}
	}
	int mx = -1;
	// for(int i = 0; i < n; i++)
	// 	mx = max(mx, dp[n-1][i]);
	// for(int i= 0 ; i < n; i++) {
	// 	for(int j= 0; j < n; j++ )
	// 		cout << dp[i][j] << " ";
	// 	cout << endl;	
	// }
	return mx;
}
//====================================================================================
int maximum_sum_no2_are_adjacent(int n, int a[]) {
	if(n <= 0)
		return 0;
	return max(maximum_sum_no2_are_adjacent(n-1,a), a[n-1] + 
		maximum_sum_no2_are_adjacent(n-2,a));
}
//====================================================================================
int maximmum_sum_no2_are_adjacent_dp(int n, int a[]) {
	int dp[n] = {0};
	dp[0] = a[0];
	if(n > 0) 
		dp[1] = max(dp[0], dp[1]);
	for(int i = 2; i < n; i++) 
		dp[i] = max(dp[i-1], dp[i-2] + a[i]);
	return dp[n-1];
	// PUSH DP :
	// for(int i = 0; i + 2 < n; i++) 
	// 	dp[i+2] = max(dp[i] + a[i+2], dp[i+1]);
	// return dp[n-1];
}
//====================================================================================
int maxsum_2XN_grid(int n, vector<vector<int>> & mat) {
	if(n <= 0) 
		return 0;
	return max(max(maxsum_2XN_grid(n-1, mat), maxsum_2XN_grid(n-2, mat) + mat[0][n-1]),
				max(maxsum_2XN_grid(n-1, mat), maxsum_2XN_grid(n-2, mat) + mat[1][n-1]));
}
//====================================================================================
int maxsum_2XN_grid_DP(int n, vector<vector<int>> & mat) {
	int dp[n] = {0};
	dp[0] = max(mat[0][0], mat[1][0]);
	if(n > 1)
		dp[1] = max(dp[0], max(mat[0][1], mat[1][1]));
	for(int i = 2; i < n; i++) {
		dp[i] = max((dp[i-2] + max(mat[0][i], mat[1][i])), dp[i-1]);
	}
	return dp[n-1];
}
//====================================================================================
void maxsum_path_divisiblty(int n, int a[]) { // time : O(n * sqrt(n))
	int dp[n] = {0};
	for(int i = 0; i < n; i++) dp[i] = a[i];

	for(int i = 1; i < n; i++) {
		for(int j = 0; j * j < i; j++) {	// check this 
			if((i+1) % (j+1) == 0 && i > j) 
				dp[i] = max(dp[i] , a[i] + dp[j]);
		}
	}
	for(int i = 0; i < n; i++) 
		cout << dp[i] << " ";
	cout << endl;
}
//====================================================================================
int maxdiff_0_1(string s) {
	int n = s.length();
	int dp[n];
	for(int i = 0; i < n; i++) 
		dp[i] = s[i] == '0' ? 1 : -1;

	int max_ending_here = 0;
	int max_so_far = 0;

	for(int i = 0; i < n; i++) {
		max_ending_here += dp[i];
		if(max_ending_here < 0) 
			max_ending_here = 0;
		if(max_ending_here > max_so_far)
			max_so_far = max_ending_here;
	}
	return max_so_far > 0 ? max_so_far : -1;
}
//====================================================================================
int max_subarray_sum_K_times(int n, int a[], int k) {
	int max_so_far = 0; int max_ending_here = 0;

	for(int i = 0; i < n * k; i++){
		max_ending_here += a[i % n];	// like a circular array
		if(max_ending_here < 0) 
			max_ending_here = 0;
		if(max_ending_here > max_so_far)
			max_so_far = max_ending_here;
	}
	return max_so_far;
}
//====================================================================================
int maxpath_sum(vector<vector<int>> & mat){ 
	int r = mat.size();
	int c = mat[0].size();
	int dp[r][c];
	for(int i = 0; i < c; i++)
		dp[0][i] = mat[0][i];
	for(int i = 1; i < r; i++) {
		for(int j = 0; j < c; j++) {
			if(j > 0 && j + 1 < c)
				dp[i][j] = max(dp[i-1][j], max(dp[i-1][j-1], dp[i-1][j+1])) + mat[i][j];
			else if(j == 0)
				dp[i][j] = max(dp[i-1][j], dp[i-1][j+1]) + mat[i][j];
			else if(j == c-1) 
				dp[i][j] = max(dp[i-1][j], dp[i-1][j-1]) + mat[i][j];
		}
	}
	int ans = dp[r-1][0];
	for(int i = 1; i < c; i++) 
		ans = max(ans, dp[r-1][i]);
	return ans;
}
//====================================================================================
int minimum_jumps_to_end(int n, int a[]) {
	int dp[n];
	if(n == 0 || a[0] == 0) 
		return -1;
	// NOTE : if dp[n] = {INF} and then dp[0] = 0 then all dp[i] is automatically 
	// initialized to 0 . So, dp[i] = INF in the first loop is neccessary
	dp[0] = 0;
	// for(int i = 0; i < n; i++) 
	// 	cout << dp[i] << " ";
	// cout << endl;
	for(int i = 1; i < n; i++) {
		dp[i] = INF;		
		for(int j = 0; j < i; j++) {
			if(i <= a[j] + j) 
				dp[i] = min(dp[i] , 1 + dp[j]);
		}
		// deb(dp[i])
	}
	return dp[n-1];
}
//====================================================================================
int minimum_jumps_to_end_OPTIMIZED(int n, int a[]) {
	int dp[n];

	//TO BE  CONTINUED
	
}
//====================================================================================
int mincost_wt_bag(int W, int n, int val[]) {
	int wt[n] ;
	for(int i = 0; i < n; i++) 
		wt[i] = i + 1;
	
	int dp[n+1][W+1];
	for(int i = 0; i < n + 1; i++)
		dp[i][0] = 0;
	for(int i = 1; i < W + 1; i++)
		dp[0][i] = INF;

	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= W; j++) {
			if(j < wt[i-1] || val[i-1] == -1)
				dp[i][j] = dp[i-1][j];
			else
				dp[i][j] = min(dp[i-1][j], val[i-1] + dp[i][j-wt[i-1]]);
		}
	}
	
	return dp[n][W] == INF ? -1 : dp[n][W];
}
//====================================================================================
int minimum_removals_to_max_minus_min_lesequal_K(int n, int a[], int k) {
	sort(a, a + n);
	int ans = INF;
	for(int i = 0; i < n; i++) {
		int strt = i + 1;
		int end = n - 1;
		int id = -1;
		while(strt < end) {
			int mid = (strt + end) / 2;
			if(abs(a[mid] - a[i]) <= k) {
				id = mid;
				strt = mid + 1;
			} else 
				end = id + 1;
		}
		if(id != -1)
			ans = min(ans, n - (id - i + 1));
	}
	return ans;
}
//====================================================================================
int minimum_time_to_print_character_N_times(int n, int ins, int rem, int cop) {
	int dp[n+1];
	dp[0] = 0;
	dp[1] = ins;
	for(int i = 2; i <= n; i++) {
		dp[i] = INF;
		if(i & 1) 
			dp[i] = min(ins + dp[i-1], cop + rem + dp[(i+1)/2]);
		else
			dp[i] = min(ins + dp[i-1], cop + dp[i/2]);
	}
	return dp[n];
}
//====================================================================================
int minimum_steps(int n) {
	if(n == 1) return 0;
	int dp[n+1];
	for(int i = 0; i <= n; i++) 
		dp[i] = 0;

	for(int i = 2; i <= n; i++) {
	    if(i % 2 == 0 && i % 3 == 0)
	        dp[i] = 1 + min(dp[i-1], min(dp[i/2], dp[i/3]));
		else if(i % 2 == 0)
			dp[i] = 1 + min(dp[i-1], dp[i/2]);
		else if(i % 3 == 0)
			dp[i] = 1 + min(dp[i-1], dp[i/3]);
		else 
			dp[i] = 1 + dp[i-1];
	}
	return dp[n];
}
//====================================================================================
int longest_commong_substring(string a, string b) {
	int n = a.length(); int m = b.length();
	int dp[n+1][m+1];

	for(int i = 0; i < n + 1; i++)
		dp[i][0] = 0;
	for(int i = 0; i < m + 1; i++)
		dp[0][i] = 0;

	int result = 0;
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= m; j++) {
			if(a[i-1] == b[j-1])
				dp[i][j] = 1 + dp[i-1][j-1];
			else
				dp[i][j] = 0;
			result = max(result, dp[i][j]);
		}
	}
	// return dp[n][m];
	return result;
}
//====================================================================================
int moves(int n){
	int dp[n+1] = {0};
	dp[0] = 1;
	for(int i = 3; i <= n; i++)
		dp[i] += dp[i-3];
	for(int i = 5; i <= n; i++)
		dp[i] += dp[i-5];
	for(int i = 10; i <= n; i++)
		dp[i] += dp[i-10];
	return dp[n];
}
//====================================================================================
int staircase(int n) {
	// if(n <= 0)
	// 	return 0;
	// if(n == 1)
	// 	return 1;
	// if(n == 2)
	// 	return 2;
	// if(n == 3)
	// 	return 4;
	// return  staircase(n-1) + staircase(n-2) + staircase(n-3);
	int dp[n+1] = {0};
	dp[0] = 1; dp[1] = 1; dp[2] = 2; //dp[3] = 4;
	for(int i = 3; i <= n; i++)
		dp[i] += dp[i-1] + dp[i-2] + dp[i-3];
	return dp[n];
}
//====================================================================================
// VERY CONCEPTUAL 
// https://www.geeksforgeeks.org/count-ways-build-street-given-constraints/
ll build_street(ll n) {
	ll dp[n+1] = {0};
	dp[1] = 3;
	dp[2] = 7;
	for(ll i = 3; i <= n; i++) 
		dp[i] = 2 * dp[i-1] + dp[i-2];
	return dp[n];
}
//====================================================================================
int tiling_ways(int n, int m) {
	// ways to tile n X m floor by 1 X m tiles
	int dp[n+1] = {0};
	for(int i = 1; i <= n; i++){
		if(i > m)
			dp[i] = dp[i-1] + dp[i-m];
		else if(i < m)
			dp[i] = 1;
		else if(i == m)
			dp[i] = 2;
	}
	return dp[n];
}
//====================================================================================
int smallest_contigous_sum_subarray(int n, int a[]) {
	int minsofar = INF;
	int minendinghere = INF;
	for(int i = 0; i < n; i++) {
		if(minendinghere > 0)
			minendinghere = a[i];
		else 
			minendinghere += a[i];

		if(minsofar > minendinghere)
			minsofar = minendinghere;
	}
	return minsofar;
}
//====================================================================================
int repeated_deletion_lis_size(int n, int a[]) {
	bool vis[n];
	for(int i = 0; i < n; i++) 
		vis[i] = false;
	for(int i = 1; i < n; i++ ){
		for(int j = 0; j < i; j++ ){
			if(a[j] < a[i]) {
				vis[i] = true;
				vis[j] = true;
			}
		}
	}
	int cnt = 0;
	for(int i = 0; i < n; i++) 
		if(!vis[i])
			cnt++;
	return cnt == 0 ? -1 : cnt;
}
//====================================================================================
int utility(int dp[100][100], int a[], int low, int high, int turn) {
	if(low == high)
		return a[low] * turn;
	if(dp[low][high] != 0)
		return dp[low][high];

	dp[low][high] = max(a[low] * turn + utility(dp,a,low+1,high,turn+1), 
					a[high] * turn + utility(dp,a,low,high-1,turn+1));
	return dp[low][high];
}
int remove_array_element_to_maximize_product_of_sum(int n, int a[]) {
	int low = 0; int high = n-1; int turn = 1; int dp[100][100];
	memset(dp, 0 , sizeof(dp));
	return utility(dp,a,low,high,turn);
}
//====================================================================================
int longest_alternating_subarray(int n, int a[]) {
	int dp[n];
	dp[n-1] = 1;
	for(int i = n-2; i >= 0; i--) {
		if(a[i] * a[i+1] < 0)
			dp[i] = 1 + dp[i+1];
		else
			dp[i] = 1;
	}
	int mx = 1;
	for(int i = 0; i < n; i++) 
		mx = max(mx, dp[i]);
	return mx;
}
//====================================================================================
int sum_ways(int N, vector<int> & a) {
	int n = a.size();
	int dp[N+1];
	memset(dp, 0 , sizeof(dp));
	dp[0] = 1;
	for(int i = 1; i <= N; i++) {
		for(int j = 0; j < n; j++) {
			if(i >= a[j])
				dp[i] += dp[i-a[j]];
		}
	}
	return dp[N];
}
//====================================================================================
int ways_toget_sum_n_with_greater_than_equal_m(int n, int m) {
	if(n < m || n < 0)
		return 0;
	if(n == m) 
		return 1;
	// chose mth then recur for n-m, OR recur for n using greater than m
	return ways_toget_sum_n_with_greater_than_equal_m(n-m,m) +
		ways_toget_sum_n_with_greater_than_equal_m(n,m+1);
}
//====================================================================================

//====================================================================================

//======== DRIVER FUNCTION =======================================================================================
int32_t main(){
	// int s[] = {1,2,3};
	// int sum = 4;
	// deb(coin_chnge_1(s,3,sum));
	// string a = "aggtab";
	// int n = a.length();
	// string b = "gxtxayb";
	// int m = b.length();
	// deb(longest_common_subsequence(n,a,m,b));

	// int a[] =  {1,2,3,4};
	// int k = 10;
	// int n = sizeof(a)/sizeof(int);
	// deb(longest_increasing_subsequence_1(n,a));
	// deb(longest_increasing_subsequence_optimized(n,a));
	// deb(longest_bitonic_subsequence(n,a))
	// deb(maximum_bitonic_sum(n,a))
	// deb(max_lis_sum(n,a))
	// deb(maximum_product_increasing_subsequence(n,a))
	// deb(subsequences_product_less_than_K_recursion(n,a,k))
	// deb(subsequences_product_less_than_K(n,a,k))

	// int a[] =   {-2, -3, 4, -1, -2, 1, 5, -3}; 
	// int n = sizeof(a)/sizeof(int);
	// deb(longest_subsequence_difference_adjacents_is_1(n,a))
	// deb(maximum_len_subseq_diff_bt_adj_0_or_1(n,a))

	// vector<pair<int,int>> vp = {{5, 24}, {39, 60}, {15, 28}, {27, 40}, {50, 90} };
	// deb(maximum_len_chain_pairs(vp))
	// int cost[][3] = {{1,2,3},{4,8,2},{1,5,3}};
	// int m = 3; int n = 3;
	// deb(minimum_cost_path_matrix(m,n,cost))

	// deb(kadanes(n,a))
	// deb(maxsum_subarray(n,a))
	// print_kadanes(n,a);
	// deb(maxsum_subarray_len(n,a))
	// deb(break_sum(24))
	// vector<vector<int>> a = {{ 4, 2 ,3 ,4 ,1 },		
	//                    		 { 2 , 9 ,1 ,10 ,5 },
	//                    		 {15, 1 ,3 , 0 ,20 },
	//                    		 {16 ,92, 41, 44 ,1},
	//                    		 {8, 142, 6, 4, 8} };
	//    deb(max_weight_path(a))

	// int a[] = {6, 7, 1, 3, 8, 2, 4};
	// int n = sizeof(a) / sizeof(int);
	// deb(maximum_sum_no2_are_adjacent(n,a))
	// deb(maximmum_sum_no2_are_adjacent_dp(n,a))

	// vector<vector<int>> mat = {{1,2,3,4,5},
	// 							{6,7,8,9,10}};
	// deb(maxsum_2XN_grid(5,mat))
	// deb(maxsum_2XN_grid_DP(5,mat))
	// int a[] = {2, 3, 1, 4, 6, 5};
	// int n = sizeof(a) / sizeof(int);
	// maxsum_path_divisiblty(n,a);

	// string s = "11000010001";
	// deb(maxdiff_0_1(s))

	// vector<vector<int>> mat = { {4, 2, 3, 4},
 	//                      {2, 9, 1, 10},
 	//                      {15, 1, 3, 0},
 	//                      {16, 92, 41, 44} };
 	//    deb(maxpath_sum(mat))

	// int a[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	// int n = sizeof(a) / sizeof(int);
	// deb(minimum_jumps_to_end(n,a))

	// int W = 5;
	// int val[] =  {-1, -1, 4, 5, -1};
	// int n = sizeof(val) / sizeof(int);
	// deb(mincost_wt_bag(W,n,val));

	// int a[] = {1, 3, 4, 9, 10, 11, 12, 17, 20};
	// int n = sizeof(a) / sizeof(int);	 int k = 4;
	// deb(minimum_removals_to_max_minus_min_lesequal_K(n,a,k))

	// int n = 9;
	// int insrt_time = 1;
	// int remov_time = 2;
	// int copy_time = 1;
	// deb(minimum_time_to_print_character_N_times(n,insrt_time,remov_time,copy_time))

	// int n = 6;
	// deb(minimum_steps(n))
	// string a =  "zxabcdezy";
	// string b = "yzabcdezx";
	// deb(longest_commong_substring(a,b))
	// deb(moves(13))
	// deb(staircase(4))
	// deb(build_street(34))
	// deb(tiling_ways(19,4))
	// int a[] =  {1, 2, 5, 3, 6, 4, 1};
	// int n = sizeof(a) / sizeof(int);
	// deb(repeated_deletion_lis_size(n,a))
	// int a[] = {-5, -1, -1, 2, -2, -3};
	// deb(longest_alternating_subarray(6,a))

	// vector<int> a = {1,5,6};
	// int N = 7;
	// deb(sum_ways(N, a))
	
	// deb(ways_toget_sum_n_with_greater_than_equal_m(2,1))

	return 0;
}