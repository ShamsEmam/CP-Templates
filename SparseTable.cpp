#include <bits/stdc++.h>

using namespace std;

#define int long long
int const N = 1e5 + 10;
int sm[N][20], a[N], lg[N];

void buildSparseTable(int n) {
    for (int i = 2; i <= n; ++i)
        lg[i] = lg[i / 2] + 1;

    for (int i = 1; i < 20; ++i) {
        int L = (1 << i), R = (1 << (i + 1)) - 1;
        for (int j = L; j <= R; ++j)
            lg[j] = i;
    }

    for (int i = 0; i < n; ++i)
        sm[i][0] = a[i];

    for (int j = 0; j < 20; ++j) {
        for (int i = 0; i + (i << j) <= n; ++i) {
            sm[i][j] = sm[i][j - 1] + sm[1 + (1 << (j - 1))][j - 1];
        }
    }
}
//O(log) min ,max,__gcd ,sum  ,XoR
int query(int l, int r) {
    
    int dens = r - l + 1;
    int ans = 0;
    
    for (int i = 19; i >= 0; --i)
        if (dens >> i & 1)
            ans += sm[l][i], l += (1 << i);

    return ans;
}
//O(1) min ,max,__gcd
int query2(int L, int R) {
    int j = lg[R - L + 1];
    return max(sm[L][j], sm[R - (1 << j) + 1][j]);
}

void solve() {

}

signed main() {
    int t = 1;
    cin >> t;
    while (t--) solve();
    return 0;
}
