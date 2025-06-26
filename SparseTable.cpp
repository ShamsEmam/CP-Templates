#include <bits/stdc++.h>

using namespace std;

#define int long long
const int N = 1e6 + 10;
int sm[N][20], a[N], lg[N];

void buildSparseTable(int n) {

    for (int i = 2; i <= n; ++i)
        lg[i] = lg[i / 2] + 1;

    for (int i = 0; i < n; ++i)
        sm[i][0] = a[i];

    for (int j = 1; (1 << j) <= n; ++j) {
        for (int i = 0; i + (1 << j) <= n; ++i) {
            sm[i][j] = max(sm[i][j - 1], sm[i + (1 << (j - 1))][j - 1]);
        }
    }
}

// O(1) max query
int query2(int L, int R) {
    int j = lg[R - L + 1];
    return max(sm[L][j], sm[R - (1 << j) + 1][j]);
}

void solve() {
 
}


void speed() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
}

signed main() {
    speed();
    int t = 1;
    while (t--) solve();
    return 0;
}
