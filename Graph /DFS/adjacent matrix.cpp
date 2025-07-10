#include <bits/stdc++.h>

#define int long long
using namespace std;
const int MOD = 1e9 + 7;
const int N = 1e3 + 7;

//Memory o(n^2)
//tie to check is connected o(1) 
int g[N][N];

bool CheckIsConnected(int i, int j) {
    return (g[i][j] || g[j][i]);
}

void solve() {
    int n, m;
    cin >> n >> m;

    int u, v;
    for (int i = 1; i <= m; ++i) {
        cin >> u >> v;
        g[u][v] = 1;
        g[v][u] = 1;
    }
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            cout << g[i][j] << ' ';
        }
        cout << '\n';
    }
    cin >> u >> v;
    if (CheckIsConnected(u, v))cout << "YES\n";
    else cout << "NO\n";
}

void speed() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
}

void setupIO() {
#ifndef ONLINE_JUDGE
    freopen("in.txt", "r", stdin);
    freopen("out.txt", "w", stdout);
#endif
}

signed main() {
    setupIO();
    speed();
    int t = 1;
//    cin >> t;
    while (t--) solve();
    return 0;
}
