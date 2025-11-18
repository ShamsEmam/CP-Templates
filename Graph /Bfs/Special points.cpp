#include <bits/stdc++.h>

using namespace std;

#define ll long long
#define all(s) s.begin(), s.end()
#define sz(s) (int)(s).size()
int const N = 1e6 + 5;

vector<int> adj[N];
bool vis[N];
map<int, int> is_special;
int n, m;
vector<int> lvl, InDegree, Topological;

void bfs() {
    queue<int> q;
    for (int start = 1; start <= n; start++) {
        if (is_special[start]) {
            q.push(start);
            lvl[start] = 0;
        }
    }


    while (!q.empty()) {
        int node = q.front();
        q.pop();

        for (auto child: adj[node]) {
            InDegree[child]--;
            if (is_special[child]) {
                lvl[child] = lvl[node] + 1;
                q.push(child);

            }
        }
    }
}


void solve() {
    cin >> n >> m;

    for (int i = 1; i <= n; i++) {
        adj[i].clear();
    }

    lvl.assign(n + 1, 0);
    InDegree.assign(n + 1, 0);
    Topological.clear();

    int u, v;
    for (int i = 0; i < m; ++i) {
        cin >> u >> v;
        adj[u].push_back(v);
        InDegree[v]++;
    }

    bfs();
    int szSpecial;
    cin >> szSpecial;
    for (int i = 0, sp; i < szSpecial; ++i) {
        cin >> sp;
        is_special[sp] = 1;
    }
    int q;
    cin >> q;
    while (q--) {
        int x;
        cin >> x;
        cout<<lvl[x]<<'\n';
    }

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
