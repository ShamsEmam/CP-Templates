#include <bits/stdc++.h>

using namespace std;

#define ll long long
#define all(s) s.begin(), s.end()
#define sz(s) (int)(Qry).size()
const int B = 337, N = 1e6 + 5;

struct Qry {
    int l, r, i;

    bool operator<(const Qry &other) const {
        // index block = l / sz (block)
        return make_pair(l / B, r) < make_pair(other.l / B, other.r);

    };
};

void solve() {
    int n, q;
    cin >> n;
    vector<int> v(n);
    for (auto &i: v)cin >> i;
    cin >> q;
    vector<Qry> queries(q);
    for (int i = 0; i < q; ++i) {
        int l, r;
        cin >> l >> r;
        queries[i] = {--l, --r, i};
    }
    sort(all(queries));

    vector<int> frq(N), sum(N), frq_frq(N);
    int mx = -1;
    auto add = [&](int idx) {
        int x = v[idx];
        mx = max(mx, ++frq[x]);
        frq_frq[frq[x]]++;
        sum[frq[x]] += x;
    };
    auto erase = [&](int idx) {
        int x = v[idx];
        sum[frq[x]] -= x;
        frq_frq[frq[x]]--;
        --frq[x];
        if (!sum[mx])--mx;
    };
    auto get_answer = [&]() {
        return sum[mx];
    };
    int l = 0, r = -1;
    vector<int> res(q);
    for (const auto &[lq, rq, iq]: queries) {
        while (lq < l) add(--l);
        while (r < rq) add(++r);
        while (l < lq) erase(l++);
        while (rq < r) erase(r--);

        res[iq] = get_answer();
    }
    for (auto &x: res)
        cout << x << '\n';

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

#undef int

int main() {
    setupIO();
    speed();
    int t = 1;
//    cin >> t;
    while (t--) solve();
    return 0;
}
