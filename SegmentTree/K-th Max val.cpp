#include <bits/stdc++.h>

using namespace std;

#define ll long long
#define ld long double
#define all(s) s.begin(), s.end()
#define sz(s) (int)(s).size()
ll const oo = 1e15;
int const K = 10;

struct Node {
    vector<int> v;

    Node() {}

    Node(int val) {
        v = {val};
    }
};

struct SegTree {
#define mid ((lx+rx)/2)
#define lf (2*x+1)
#define rt (2*x+2)
#define neutral 0
    vector<Node> seg;
    int n;

    void init(vector<int> &v) {
        n = sz(v);
        seg.assign(2 * n, {neutral});
        build(v, 0, 0, n - 1);
    }

    Node mrg(const Node &Lf, const Node &Rt) {
        int i = 0, j = 0;
        vector<int> mxs;
        int limit = min({K, sz(Lf.v) + sz(Rt.v)});

        auto &a = Lf.v;
        auto &b = Rt.v;

        vector<int> ans;
        while (sz(ans) < limit) {
            if (j < sz(b) and (i == sz(a) or b[j] > a[i])) {
                ans.push_back(b[j++]);
            } else {
                ans.push_back(a[i++]);
            }
        }
        Node ret;
        ret.v = ans;
        return ret;

    }

    void build(vector<int> &v, int x, int lx, int rx) {
        if (lx == rx) {
            seg[x] = Node(v[lx]);
            return;
        }
        build(v, lf, lx, mid);
        build(v, rt, mid + 1, rx);
        seg[x] = mrg(seg[lf], seg[rt]);
    }

    void update(int idx, int val, int x, int lx, int rx) {
        if (lx == rx) {
            seg[x] = {val};
            return;
        }
        if (idx <= mid)
            update(idx, val, lf, lx, mid);
        else
            update(idx, val, rt, mid + 1, rx);

        seg[x] = mrg(seg[lf], seg[rt]);
    }

    void update(int idx, int val) {
        return update(idx, val, 0, 0, n - 1);
    }

    Node query(int lq, int rq, int x, int lx, int rx) {
        if (rx < lq || lx > rq)
            return {};
        if (lx >= lq && rx <= rq)
            return seg[x];

        query(lq, rq, lf, lx, mid);
        query(lq, rq, rt, mid + 1, rx);

        return mrg(query(lq, rq, lf, lx, mid),
                   query(lq, rq, rt, mid + 1, rx));
    }

    Node query(int lq, int rq) {
        return query(lq, rq, 0, 0, n - 1);
    }

#undef mid
#undef lf
#undef rt
#undef neutral
};

void solve() {

    int n, q;
    cin >> n >> q;
    vector<int> v(n);
    SegTree st;
    for (int i = 0; i < n; ++i) {
        cin >> v[i];
    }
    st.init(v);
    int op, val, idx, l, r, k;
    while (q--) {
        cin >> op;
        if (op == 1) {
            cin >> idx >> val;
            st.update(idx, val);
        } else {
            cin >> l >> r >> k;
            cout << st.query(l, r - 1).v[--k] << '\n';
        }
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
//    cout << fixed << setprecision(15);
//    cin >> t;
    while (t--) solve();
    return 0;
}
