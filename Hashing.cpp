#include <bits/stdc++.h>
using namespace std;

#define ll long long

const int N = 2005;
const int MOD1 = 1e9 + 7;
const int MOD2 = 1e9 + 9;

int pw1[N], pw2[N];
int inv1[N], inv2[N];
int BASE;

// ------------------ Utilities ------------------

int fix(ll x, int mod) {
    return (x % mod + mod) % mod;
}

int fpow(int a, int b, int mod) {
    if (b == 0) return 1;
    int ret = fpow(a, b / 2, mod);
    ret = fix(1LL * ret * ret, mod);
    if (b & 1) ret = fix(1LL * ret * a, mod);
    return ret;
}

bool isPrime(int x) {
    for (int i = 2; i * i <= x; i++) {
        if (x % i == 0) return false;
    }
    return x > 1;
}

// ------------------ Init ------------------

void init() {
    static bool done = false;
    if (done) return;
    done = true;

    mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());
    uniform_int_distribution<int> dist(257, 1000);

    do {
        BASE = dist(rng);
    } while (!isPrime(BASE));

    pw1[0] = pw2[0] = 1;
    inv1[0] = inv2[0] = 1;

    int invBase1 = fpow(BASE, MOD1 - 2, MOD1);
    int invBase2 = fpow(BASE, MOD2 - 2, MOD2);

    for (int i = 1; i < N; i++) {
        pw1[i] = fix(1LL * pw1[i - 1] * BASE, MOD1);
        pw2[i] = fix(1LL * pw2[i - 1] * BASE, MOD2);

        inv1[i] = fix(1LL * inv1[i - 1] * invBase1, MOD1);
        inv2[i] = fix(1LL * inv2[i - 1] * invBase2, MOD2);
    }
}

// ------------------ Hash Struct ------------------

struct Hash {
    vector<pair<int,int>> pre;

    Hash(const string &s) {
        init();
        int n = s.size();
        pre.assign(n + 1, {0, 0});

        for (int i = 0; i < n; i++) {
            pre[i + 1].first = fix(pre[i].first + 1LL * pw1[i] * s[i], MOD1);
            pre[i + 1].second = fix(pre[i].second + 1LL * pw2[i] * s[i], MOD2);
        }
    }

    // get hash of substring [l, r]
    pair<int,int> get(int l, int r) {
        ll h1 = pre[r + 1].first - pre[l].first;
        ll h2 = pre[r + 1].second - pre[l].second;

        h1 = fix(h1, MOD1);
        h2 = fix(h2, MOD2);

        h1 = fix(h1 * inv1[l], MOD1);
        h2 = fix(h2 * inv2[l], MOD2);

        return {h1, h2};
    }
};
