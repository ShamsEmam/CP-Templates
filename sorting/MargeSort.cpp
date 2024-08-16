#include <bits/stdc++.h>

#define fast ios::sync_with_stdio(0), cin.tie(0), cout.tie(0);
#define int long long
#define nl "\n"
using namespace std;

void marg(vector<int> &vec, int l, int r) {
    queue<int> left, right;
    int mid = (l + r) / 2;

    for (int i = l; i <= mid; ++i) {
        left.push(vec[i]);
    }
    for (int i = mid + 1; i <= r; ++i) {
        right.push(vec[i]);
    }
    for (int i = l; i < r + 1; ++i) {

        if (left.empty()) {
            vec[i] = right.front();
            right.pop();
        } else if (right.empty()) {
            vec[i] = left.front();
            left.pop();
        } else if (left.front() < right.front()) {
            vec[i] = left.front();
            left.pop();
        } else {
            vec[i] = right.front();
            right.pop();
        }
    }
}

void margSort(vector<int> &vec, int l, int r) {
    if (l >= r)
        return;
    int mid = (l + r) / 2;
    margSort(vec, l, mid);
    margSort(vec, mid + 1, r);
    marg(vec, l, r);

}

void solve() {
    int n;
    cin >> n;
    vector<int> v(n);
    for (int i = 0; i < n; ++i) {
        cin >> v[i];
    }
    margSort(v, 0, n - 1);
    for (auto x: v)
        cout << x << " ";
}


signed main() {
 
    fast;
    int t;
    t = 1;
//    cin >> t;
    while (t--) {
        solve();
    }
    return 0;
}
