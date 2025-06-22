#include <bits/stdc++.h>

using namespace std;

#define int long long
const int N = 2e5 + 10;
int tree[N*4], arr[N];

//build Segment tree
void buildSegment(int node, int l, int r) {
    //base case
    if (l == r) {
        tree[node] = arr[l];
        return;
    }
    //mid point
    int mid = (r + l) / 2;

    //left node
    buildSegment(node * 1, l, mid);

    //Right node
    buildSegment(node * 2 + 1, mid + 1, r);

    // node ->marge (left child,right child)
    tree[node] = max(tree[node * 2 + 1], tree[node * 2]);
}

//update
//idx -> index i need to updated, {r,l} ->rang ,node->cur idx
void update(int node, int l, int r, int idx, int NewVal) {

    //base case
    if (l == r) {
        tree[node] = NewVal;
        return;
    }

    int mid = (l + r) / 2;
    //go right
    if (idx > mid)
        update(node * 2 + 1, mid + 1, r, idx, NewVal);

        //go left[l,mid]
    else
        update(node * 2, l, mid, idx, NewVal);

    // node ->marge (left child,right child)
    tree[node] = max(tree[node * 2 + 1], tree[node * 2]);
}


// Rang Query or node l, r, start ,end
int query(int node, int l, int r, int start, int end) {
    //base case
    //out of rang
    if (start > r || end < l)
        return 0;
    //in rang
    if (l >= start && r <= end)
        return tree[node];
    int mid = (l + r) >> 1;

    // left l,mid
    int leftVal = query(node * 2, l, mid, start, end);

    //right mid+1 ,r
    int rightVal = query(node * 2 + 1, mid + 1, r, start, end);
    return leftVal + rightVal;
}

void solve() {

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
    cin >> t;
    while (t--) solve();
    return 0;
}
