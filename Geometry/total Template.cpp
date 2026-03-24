

// =========================================
//           GEOMETRY LIBRARY
//      Organized by Shams Emam
// =========================================

#include <bits/stdc++.h>
using namespace std;

#define EPS 1e-9
#define ld long double

// ==================== BASICS ====================
typedef ld T;
typedef complex<T> pt;
#define x real()
#define y imag()

bool operator==(pt a, pt b) { return fabs(a.x - b.x) < EPS && fabs(a.y - b.y) < EPS; }
bool operator!=(pt a, pt b) { return !(a == b); }

T sq(pt p) { return p.x*p.x + p.y*p.y; }
ld abs(pt p) { return sqrt(sq(p)); }
pt perp(pt p) { return {-p.y, p.x}; }

T dot(pt v, pt w) { return v.x*w.x + v.y*w.y; }
T cross(pt v, pt w) { return v.x*w.y - v.y*w.x; }

int sgn(T val) { return (val > EPS) - (val < -EPS); }
bool isPerp(pt v, pt w) { return fabs(dot(v, w)) < EPS; }

// ==================== TRANSFORMATIONS ====================
pt translate(pt v, pt p) { return p + v; }

pt scale(pt c, T factor, pt p) {
    return c + (p-c)*factor;
}

pt rot(pt p, pt c, ld a) {
    pt v = p - c;
    return {c.x + v.x*cos(a) - v.y*sin(a), c.y + v.x*sin(a) + v.y*cos(a)};
}

pt linearTransfo(pt p, pt q, pt r, pt fp, pt fq) {
    pt pq = q - p, num{cross(pq, fq-fp), dot(pq, fq-fp)};
    return fp + pt{cross(r-p, num), dot(r-p, num)} / sq(pq);
}

// ==================== ANGLES ====================
T orient(pt a, pt b, pt c) { return cross(b-a, c-a); }

ld angle(pt v, pt w) {
    return acos(clamp(dot(v,w)/abs(v)/abs(w), (T)-1.0, (T)1.0));
}

ld orientedAngle(pt a, pt b, pt c) {
    if (orient(a,b,c) >= 0)
        return angle(b-a, c-a);
    else
        return 2*M_PI - angle(b-a, c-a);
}

ld angleTravelled(pt a, pt p, pt q) {
    double ampli = angle(p-a, q-a);
    if (orient(a,p,q) > 0) return ampli;
    else return -ampli;
}

bool inAngle(pt a, pt b, pt c, pt p) {
    T abp = orient(a,b,p), acp = orient(a,c,p), abc = orient(a,b,c);
    if (abc < 0) swap(abp, acp);
    return (abp >= 0 && acp <= 0) ^ (abc < 0);
}

bool half(pt p) { return p.y > 0 || (p.y == 0 && p.x < 0); }

// ==================== LINES ====================
struct line {
    pt v; T c;
    line(pt v, T c) : v(v), c(c) {}
    line(T a, T b, T _c) { v = {b,-a}; c = _c; }
    line(pt p, pt q) { v = q-p; c = cross(v,p); }

    T side(pt p) { return cross(v,p)-c; }
    double dist(pt p) { return abs(side(p))/abs(v); }
    double sqDist(pt p) { return side(p)*side(p)/sq(v); }
    line perpThrough(pt p) { return {p, p + perp(v)}; }
    bool cmpProj(pt p, pt q) { return dot(v,p) < dot(v,q); }
    line translate(pt t) { return {v, c + cross(v,t)}; }
    line shiftLeft(double dist) { return {v, c + dist*abs(v)}; }
    pt proj(pt p) { return p - perp(v)*side(p)/sq(v); }
    pt refl(pt p) { return p - perp(v) * (T)2.0 * side(p)/sq(v); }
};

bool inter(line l1, line l2, pt &out) {
    T d = cross(l1.v, l2.v);
    if (fabs(d) <= EPS) return false;
    out = (l2.v*l1.c - l1.v*l2.c)/d;
    return true;
}

line bisector(line l1, line l2, bool interior) {
    assert(cross(l1.v, l2.v) != 0);
    T sign = interior ? 1 : -1;
    return {l2.v/abs(l2.v) + l1.v/abs(l1.v) * sign, l2.c/abs(l2.v) + l1.c/abs(l1.v) * sign};
}

// ==================== SEGMENTS ====================
bool inDisk(pt a, pt b, pt p) { return dot(a-p, b-p) <= EPS; }
bool onSegment(pt a, pt b, pt p) { return fabsl(orient(a,b,p)) <= EPS && inDisk(a,b,p); }

bool properInter(pt a, pt b, pt c, pt d, pt &out) {
    T oa = orient(c,d,a), ob = orient(c,d,b), oc = orient(a,b,c), od = orient(a,b,d);
    if (sgn(oa)*sgn(ob) < 0 && sgn(oc)*sgn(od) < 0) {
        out = (a*ob - b*oa)/(ob-oa);
        return true;
    }
    return false;
}

set<pair<ld,ld>> inters(pt a, pt b, pt c, pt d) {
    set<pair<ld,ld>> s; pt out;
    if(a == c || a == d) s.insert({a.x, a.y});
    if(b == c || b == d) s.insert({b.x, b.y});
    if(s.size()) return s;
    if(properInter(a,b,c,d,out)) return {{out.x, out.y}};
    if(onSegment(c,d,a)) s.insert({a.x, a.y});
    if(onSegment(c,d,b)) s.insert({b.x, b.y});
    if(onSegment(a,b,c)) s.insert({c.x, c.y});
    if(onSegment(a,b,d)) s.insert({d.x, d.y});
    return s;
}

ld segPoint(pt a, pt b, pt p) {
    if(a != b) { line l(a,b); if(l.cmpProj(a,p) && l.cmpProj(p,b)) return l.dist(p); }
    return min(abs(p-a), abs(p-b));
}

ld segSeg(pt a, pt b, pt c, pt d) {
    pt dummy; if(properInter(a,b,c,d,dummy)) return 0;
    return min({segPoint(a,b,c), segPoint(a,b,d), segPoint(c,d,a), segPoint(c,d,b)});
}

// ==================== POLYGONS ====================
bool isConvex(vector<pt> p) {
    bool hasPos=false, hasNeg=false;
    for(int i=0,n=p.size();i<n;i++){
        int o = sgn(orient(p[i], p[(i+1)%n], p[(i+2)%n]));
        if(o>0) hasPos=true; if(o<0) hasNeg=true;
    }
    return !(hasPos && hasNeg);
}

ld areaTriangle(pt a, pt b, pt c) { return abs(cross(b-a, c-a))/2.0; }

ld areaPolygon(vector<pt> p){
    ld area=0;
    for(int i=0,n=p.size();i<n;i++) area += cross(p[i], p[(i+1)%n]);
    return abs(area)/2.0;
}

bool above(pt a, pt p) { return p.y >= a.y; }
bool crossesRay(pt a, pt p, pt q) { return (above(a,q)-above(a,p))*orient(a,p,q) > 0; }

bool inPolygon(vector<pt> p, pt a, bool strict=true) {
    int numCrossings=0;
    for(int i=0,n=p.size();i<n;i++){
        if(onSegment(p[i],p[(i+1)%n],a)) return !strict;
        numCrossings += crossesRay(a,p[i],p[(i+1)%n]);
    }
    return numCrossings&1;
}

// ==================== CONVEX HULL ====================
bool cw(pt a, pt b, pt c, bool include_collinear) {
    int o = sgn(orient(a,b,c));
    return o<0 || (include_collinear && o==0);
}

bool collinear(pt a, pt b, pt c) { return sgn(orient(a,b,c))==0; }

void convex_hull(vector<pt> &a, bool include_collinear=false) {
    pt p0 = *min_element(a.begin(),a.end(),[](pt a, pt b){return make_pair(a.y,a.x) < make_pair(b.y,b.x);});
    sort(a.begin(), a.end(), [&p0](const pt &a, const pt &b){
        int o = sgn(orient(p0,a,b));
        if(o==0) return (p0.x-a.x)*(p0.x-a.x)+(p0.y-a.y)*(p0.y-a.y) < (p0.x-b.x)*(p0.x-b.x)+(p0.y-b.y)*(p0.y-b.y);
        return o<0;
    });
    if(include_collinear){
        int i=(int)a.size()-1;
        while(i>=0 && collinear(p0,a[i],a.back())) i--;
        reverse(a.begin()+i+1,a.end());
    }
    vector<pt> st;
    for(int i=0;i<(int)a.size();i++){
        while(st.size()>1 && !cw(st[st.size()-2], st.back(), a[i], include_collinear)) st.pop_back();
        if(st.empty() || a[i]!=st.back()) st.push_back(a[i]);
    }
    if(!include_collinear && st.size()==2 && st[0]==st[1]) st.pop_back();
    a = st;
}

// ==================== MINKOWSKI SUM ====================
void reorder_polygon(vector<pt> &P){
    size_t pos=0;
    for(size_t i=1;i<P.size();i++) if(P[i].y<P[pos].y || (P[i].y==P[pos].y && P[i].x<P[pos].x)) pos=i;
    rotate(P.begin(),P.begin()+pos,P.end());
}

vector<pt> minkowski(vector<pt> P, vector<pt> Q){
    reorder_polygon(P); reorder_polygon(Q);
    P.push_back(P[0]); P.push_back(P[1]);
    Q.push_back(Q[0]); Q.push_back(Q[1]);
    vector<pt> result; size_t i=0,j=0;
    while(i<P.size()-2 || j<Q.size()-2){
        result.push_back(P[i]+Q[j]);
        auto crossV = cross(P[i+1]-P[i], Q[j+1]-Q[j]);
        if(crossV>=0 && i<P.size()-2) i++;
        if(crossV<=0 && j<Q.size()-2) j++;
    }
    return result;
}

// ==================== HALF-PLANE INTERSECTION ====================
struct Halfplane {
    pt p, pq; ld angle;
    Halfplane() {}
    Halfplane(const pt &a, const pt &b) : p(a), pq(b-a) { angle=atan2l(pq.y,pq.x); }
    bool out(const pt &r){ return cross(pq,r-p)<-EPS; }
    bool operator<(const Halfplane &e) const { return angle<e.angle; }
    friend pt inter(const Halfplane &s,const Halfplane &t){ return s.p + s.pq * (cross(t.p-s.p,t.pq)/cross(s.pq,t.pq)); }
};

vector<pt> hp_intersect(vector<Halfplane> &H){
    const int inf = 1e9;
    pt box[4] = { {inf,inf},{-inf,inf},{-inf,-inf},{inf,-inf} };
    for(int i=0;i<4;i++) H.push_back(Halfplane(box[i],box[(i+1)%4]));
    sort(H.begin(),H.end());
    deque<Halfplane> dq; int len=0;
    for(int i=0;i<(int)H.size();i++){
        while(len>1 && H[i].out(inter(dq[len-1],dq[len-2]))){ dq.pop_back(); --len; }
        while(len>1 && H[i].out(inter(dq[0],dq[1]))){ dq.pop_front(); --len; }
        if(len>0 && fabsl(cross(H[i].pq,dq[len-1].pq))<EPS){
            if(dot(H[i].pq,dq[len-1].pq)<0.0) return vector<pt>();
            if(H[i].out(dq[len-1].p)){ dq.pop_back(); --len; } else continue;
        }
        dq.push_back(H[i]); ++len;
    }
    while(len>2 && dq[0].out(inter(dq[len-1],dq[len-2]))) { dq.pop_back(); --len; }
    while(len>2 && dq[len-1].out(inter(dq[0],dq[1]))) { dq.pop_front(); --len; }
    if(len<3) return vector<pt>();
    vector<pt> ret(len);
    for(int i=0;i+1<len;i++) ret[i]=inter(dq[i],dq[i+1]);
    ret.back()=inter(dq[len-1],dq[0]);
    return ret;
}

// ==================== ANTI-PODAL PAIRS ====================
vector<pair<int,int>> all_anti_podal(int n, vector<pt> &p){
    vector<pair<int,int>> result;
    auto nx = [&](int i){return (i+1)%n;};
    auto pv = [&](int i){return (i-1+n)%n;};
    vector<bool> vis(n,false);
    for(int p1=0,p2=0;p1<n;p1++){
        pt base = p[nx(p1)]-p[p1];
        while(p2==p1 || p2==nx(p1) || sgn(cross(base,p[nx(p2)]-p[p2]))==sgn(cross(base,p[p2]-p[pv(p2)]))) p2=nx(p2);
        if(vis[p1]) continue;
        vis[p1]=true;
        result.push_back({p1,p2});
        result.push_back({nx(p1),p2});
        if(sgn(cross(base,p[nx(p2)]-p[p2]))==0){ result.push_back({p1,nx(p2)}); result.push_back({nx(p1),nx(p2)}); vis[p2]=true; }
    }
    return result;
}

// ==================== MAX DISTANCE BETWEEN CONVEX POLYGONS ====================
double maximum_dist_from_polygon_to_polygon(vector<pt> &u, vector<pt> &v){
    int n=(int)u.size(), m=(int)v.size(); ld ans=0;
    if(n<3 || m<3){ for(int i=0;i<n;i++) for(int j=0;j<m;j++) ans=max(ans,sq(u[i]-v[j])); return sqrt(ans);}
    if(u[0].x>v[0].x) swap(n,m), swap(u,v);
    int i=0,j=0,step=n+m+10;
    while(j+1<m && v[j].x<v[j+1].x) j++;
    while(step--){ if(cross(u[(i+1)%n]-u[i],v[(j+1)%m]-v[j])>=0) j=(j+1)%m; else i=(i+1)%n; ans=max(ans,sq(u[i]-v[j])); }
    return sqrt(ans);
}
// ==================== CIRCLES ====================
struct circle {
    pt p; T r;
    circle() {}
    circle(pt p, T r) : p(p), r(r) {}
};

bool interCC(circle C1, circle C2, pt &p1, pt &p2) {
    T d = abs(C1.p - C2.p);
    if (d > C1.r + C2.r + EPS || d + min(C1.r, C2.r) + EPS < max(C1.r, C2.r)) return false;
    T a = (sq(C1.r) - sq(C2.r) + sq(d)) / (2*d);
    T h = sqrt(max((T)0, sq(C1.r) - sq(a)));
    pt P0 = C1.p + (C2.p - C1.p) * (a/d);
    pt offset = perp(C2.p - C1.p) * (h/d);
    p1 = P0 + offset;
    p2 = P0 - offset;
    return true;
}

bool interCL(circle C, line L, pt &p1, pt &p2) {
    pt P = L.proj(C.p);
    T h2 = sq(C.r) - sq(P - C.p);
    if(h2 < -EPS) return false;
    T h = sqrt(max((T)0, h2));
    pt dir = L.v / abs(L.v) * h;
    p1 = P + dir; 
    p2 = P - dir;
    return true;
}

vector<pt> tangents(circle C1, circle C2) {
    vector<pt> res;
    for(int sign1 = +1; sign1 >= -1; sign1 -= 2) {
        T r = C2.r * sign1 - C1.r;
        pt c2c1 = C2.p - C1.p;
        T z = sq(c2c1) - r*r;
        if(z < -EPS) continue;
        z = sqrt(max((T)0, z));
        pt v = (c2c1*r + pt{-c2c1.y, c2c1.x}*z)/sq(c2c1);
        res.push_back(C1.p + C1.r*v);
    }
    return res;
}
