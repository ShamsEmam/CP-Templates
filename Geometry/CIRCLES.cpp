#define ld long double

typedef ld T;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
typedef complex<T> pt;
#define x real()
#define y imag()
T const EPS = 1e-12;


//////////////////////////////////////////  CIRCLES   //////////////////////////////////////////
struct Circle {
    pt center;
    ld radius;
};

// Circumcircle of 3 points
Circle findCircle(pt a, pt b, pt c) {
    pt mid1 = (a + b) / (T)2;
    pt mid2 = (b + c) / (T)2;
    pt d1 = pt(-(b - a).y, (b - a).x); // perpendicular to (b-a)
    pt d2 = pt(-(c - b).y, (c - b).x); // perpendicular to (c-b)

    // mid1 + t*d1 = mid2 + s*d2  => solve for t
    // cross both sides with d2
    T denom = cross(d1, d2);
    if (fabs(denom) < EPS) {
        // collinear: return circle around the two farthest points
        ld d1_ = len(a - b), d2_ = len(b - c), d3_ = len(a - c);
        if (d1_ >= d2_ && d1_ >= d3_) return { (a+b)/(T)2, d1_/(T)2 };
        if (d2_ >= d1_ && d2_ >= d3_) return { (b+c)/(T)2, d2_/(T)2 };
        return { (a+c)/(T)2, d3_/(T)2 };
    }

    T t = cross(mid2 - mid1, d2) / denom;
    pt center = mid1 + d1 * t;
    return { center, len(a - center) };
}

// ==================== WELZL ====================
const int MAX = 305;
pt pnts[MAX], r[3], cen;
ld rad;
int ps, rs;

void MEC() {
    if (ps == 0 && rs == 0) { rad = 0; return; }
    if (ps == 0 && rs == 1) { cen = r[0]; rad = 0; return; }
    if (ps == 0 && rs == 2) {
        cen = (r[0] + r[1]) / (T)2;
        rad = len(r[0] - cen);
        return;
    }
    if (rs == 3) {
        Circle c = findCircle(r[0], r[1], r[2]);
        cen = c.center;
        rad = c.radius;
        return;
    }

    ps--;
    MEC();
    if (len(pnts[ps] - cen) > rad + EPS) {
        r[rs++] = pnts[ps];
        MEC();
        rs--;
    }
    ps++;
}

Circle minCircle(vector<pt>& pts) {
    if (pts.empty()) return { {0,0}, 0 };
    ps = sz(pts); rs = 0;
    for (int i = 0; i < ps; i++) pnts[i] = pts[i];
    shuffle(pnts, pnts + ps, mt19937{random_device{}()});
    MEC();
    return { cen, rad };
}
pair<pt, T> circumCircle(pt a, pt b, pt c) {
    b = b - a, c = c - a; // consider coordinates relative to A
    assert(cross(b,c) != 0); // no circumcircle if A,B,C aligned
    return {
        a + perp(b * sq(c) - c * sq(b)) / cross(b, c) / (T) 2, abs(perp(b * sq(c) - c * sq(b)) / cross(b, c) / (T) 2)
    };
}

int circleLine(pt o, double r, line l, pair<pt, pt> &out) {
    double h2 = r * r - l.sqDist(o);
    if (h2 >= 0) {
        // the line touches the circle
        pt p = l.proj(o); // point P
        pt h = l.v * (T) (sqrt(h2) / abs(l.v)); // vector parallel to l, of length h
        out = {p - h, p + h};
    }
    return 1 + sgn(h2);
}


int circleCircle(pt o1, T r1, pt o2, T r2, pair<pt, pt> &out) {
    pt d = o2 - o1;
    T d2 = sq(d);
    if (d2 == 0) {
        assert(r1 != r2);
        return 0;
    } // concentric circles
    T pd = (d2 + r1 * r1 - r2 * r2) / 2; // = |O_1P| * d
    T h2 = r1 * r1 - pd * pd / d2; // = hˆ2
    if (h2 >= 0) {
        pt p = o1 + d * pd / d2, h = perp(d) * sqrt(h2 / d2);
        out = {p - h, p + h};
    }
    return 1 + sgn(h2);
}

int tangents(pt o1, T r1, pt o2, T r2, bool inner, vector<pair<pt, pt> > &out) {
    if (inner) r2 = -r2;
    pt d = o2 - o1;
    T dr = r1 - r2, d2 = sq(d), h2 = d2 - dr * dr;
    if (d2 == 0 || h2 < 0) {
        assert(h2 != 0);
        return 0;
    }
    for (T sign: {-1, 1}) {
        pt v = (d * dr + perp(d) * sqrt(h2) * sign) / d2;
        out.push_back({o1 + v * r1, o2 + v * r2});
    }
    return 1 + (h2 > 0);
}

//B not included
/// no of point on line
int laticeAB(pt A, pt B) {
    return __gcd(abs(A.x - B.x), abs(A.y - B.y));
}

// the points on polygon not out not inside
int calcBoundary(vector<pt> &p) {
    int n = p.size(), ans = 0;
    for (int i = 0; i < n; ++i) {
        ans += laticeAB(p[i], p[(i + 1) % n]);
    }
    return ans;
}

T circleIntersectionArea(pt c1, T r1, pt c2, T r2) {
    T d = abs(c1 - c2); // distance between centers

    // No intersection or external tangent
    if (d > r1 + r2 - EPS) return 0.0;

    // One circle is completely inside the other (or internal tangent)
    if (d < abs(r1 - r2) + EPS) {
        T r = min(r1, r2);
        return M_PI * r * r;
    }

    // Clamp value to avoid precision issues with acos
    auto clamp = [&](T val) {
        return max((T)-1, min((T)1, val));
    };

    // Angles of the circular segments
    T alpha = acos(clamp((r1*r1 + d*d - r2*r2) / (2*r1*d)));
    T beta  = acos(clamp((r2*r2 + d*d - r1*r1) / (2*r2*d)));

    // Area of intersection = sum of sector areas - triangle areas
    T area =
        r1*r1*alpha +
        r2*r2*beta -
        r1*r1*sin(2*alpha)/2 -
        r2*r2*sin(2*beta)/2;

    return area;
}
