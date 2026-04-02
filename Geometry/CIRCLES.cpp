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
pair<double, point> findCircle(point a, point b, point c) {
	//create median, vector, its prependicular
	point m1 = (b+a)*0.5, v1 = b-a, pv1 = point(v1.Y, -v1.X);
	point m2 = (b+c)*0.5, v2 = b-c, pv2 = point(v2.Y, -v2.X);
	point end1 = m1+pv1, end2 = m2+pv2, center;
	intersectSegments(m1, end1, m2, end2, center);
	return make_pair( length(a-center), center );  
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
