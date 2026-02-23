typedef ld T;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
typedef complex<T> pt;
#define x real()
#define y imag()
T const EPS = 1e-12;


// r^2
T sq(pt p) {
    return p.y * p.y + p.x * p.x;
}

T dot(pt v, pt w) {
    return v.x * w.x + v.y * w.y;
}

T cross(pt v, pt w) {
    return v.x * w.y - v.y * w.x;
}

// Precisions and Errors
int sgn(T val) {
    if (val > EPS) return 1;
    if (val < -EPS) return -1;
    return 0;
}

// cos(90)=0
bool isPerp(pt v, pt w) {
    return fabs(dot(v, w)) < EPS;
}

// to get the prep vector
pt prep(pt p) {
    return {-p.y, p.x};
}

//-------------------- Transformation -----------------------

pt translate(pt v, pt p) {
    return p + v;
}

//  factor (scaling value)
pt scale(pt c, ld factor, pt p) {
    return c + (p - c) * factor;
}

// angle in red
// rotate p around c by angle a
pt rot(pt p, pt c, ld a) {
    pt v = p - c;
    pt rotate = {cos(a), sin(a)};
    return c + rotate * v;
}

// ---general transformation
pt linearTransfo(pt p, pt q, pt r, pt fp, pt fq) {
    return fp + (r - p) * (fq - fp) / (q - p);
}

//------------------- angle -------------------

// make line btw 2 point and check if the other point at right or left
T orient(pt a, pt b, pt c) {
    return cross(b - a, c - a);
}

// calc angle btw two vectors from org point
// V . W =|V| * |W| * cos(sita)
T angle(pt v, pt w) {
    return acos(clamp(dot(v, w) / abs(v) / abs(w), (T) -1.0, (T) 1.0));
}

// calc right angle from user point
// if is acute or obtuse as we need
T orientedAngle(pt a, pt b, pt c) {
    ld ampli = angle(b - a, c - a);
    if (orient(a, b, c) > 0) return ampli;
    else return 2 * M_PI - ampli;
}

// small angle
// always acute angle
T angleTravelled(pt a, pt b, pt c) {
    ld ampli = angle(b - a, c - a);
    if (orient(a, b, c) > 0) return ampli;
    else return -ampli;
}

//check p in between angle(bac) counter clockwise
bool inAngle(pt a, pt b, pt c, pt p) {
    T abp = orient(a, b, p), acp = orient(a, c, p), abc = orient(a, b, c);
    if (abc < 0) swap(abp, acp);
    return (abp >= 0 && acp <= 0) ^ (abc < 0);
}

struct line {
    pt v;
    T c;

    line(pt v, T c) : v(v), c(c) {}

    // from equation ax+by = c
    line(T a, T b, T _c) {
        // parallel vector
        v = {b, -a};
        c = _c;
    }

    //line from two points
    line(pt p, pt q) {
        // parallel vector btw two points
        v = q - p;
        // v * (x , y) = c -->line equation
        c = cross(v, p);
    }

    // to check if the point lies on
    // (b,-a)-->parallel vector
    // (b,-a)*(x,y)=by+ax=c-->this is the line equation
    T side(pt p) { return cross(v, p) - c; }

    // distance btw point and line
    ld dist(pt p) { return abs(side(p)) / abs(v); }

    T distPointLine(T X, T Y, T A, T B, T C) {
        return abs(A * X + B * Y + C) / sqrt(A * A + B * B);
    }

    double sqDist(pt p) { return side(p) * side(p) / (T) sq(v); }

    //line throught point and prep on vector
    line prepThrought(pt p) { return {p, p + prep(p)}; }

    bool cmpProj(pt p, pt q) {
        return dot(v, p) < dot(v, q);
    }

    line translate(pt t) { return {v, c + cross(v, t)}; }

    line shiftLeft(T dist) { return {v, c + dist * abs(v)}; }

    pt proj(pt p) { return p - prep(v) * side(p) / sq(v); }

    pt refl(pt p) { return p - prep(v) * (T) 2.0 * side(p) / sq(v); }
};

//intersection btw two line
bool inter(line l1, line l2, pt &out) {
    T d = cross(l1.v, l2.v);
    if (fabs(d) < EPS) return false;
    out = (l2.v * l1.c - l1.v * l2.c) / d; // requires floating-point coordinates
    return true;
}

// line that divided the in angel btw two line
line bisector(line l1, line l2, bool interior) {
    assert(cross(l1.v, l2.v) != 0); // l1 and l2 cannot be parallel!
    ld sign = interior ? 1 : -1;
    return {l2.v / abs(l2.v) + l1.v / abs(l1.v) * sign,
            l2.c / abs(l2.v) + l1.c / abs(l1.v) * sign};
}

//////////////////////////////////////////  SEGMENTS   //////////////////////////////////////////


bool inDisk(pt a, pt b, pt p) {
    return dot(a - p, b - p) <= EPS;
}

bool onSegment(pt a, pt b, pt c) {
    return orient(a, b, c) == 0 && inDisk(a, b, c);
}

bool properInter(pt a, pt b, pt c, pt d, pt &out) {
    T oa = orient(c, d, a),
            ob = orient(c, d, b),
            oc = orient(a, b, c),
            od = orient(a, b, d);
// Proper intersection exists iff opposite signs
    if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0) {
        out = (a * ob - b * oa) / (ob - oa);
        return true;
    }
    return false;
}

set<pair<ld, ld>> inters(pt a, pt b, pt c, pt d) {
    set<pair<ld, ld>> s;
    pt out;
    if (a == c || a == d) {
        s.insert(make_pair(a.x, a.y));
    }
    if (b == c || b == d) {
        s.insert(make_pair(b.x, b.y));
    }
    if (s.size()) return s;

    if (properInter(a, b, c, d, out)) return {make_pair(out.x, out.y)};
    if (onSegment(c, d, a)) s.insert(make_pair(a.x, a.y));
    if (onSegment(c, d, b)) s.insert(make_pair(b.x, b.y));
    if (onSegment(a, b, c)) s.insert(make_pair(c.x, c.y));
    if (onSegment(a, b, d)) s.insert(make_pair(d.x, d.y));

    return s;
}

ld segPoint(pt a, pt b, pt p) {
    if (a != b) {
        line l(a, b);
        if (l.cmpProj(a, p) && l.cmpProj(p, b)) // if closest to projection
            return l.dist(p); // output distance to line
    }
    return min(abs(p - a), abs(p - b)); // otherwise distance to A or B
}

ld segSeg(pt a, pt b, pt c, pt d) {
    pt dummy;
    if (properInter(a, b, c, d, dummy))
        return 0;
    return min({segPoint(a, b, c), segPoint(a, b, d),
                segPoint(c, d, a), segPoint(c, d, b)});
}

ld areaTriangle(pt a, pt b, pt c) {
    return abs(cross(b - a, c - a)) / 2.0;
}

ld areaPolygon(vector<pt> p) {
    ld area = 0.0;
    for (int i = 0, n = sz(p); i < n; i++) {
        area += cross(p[i], p[(i + 1) % n]); // wrap back to 0 if i == n - 1
    }
    return abs(area) / 2.0;
}

bool above(pt a, pt p) {
    return p.y >= a.y;
}

// check if [PQ] crosses ray from A
bool crossesRay(pt a, pt p, pt q) {
    return (above(a, q) - above(a, p)) * orient(a, p, q) > 0;
}

bool inPolygon(vector<pt> p, pt a, bool strict = true) {
    int numCrossings = 0;
    for (int i = 0, n = p.size(); i < n; i++) {
        if (onSegment(p[i], p[(i + 1) % n], a))
            return !strict;
        numCrossings += crossesRay(a, p[i], p[(i + 1) % n]);
    }
    return numCrossings & 1; // inside if odd number of crossings
}

ld rayPoint(pt A, pt B, pt P) {
    pt AB = B - A;
    pt AP = P - A;

    if (dot(AP, AB) >= 0) {
        T t = dot(AP, AB) / sq(AB);
        pt Q = A + AB * t;
        return abs(P - Q);
    } else {
        return abs(P - A);
    }
}
ld rayRay(pt A, pt B, pt C, pt D) {
    pt AB = B - A;
    pt CD = D - C;

    // check if they are on same line
    if (fabs(cross(AB, CD)) < EPS && fabs(cross(AB, C - A)) < EPS) {
        // نفس الخط
        if (dot(C - A, AB) >= 0 || dot(A - C, CD) >= 0) {
            return 0;
        }
    }

    line l1(A, B);
    line l2(C, D);

    pt interPoint;
    if (inter(l1, l2, interPoint)) {
        if (dot(interPoint - A, AB) >= 0 &&
            dot(interPoint - C, CD) >= 0) {
            return 0;
        }
    }

    return min({
                       rayPoint(A, B, C),
                       rayPoint(A, B, D),
                       rayPoint(C, D, A),
                       rayPoint(C, D, B)
               });
}
void takePoint(pt &p) {
    T _, __;
    cin >> _ >> __;
    p = {_, __};
}
