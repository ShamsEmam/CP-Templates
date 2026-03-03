#define ld long double

typedef ld T;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
typedef complex<T> pt;
#define x real()
#define y imag()
T const EPS = 1e-12;

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
