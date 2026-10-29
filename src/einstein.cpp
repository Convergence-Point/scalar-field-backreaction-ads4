#include "einstein.hpp"
#include "derivative.hpp"

#include <cmath>

using std::vector;

// Cumulative integral out[i] = int_0^{x_i} f dx on a uniform grid, 4th order accurate.
// Each interval [x_i, x_{i+1}] is integrated exactly for the cubic through four neighbouring
// points (one-sided at the first and last interval). Local error O(h^5), global O(h^4).
static void cumulativeIntegral(const vector<double> &f, double h, vector<double> &out)
{
    int N = f.size() - 1;
    out.assign(N + 1, 0.);
    for (int i = 0; i < N; i++) {
        double w;
        if (i == 0)
            w = 9. * f[0] + 19. * f[1] - 5. * f[2] + f[3];
        else if (i == N - 1)
            w = f[N - 3] - 5. * f[N - 2] + 19. * f[N - 1] + 9. * f[N];
        else
            w = -f[i - 1] + 13. * f[i] + 13. * f[i + 1] - f[i + 2];
        out[i + 1] = out[i] + h * w / 24.;
    }
}

// The constraints in closed form (delta(0) = 0, A(0) = 1):
//
//   delta(x) = - int_0^x sin y cos y S(y) dy,                S = Phi^2 + Pi^2
//
// The A equation is linear, A' + (Q + sin cos S) A = Q with Q = (1+2 sin^2)/(sin cos).
// Its integrating factor is mu = (sin x / cos^3 x) e^{-delta} (note exp(int sin cos S) = e^{-delta}),
// and after an integration by parts
//
//   1 - A(x) = e^{delta(x)} (cos^3 x / sin x) int_0^x e^{-delta(y)} tan^2 y S(y) dy.
//
// Both integrands are smooth, also at x = pi/2 where S ~ rho^4 and tan^2 ~ rho^-2.
// Computing 1 - A directly avoids the cancellation between A and 1 near x = pi/2 and the
// 1/(sin x cos x) singularities of the differential form.
void solveGeometry(
    const vector<double> &Pi,
    const vector<double> &Phi,
    double dr,
    vector<double> &A,
    vector<double> &delta)
{
    int N = Pi.size() - 1;
    vector<double> S(N + 1), f(N + 1), cum;

    for (int i = 0; i <= N; i++) S[i] = Phi[i] * Phi[i] + Pi[i] * Pi[i];

    // delta
    for (int i = 0; i <= N; i++) {
        double x = dr * i;
        f[i] = -std::sin(x) * std::cos(x) * S[i];
    }
    cumulativeIntegral(f, dr, delta);

    // I(x) = int_0^x e^{-delta} tan^2 S dy   (integrand -> 0 at x = pi/2)
    for (int i = 0; i < N; i++) {
        double t = std::tan(dr * i);
        f[i] = std::exp(-delta[i]) * t * t * S[i];
    }
    f[N] = 0.;
    cumulativeIntegral(f, dr, cum);

    A.assign(N + 1, 1.);
    for (int i = 1; i <= N; i++) {
        double x = dr * i;
        double c = std::cos(x);
        A[i] = 1. - std::exp(delta[i]) * c * c * c / std::sin(x) * cum[i];
    }
}

double minimumA(
    const vector<double> &Pi,
    const vector<double> &Phi,
    double dr,
    int *where)
{
    vector<double> A, delta;
    solveGeometry(Pi, Phi, dr, A, delta);
    double m = A[0];
    int im = 0;
    for (size_t i = 0; i < A.size(); i++) {
        if (!std::isfinite(A[i])) {
            if (where) *where = i;
            return std::nan("");
        }
        if (A[i] < m) { m = A[i]; im = i; }
    }
    if (where) *where = im;
    return m;
}

// Right-hand sides of
//     Phi_t = (m Pi)_x ,        Pi_t = (1/tan^2 x) (tan^2 x m Phi)_x ,        m = A e^{-delta}
// in a form that conserves a discrete energy exactly (no growing modes at the origin).
//
// Why not the conservative form (tan^2 x m Phi)_x / tan^2 x: its discrete operator cannot be skew with
// respect to the energy weight tan^2 x. The origin node Pi_0 has weight tan^2(0) = 0 but still drives
// Phi_1 and Phi_2 through the stencil, and its own equation (3 Phi_x) is not the conjugate of those
// couplings. The result is a growing mode at i = 1 that does not shrink with n (Re lambda about 0.2;
// roundoff grows into a blow-up after t of order 100-300, and for the Einstein system it looks like a
// collapse before it turns into NaN).
//
// Cure (as in derivative.cpp for m = 1): absorb the weight into the fields. With a = tan x, p = sqrt(m),
//     Theta = p a Pi ,     Psi = p a Phi ,     b = 1/(sin x cos x)
// the system becomes
//     Theta_t = S Psi + m b Psi ,        Psi_t = S Theta - m b Theta ,        S = (m d/dx + d/dx m) / 2
// S is skew-adjoint and the b terms cancel pointwise in Theta Psi_t + Psi Theta_t, so
// E = int (Theta^2 + Psi^2) is conserved: discretely, with the 7-point antisymmetric stencil and the
// reflection rules below, the matrix is skew in the reflection-extended inner product, for ANY profile m
// (checked: purely imaginary eigenvalues for m = 1, a bump, a steep dip and a profile with an odd rho^3
// part at pi/2). Going back to Pi and Phi is a diagonal change of variables:
//     Pi_t  = S(p a Phi) / (p a) + m b Phi ,        Phi_t = S(p a Pi) / (p a) - m b Pi .
// Theta and Psi are never stored: q_j = p_j a_j f_j is formed on the fly at the stencil points.
//
// Parities (phi = f rho^3 + O(rho^5), rho = pi/2 - x): the products q = p a f are
//   f = Phi: even at x = 0, odd at pi/2 ;    f = Pi: odd at x = 0, even at pi/2 ;   zero at both end nodes.
// m is even at x = 0 and is mirrored at pi/2 (its odd rho^3 part makes that mirror accurate only to
// O(rho^3); the division by p a ~ 1/rho brings it to O(dr^4) in the equation, as before).
// End values: Pi_t(0) = 3 (m Phi)_x(0) (Pi_0 does not enter any other equation), Pi_t(pi/2) = 0,
// Phi_t = 0 at both ends.
static void computeRHS(
    const vector<double> &Pi,
    const vector<double> &Phi,
    double dr,
    const vector<double> &A,
    const vector<double> &delta,
    vector<double> &kPi,
    vector<double> &kPhi)
{
    int N = Pi.size() - 1;

    vector<double> m(N + 1), p(N + 1), a(N + 1, 0.), H(N + 1);
    for (int i = 0; i <= N; i++) {
        m[i] = A[i] * std::exp(-delta[i]);
        p[i] = std::sqrt(m[i]);
        H[i] = m[i] * Phi[i];
        if (i > 0 && i < N) a[i] = std::tan(dr * i);   // a = 0 at both end nodes: the products vanish there
    }

    // 7-point antisymmetric weights (offsets +-k, k = 1..3)
    static const array<double, 4> c {0., 3./4., -3./20., 1./60.};

    // q_j = p_j a_j f_j for any j in [-3, N+3]; ghost points by the parity of the product
    auto q = [&](const vector<double> &f, int j, double leftParity, double rightParity){
        double sign = 1.;
        int mj = j;
        if (j < 0)      { mj = -j;        sign = -leftParity;  }
        else if (j > N) { mj = 2 * N - j; sign = -rightParity; }
        return sign * p[mj] * a[mj] * f[mj];             // a[0] = a[N] = 0
    };
    auto mext = [&](int j){ return m[j < 0 ? -j : (j > N ? 2 * N - j : j)]; };

    // S q at node i, divided by p_i a_i:  (m_i (D q)_i + (D (m q))_i) / (2 p_i a_i)
    auto sOverPa = [&](const vector<double> &f, int i, double leftParity, double rightParity){
        double dq = 0., dmq = 0.;
        for (int k = 1; k <= 3; k++) {
            double qp = q(f, i + k, leftParity, rightParity);
            double qm = q(f, i - k, leftParity, rightParity);
            dq  += c[k] * (qp - qm);
            dmq += c[k] * (mext(i + k) * qp - mext(i - k) * qm);
        }
        return 0.5 * (m[i] * dq + dmq) / dr / (p[i] * a[i]);
    };

    kPi.assign(N + 1, 0.);
    kPhi.assign(N + 1, 0.);
    kPi[0] = 3. * derP(H, 0, dr, -1., 1.);                // limit of (tan^2 H)'/tan^2 for H ~ x
    for (int i = 1; i < N; i++) {
        double bm = m[i] / (std::sin(dr * i) * std::cos(dr * i));      // m b
        kPi[i]  = sOverPa(Phi, i, -1., 1.) + bm * Phi[i];
        kPhi[i] = sOverPa(Pi,  i, 1., -1.) - bm * Pi[i];
    }
}

// Full system: A and delta follow from the constraints at every call.
void einsteinRHS(
    const vector<double> &Pi,
    const vector<double> &Phi,
    double dr,
    vector<double> &kPi,
    vector<double> &kPhi)
{
    vector<double> A, delta;
    solveGeometry(Pi, Phi, dr, A, delta);
    computeRHS(Pi, Phi, dr, A, delta, kPi, kPhi);
}

// Fixed AdS background (no backreaction): A = 1, delta = 0 for all x, i.e. the linear
// wave equation Phi_t = Pi_x, Pi_t = (1/tan^2 x)(tan^2 x Phi)'. Same code path as above.
void adsRHS(
    const vector<double> &Pi,
    const vector<double> &Phi,
    double dr,
    vector<double> &kPi,
    vector<double> &kPhi)
{
    vector<double> A(Pi.size(), 1.), delta(Pi.size(), 0.);
    computeRHS(Pi, Phi, dr, A, delta, kPi, kPhi);
}

// Fixed AdS background, equations in the symmetrized "expanded" form (see derivative.cpp):
//   Pi_t = D(tan x Phi)/tan x + Phi/(sin x cos x),  Phi_t = D(tan x Pi)/tan x - Pi/(sin x cos x).
void adsExpandedRHS(
    const vector<double> &Pi,
    const vector<double> &Phi,
    double dr,
    vector<double> &kPi,
    vector<double> &kPhi)
{
    adsDPiExpanded(Pi, dr, kPhi);
    adsDivergenceExpanded(Phi, dr, kPi);
}
