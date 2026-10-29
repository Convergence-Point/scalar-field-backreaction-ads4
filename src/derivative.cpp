#include "derivative.hpp"

// 4th order centred first derivative at point p.
// Ghost points beyond either end are filled by reflection:  v(-k) = leftParity * v(k),
// v(N+k) = rightParity * v(N-k)   (+1 even, -1 odd).
// 4th order centred stencil applied to five values f = (f(i-2), f(i-1), f(i), f(i+1), f(i+2)).
static double stencil5(const array<double, 5> &f, double dx){
    const array<double, 5> factors {1./12., -2./3., 0, 2./3., -1./12.};
    return std::inner_product(f.begin(), f.end(), factors.begin(), 0.) / dx;
}

double derP(
    const vector<double> &v,
    int p,
    double dx,
    double leftParity,
    double rightParity
){
    int N = v.size() - 1;

    // fixed-size arrays: no heap allocation on this hot path (it is called ~10^6 times per run)
    const array<double, 5> f {
        (p < 2? leftParity : 1.) * v.at(abs(p - 2)),
        (p < 1? leftParity : 1.) * v.at(abs(p - 1)),
        v.at(p),
        (p > N-1? rightParity : 1.) * v.at(N-abs(N - p - 1)),
        (p > N-2? rightParity : 1.) * v.at(N-abs(N - p - 2))
    };
    return stencil5(f, dx);
}

double der(
    const vector<double> &v,
    int p,
    double dx,
    Field field
){
    double parity = (field == PhiField ? -1. : 1.);
    return derP(v, p, dx, parity, parity);
}

// AdS4 boundary behaviour: phi = f(t) rho^3 + O(rho^5), rho = pi/2 - x, so
//   Pi  : even at x = 0, odd  at x = pi/2
//   Phi : odd  at x = 0, even at x = pi/2   (Phi ~ -3 f rho^2)
//   g = tan^2 x Phi : odd at 0, even at pi/2 with finite limit g(pi/2) = -3 f
void adsDPi(const vector<double> &Pi, double dr, vector<double> &out){
    out.resize(Pi.size());
    for (size_t i = 0; i < Pi.size(); i++) out[i] = derP(Pi, i, dr, 1., -1.);
}

// Right-hand side of the Pi equation on AdS4:
//
//     out(x) = (1/tan^2 x) d(tan^2 x Phi)/dx  =  Phi_x + 2 Phi / (sin x cos x),
//
// with tan^2 x = sin^2 x / cos^2 x.
//
// Input : Phi[i] at x_i = i*dr, i = 0..N (x_N = pi/2)
// Output: out[i] = Pi_t at every grid point
//
// Why this form: the expanded version Phi_x + 2 Phi/(sin x cos x) (see adsDivergenceExpanded) is
// unstable in the evolution for any time step: it starts to grow near the outer boundary and,
// if that is cured, near the origin; only the middle of the grid can use it. Differentiating the
// product g = (sin^2 x / cos^2 x) Phi with the same 4th order stencil and multiplying by
// cos^2 x / sin^2 x afterwards (conservative form) is stable.
//
// NOTE: this conservative operator has a growing mode at the first grid points (Re lambda about 0.2, independent of n;
// roundoff grows into a blow-up after t of order 100-300) because it is not skew in the energy weight tan^2 x. It is
// kept for test_convergence_linear.cpp (T <= pi). For long runs use adsDivergenceExpanded / adsDPiExpanded, or ads/einstein, which
// use the symmetrized right-hand side of einstein.cpp.
//
// Parities (see the comment above adsDPi):
//   Phi is odd at x = 0 and even at x = pi/2  (Phi ~ -3 f rho^2),
//   sin^2 x / cos^2 x is even about both ends, so g has the same parities as Phi,
//   and g(pi/2) is finite and generally nonzero (g -> -3 f).
void adsDivergence(const vector<double> &Phi, double dr, vector<double> &out){
    int N = Phi.size() - 1;
    vector<double> g(N + 1);
    out.resize(N + 1);

    // Step 1: g = (sin^2 x / cos^2 x) * Phi at every point except the last one, where
    // cos(pi/2) = 0 (the finite limit is filled in by step 2).
    for (int i = 0; i < N; i++){
        double s = std::sin(dr * i);
        double c = std::cos(dr * i);
        g[i] = (s * s) / (c * c) * Phi[i];
    }

    // Step 2: g(pi/2). g is even about pi/2, so with rho = pi/2 - x
    //     g = g0 + g2 rho^2 + g4 rho^4 + ...
    // and g0 = g[N] is unknown (Phi[N] = 0 but sin^2/cos^2 = infinity, so 0 * inf).
    // Fit a polynomial in rho^2 through g at rho = dr, 2dr, 3dr (rho^2 = 1, 4, 9 in
    // units of dr^2) and evaluate it at rho = 0. The Lagrange weights are
    //     (4*9)/((4-1)(9-1)) = 1.5,  (1*9)/((1-4)(9-4)) = -0.6,  (1*4)/((1-9)(4-9)) = 0.1
    // (they sum to 1). The truncation error is O(rho^6) = O(dr^6), smaller than the
    // O(dr^4) error of the stencil.
    g[N] = 1.5 * g[N-1] - 0.6 * g[N-2] + 0.1 * g[N-3];

    // Step 3: apply the operator point by point.
    out[0] = 3. * derP(Phi, 0, dr, -1., 1.);
    out[N] = 0.;
    double s = std::sin(dr);
    double c = std::cos(dr);
    out[1] = derP(Phi, 1, dr, -1., 1.) + 2. * Phi[1] / (s * c);
    for (int i = 2; i < N; i++){
        // Interior: stencil on g, then multiply by cos^2 x / sin^2 x.
        // Ghost points beyond the ends are reflections of g: odd at x = 0 (-1),
        // even at x = pi/2 (+1).
        double s = std::sin(dr * i);
        double c = std::cos(dr * i);
        out[i] = derP(g, i, dr, -1., 1.) * (c * c) / (s * s);
    }
}

// ---------------------------------------------------------------------------------------------
// "Expanded" form of the same equations, made stable by symmetrizing them.
//
// The equations are  Pi_t = Phi_x + 2 Phi/(sin x cos x)  and  Phi_t = Pi_x.
// Writing the first one directly with a centred stencil, Phi_x + 2 Phi/(sin x cos x), is unstable
// (growth rate ~ 0.39/dr near x = pi/2, ~ 0.1/dr near x = 0, so a finer grid blows up sooner and no
// time step helps). The reason is an energy argument. The continuous system conserves
//     E = int tan^2 x (Pi^2 + Phi^2) dx
// only because  Pi (w Phi)_x + w Phi Pi_x = (w Pi Phi)_x  with w = tan^2 x, which uses the product
// rule  (w Phi)_x = w Phi_x + w_x Phi. A finite-difference stencil does not obey the product rule:
// for smooth data the defect is tiny (1e-5 or less), but for grid-scale oscillations it is O(1)
// (0.6 to 2.2 relative), which feeds energy into exactly the modes that grow.
//
// Cure: absorb tan x into the fields, Theta = tan x Pi, Psi = tan x Phi. Then
//     Theta_t = Psi_x + b Psi,     Psi_t = Theta_x - b Theta,     b = 1/(sin x cos x),
// where the b terms enter with opposite signs, so in  Theta Psi_t + Psi Theta_t  they cancel
// pointwise, with no product rule. With the reflection rules used here the two stencils satisfy
// summation by parts (H D_theta + D_psi^T H = 0, H = diag(1/2, 1, ..., 1, 1/2)), so
//     E = sum H (Theta^2 + Psi^2)
// is conserved exactly by the semi-discrete system: its eigenvalues are purely imaginary
// (largest real part 1e-14 for n = 50..400, |lambda| dr = 1.44, well inside the RK4 limit).
// Going back to Pi and Phi is a diagonal change of variables (the same eigenvalues):
//     Pi_t  = D(tan x Phi)/tan x + Phi/(sin x cos x)      -> adsDivergenceExpanded
//     Phi_t = D(tan x Pi )/tan x - Pi /(sin x cos x)      -> adsDPiExpanded
// Both are needed: the cancellation only works if the two equations are written as a pair.
//
// Parities (from phi = f rho^3 + O(rho^5) at pi/2):
//   Theta: odd at 0 (x * even), even at pi/2 (~ rho^2)      Psi: even at 0 (x * odd), odd at pi/2 (~ rho)
//   Both vanish at x = 0 and x = pi/2, so those two nodes are held at zero in the pair (Theta, Psi).
// Compared with the conservative form this divides by tan x instead of tan^2 x, so the error
// amplification next to the origin is 1/dr instead of 1/dr^2. With a 5-point stencil that left O(dr^3)
// at the first grid points; the 7-point stencil of weightedStencil makes it O(dr^5) there, so the
// operators are at least 4th order at every grid point.
// Theta and Psi are only used to derive the operators: the code below never forms them as arrays
// (see weightedStencil).
// ---------------------------------------------------------------------------------------------

// (1/tan x_i) * D(tan x * f)(i): the derivative of tan x * f divided by tan x, WITHOUT building the
// product as a separate array. tan x * f is evaluated on the fly at the stencil points, i.e. the
// ordinary stencil with each neighbour weighted by tan(x_{i+k})/tan(x_i). It equals
//     f_x + f/(sin x cos x)        (since (tan x)'/tan x = 1/(sin x cos x))
// in the continuum, but the product rule is built into the weights, which is what makes the pair
// of operators below energy-conserving (see the long comment above).
//
// Stencil: 7 points (6th order, weights -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60), not the 5-point one
// of derP. The division by tan x_i ~ i dr turns the stencil error dr^p f^(p+1) into dr^(p-1)/i next
// to the origin, so with 5 points (p = 4) the operator is only O(dr^3) at the first grid points;
// with 7 points (p = 6) it is O(dr^5) there and O(dr^6) elsewhere. Energy conservation is not lost:
// it only needs the stencil to be antisymmetric and the ghost points to follow the reflection rules,
// which holds for any width (a wider stencil only near the origin would break it).
//
// Ghost points use the parity of the product tan x * f. tan x is odd about x = 0 and about x = pi/2,
// so the product has parity  -leftParity  at 0  and  -rightParity  at pi/2  (f = Phi: leftParity = -1,
// rightParity = +1; f = Pi: +1 and -1). The product is zero at the two end nodes (tan 0 = 0, and
// f = O(rho^2) against tan = O(1/rho) at pi/2 gives a product that vanishes there as well).
static double weightedStencil(
    const vector<double> &f,
    int i,
    double dr,
    double leftParity,
    double rightParity
){
    int N = f.size() - 1;
    auto product = [&](int j){              // tan(x_j) * f_j, for j in [-3, N+3]
        double sign = 1.;
        int m = j;
        if (j < 0)      { m = -j;        sign = -leftParity;  }   // reflect through x = 0
        else if (j > N) { m = 2 * N - j; sign = -rightParity; }   // reflect through x = pi/2
        if (m == 0 || m == N) return 0.;
        return sign * std::tan(dr * m) * f[m];
    };
    // antisymmetric centred weights c_k for the offsets +-k (c_0 = 0)
    static const array<double, 4> c {0., 3./4., -3./20., 1./60.};
    double sum = 0.;
    for (int k = 1; k <= 3; k++) sum += c[k] * (product(i + k) - product(i - k));
    return sum / dr / std::tan(dr * i);
}

// Pi_t = (1/tan x) D(tan x Phi) + Phi/(sin x cos x)        [ = Phi_x + 2 Phi/(sin x cos x) ]
// Ends: Pi_t(0) = 3 Phi_x (value of Pi at the origin; it does not enter any other equation of the
// pair, since Theta(0) = 0), Pi_t(pi/2) = 0 (Pi odd).
void adsDivergenceExpanded(const vector<double> &Phi, double dr, vector<double> &out){
    int N = Phi.size() - 1;
    out.assign(N + 1, 0.);

    out[0] = 3. * derP(Phi, 0, dr, -1., 1.);
    for (int i = 1; i < N; i++){
        double s = std::sin(dr * i);
        double c = std::cos(dr * i);
        out[i] = weightedStencil(Phi, i, dr, -1., 1.) + Phi[i] / (s * c);
    }
    out[N] = 0.;
}

// Phi_t = (1/tan x) D(tan x Pi) - Pi/(sin x cos x)        [ = Pi_x ]
// Pi has parity +1 at x = 0 and -1 at x = pi/2, so tan x * Pi has the opposite parity at each end
// (see weightedStencil). Ends: Phi_t = 0 at both nodes (Phi is odd at 0 and Phi ~ rho^2 at pi/2).
void adsDPiExpanded(const vector<double> &Pi, double dr, vector<double> &out){
    int N = Pi.size() - 1;
    out.assign(N + 1, 0.);

    for (int i = 1; i < N; i++){
        double s = std::sin(dr * i);
        double c = std::cos(dr * i);
        out[i] = weightedStencil(Pi, i, dr, 1., -1.) - Pi[i] / (s * c);
    }
}
