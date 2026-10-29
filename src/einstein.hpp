#ifndef EINSTEIN_HPP
#define EINSTEIN_HPP

#include <vector>

// Einstein-scalar system on global AdS4 (spherical symmetry, x in [0, pi/2], r = L tan x).
// Metric: ds^2 = L^2/cos^2 x ( -A e^{-2 delta} dt^2 + dx^2/A + sin^2 x dOmega^2 )
// Variables: Phi = phi_x, Pi = A^{-1} e^{delta} phi_t
//
//   Phi_t = (A e^{-delta} Pi)'
//   Pi_t  = (1/tan^2 x) (tan^2 x A e^{-delta} Phi)'
//   delta' = - sin x cos x (Phi^2 + Pi^2),                          delta(0) = 0
//   A'     = (1 + 2 sin^2 x)/(sin x cos x) (1 - A) - sin x cos x A (Phi^2 + Pi^2),   A(0) = 1

// Solves the two constraints for A and delta given Pi and Phi (see einstein.cpp for the
// closed integral forms that are used).
void solveGeometry(
    const std::vector<double> &Pi,
    const std::vector<double> &Phi,
    double dr,
    std::vector<double> &A,
    std::vector<double> &delta
);

// Minimum of A(x) for the given fields (solves the constraints). If `where` is not null it
// receives the grid index of the minimum. Returns NaN if A is not finite anywhere, so a
// collapse test written as  !(minimumA(...) >= threshold)  also catches a blown-up run.
double minimumA(
    const std::vector<double> &Pi,
    const std::vector<double> &Phi,
    double dr,
    int *where = nullptr
);

// Right-hand sides Pi_t and Phi_t (use with evolveRK4).
void einsteinRHS(
    const std::vector<double> &Pi,
    const std::vector<double> &Phi,
    double dr,
    std::vector<double> &kPi,
    std::vector<double> &kPhi
);

// Same right-hand sides on the fixed AdS background A = 1, delta = 0 (linear wave equation).
void adsRHS(
    const std::vector<double> &Pi,
    const std::vector<double> &Phi,
    double dr,
    std::vector<double> &kPi,
    std::vector<double> &kPhi
);

// Linear wave equation on the fixed AdS background in the symmetrized "expanded" form
// (adsDivergenceExpanded + adsDPiExpanded). Equivalent to adsRHS and stable; see derivative.cpp.
void adsExpandedRHS(
    const std::vector<double> &Pi,
    const std::vector<double> &Phi,
    double dr,
    std::vector<double> &kPi,
    std::vector<double> &kPhi
);

#endif
