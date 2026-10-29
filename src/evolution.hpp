#ifndef EVOLUTION_H
#define EVOLUTION_H
#include <vector>
#include <numeric>
#include "structure.hpp"

struct RK4 {
    std::vector<double> k1, k2, k3, k4;
};

// Optional operator applied to the whole Phi vector; its result at point i is
// passed to the callbacks as "dPhidx" (default: plain d/dx).
typedef void (*PhiOperator)(const std::vector<double> &Phi, double dr, std::vector<double> &out);

// Generic right-hand side: given Pi and Phi on the grid, return Pi_t and Phi_t.
typedef void (*RHSFunction)(
    const std::vector<double> &Pi,
    const std::vector<double> &Phi,
    double dr,
    std::vector<double> &kPi,
    std::vector<double> &kPhi
);

// One RK4 step of size sf.dt for a coupled system defined by an RHSFunction
// (used when Pi_t and Phi_t depend on both fields, e.g. einsteinRHS).
ScalarField evolveRK4(ScalarField sf, RHSFunction rhs);

typedef PhiOperator PiOperator; // same signature: whole-vector operator on Pi, result passed as "dPidx"

ScalarField evolve(
    ScalarField sf,
    double (*dPi)(double, double, double, double, double, double),
    double (*dPhi)(double, double, double, double, double, double),
    PhiOperator phiOp = nullptr,
    PiOperator piOp = nullptr
);

#endif