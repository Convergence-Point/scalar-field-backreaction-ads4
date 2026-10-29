// Fixed-time self-convergence test for the RK4 + 4th order spatial derivative.
// Pi_t = Phi_x, Phi_t = Pi_x. Every resolution is sampled at the same times.
//
// build: make build/tests/test_convergence_linear
// output (in ./convergence/):
//   error.dat : t  dr  ||Pi_n - Pi_2n||_2   (one block per resolution pair)
//   order.dat : t  order
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include "derivative.hpp"
#include "structure.hpp"
#include "evolution.hpp"

static double dPhidt(double, double, double, double, double dPidx, double) { return dPidx; }
static double dPidt(double, double, double, double, double, double dPhidx) { return dPhidx; }

// snapshots[k] = Pi at time t_k = k * dT
static std::vector<std::vector<double>> run(int n, int nSamples, int stepsPerSample)
{
    ScalarField sf;
    sf.t = 0.;
    sf.dr = M_PI_2 / n;
    sf.dt = sf.dr / 2.;
    sf.r.resize(n + 1);
    sf.Phi.assign(n + 1, 0.);
    sf.Pi.resize(n + 1);
    for (int i = 0; i <= n; i++) {
        sf.r[i] = sf.dr * i;
        sf.Pi[i] = gaussTan(sf.r[i], 1.) / M_PI;
    }

    std::vector<std::vector<double>> snaps;
    snaps.push_back(sf.Pi);
    for (int k = 1; k <= nSamples; k++) {
        for (int s = 0; s < stepsPerSample; s++) sf = evolve(sf, &dPidt, &dPhidt, &adsDivergence, &adsDPi);
        snaps.push_back(sf.Pi);
    }
    return snaps;
}

// L2 norm of the difference on the coarse grid points
static double diffNorm(const std::vector<double> &coarse, const std::vector<double> &fine, double dr)
{
    double sum = 0.;
    for (size_t i = 0; i < coarse.size(); i++) {
        double d = coarse[i] - fine[2 * i];
        sum += d * d * dr;
    }
    return std::sqrt(sum);
}

int main()
{
    const std::vector<int> ns = {250, 500, 1000, 2000};
    const int nSamples = 20;       // t_k = k * pi/20, up to t = pi
    // dt = pi/(4n)  ->  pi/20 corresponds to n/5 steps
    std::vector<std::vector<std::vector<double>>> data;
    for (int n : ns) {
        std::cerr << "running n = " << n << std::endl;
        data.push_back(run(n, nSamples, n / 5));
    }

    mkdir("convergence", 0755);
    std::ofstream err("convergence/error.dat"), ord("convergence/order.dat");
    err << std::setprecision(15);
    ord << std::setprecision(15);
    err << "# t  dr(coarse)  L2 error\n";
    ord << "# t  order\n";

    std::vector<std::vector<double>> E(ns.size() - 1, std::vector<double>(nSamples + 1));
    for (size_t j = 0; j + 1 < ns.size(); j++) {
        double dr = M_PI_2 / ns[j];
        for (int k = 0; k <= nSamples; k++) {
            double t = k * M_PI / 20.;
            E[j][k] = diffNorm(data[j][k], data[j + 1][k], dr);
            err << t << "\t" << dr << "\t" << E[j][k] << "\n";
        }
        err << "\n\n"; // gnuplot index separator
    }
    for (int k = 1; k <= nSamples; k++) {
        double t = k * M_PI / 20.;
        ord << t;
        for (size_t j = 0; j + 2 < ns.size(); j++)
            ord << "\t" << std::log2(E[j][k] / E[j + 1][k]);
        ord << "\n";
    }
    return 0;
}
