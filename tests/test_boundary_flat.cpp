// Boundary test against the exact solution of Pi_t = Phi_x, Phi_t = Pi_x on [0, pi/2]
// with Pi even / Phi odd at both ends (period-pi even extension):
//     Pi(x,t) = 1/2 [ f(x+t) + f(x-t) ],   Phi(x,t) = 1/2 [ f(x+t) - f(x-t) ]
// Initial data: Pi = f = exp(-(r-r0)^2/w^2), Phi = 0. The pulse splits in two,
// one half reflects at r = pi/2 (t ~ pi/4) and the other at r = 0.
//
// build: make build/tests/test_boundary_flat
// usage: ./build/tests/test_boundary_flat [width] [r0]         (defaults 0.1 and pi/4)
// output: convergence/boundary.dat  with columns  t  n  max|err| (all r)  max|err| (r > 1.2)
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h>
#include "derivative.hpp"
#include "structure.hpp"
#include "evolution.hpp"

static double W = 0.1, R0 = M_PI_4;

static double f(double x)
{
    double y = std::fmod(std::fabs(x), M_PI);   // even about 0, period pi
    if (y > M_PI_2) y = M_PI - y;               // even about pi/2
    return std::exp(-(y - R0) * (y - R0) / (W * W));
}

static double dPhidt(double, double, double, double, double dPidx, double) { return dPidx; }
static double dPidt(double, double, double, double, double, double dPhidx) { return dPhidx; }

int main(int argc, char **argv)
{
    if (argc > 1) W = std::atof(argv[1]);
    if (argc > 2) R0 = std::atof(argv[2]);

    const std::vector<int> ns = {250, 500, 1000, 2000};
    const int nSamples = 20; // t_k = k*pi/20
    mkdir("convergence", 0755);
    std::ofstream out("convergence/boundary.dat");
    out << std::setprecision(10);
    out << "# t  n  max_err_all  max_err_right(r>1.2)\n";

    std::vector<std::vector<double>> eAll(ns.size(), std::vector<double>(nSamples + 1));
    std::vector<std::vector<double>> eRight = eAll;

    for (size_t j = 0; j < ns.size(); j++) {
        int n = ns[j];
        ScalarField sf;
        sf.t = 0.;
        sf.dr = M_PI_2 / n;
        sf.dt = sf.dr / 2.;
        sf.r.resize(n + 1);
        sf.Phi.assign(n + 1, 0.);
        sf.Pi.resize(n + 1);
        for (int i = 0; i <= n; i++) {
            sf.r[i] = sf.dr * i;
            sf.Pi[i] = f(sf.r[i]);
        }
        for (int k = 0; k <= nSamples; k++) {
            if (k > 0)
                for (int s = 0; s < n / 5; s++) sf = evolve(sf, &dPidt, &dPhidt);
            double ea = 0., er = 0.;
            for (int i = 0; i <= n; i++) {
                double x = sf.r[i];
                double e = std::fabs(sf.Pi[i] - 0.5 * (f(x + sf.t) + f(x - sf.t)));
                e = std::max(e, std::fabs(sf.Phi[i] - 0.5 * (f(x + sf.t) - f(x - sf.t))));
                ea = std::max(ea, e);
                if (x > 1.2) er = std::max(er, e);
            }
            eAll[j][k] = ea;
            eRight[j][k] = er;
            out << sf.t << "\t" << n << "\t" << ea << "\t" << er << "\n";
        }
        out << "\n\n";
    }

    std::cout << std::setprecision(4);
    std::cout << "width=" << W << " r0=" << R0 << "\n";
    std::cout << "t\torder(all r)\torder(r>1.2)   [n=500/1000/2000 vs exact]\n";
    for (int k = 1; k <= nSamples; k++) {
        double oa = std::log2(eAll[2][k] / eAll[3][k]);
        double orr = std::log2(eRight[2][k] / eRight[3][k]);
        std::cout << k * M_PI / 20. << "\t" << oa << "\t\t" << orr << "\n";
    }
    return 0;
}
