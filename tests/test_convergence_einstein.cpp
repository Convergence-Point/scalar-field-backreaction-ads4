// Tests for the Einstein-scalar-AdS system (einstein.cpp).
//
// build: make build/tests/test_convergence_einstein
// usage: ./build/tests/test_convergence_einstein [amplitude] [sigma]       (defaults 0.5 and 1.0)
//
// 1. Constraint check: A and delta from solveGeometry are inserted into the differential
//    equations  delta' = -sin cos (Phi^2+Pi^2),
//               A' = (1+2 sin^2)/(sin cos) (1-A) - sin cos A (Phi^2+Pi^2)
//    and the residual is computed with the 4th order stencil (interior points). Expected order 4.
// 2. Self-convergence of the evolution (n = 250, 500, 1000, 2000), L2 norm of Pi, plus min A.
//    Output: convergence/einstein_order.dat  (t  order(250/500/1000)  order(500/1000/2000)  minA(n=2000))
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
#include "einstein.hpp"

static double amp = 0.5, sigma = 1.0;

static double piInit(double x) { return amp / M_PI * gaussTan(x, sigma); }
// test profile for Phi: odd at 0, even (~rho^2) at pi/2
static double phiTest(double x) { return 0.3 * x * std::cos(x) * std::cos(x); }

static void constraintCheck()
{
    std::cout << "== constraint residuals (max over interior points), expected order 4\n";
    std::cout << "n\tmax|delta' res|\tmax|A' res|\torder(delta)\torder(A)\n";
    double pd = 0., pa = 0.;
    for (int n = 250; n <= 2000; n *= 2) {
        double dr = M_PI_2 / n;
        std::vector<double> Pi(n + 1), Phi(n + 1), A, delta;
        for (int i = 0; i <= n; i++) { Pi[i] = piInit(dr * i); Phi[i] = phiTest(dr * i); }
        solveGeometry(Pi, Phi, dr, A, delta);
        double rd = 0., ra = 0.;
        for (int i = 3; i <= n - 3; i++) {
            double x = dr * i, s = std::sin(x), c = std::cos(x);
            double S = Phi[i] * Phi[i] + Pi[i] * Pi[i];
            double dd = derP(delta, i, dr, 1., 1.);
            double dA = derP(A, i, dr, 1., 1.);
            rd = std::max(rd, std::fabs(dd + s * c * S));
            ra = std::max(ra, std::fabs(dA - ((1. + 2. * s * s) / (s * c) * (1. - A[i]) - s * c * A[i] * S)));
        }
        std::cout << n << "\t" << rd << "\t" << ra;
        if (pd > 0.) std::cout << "\t" << std::log2(pd / rd) << "\t" << std::log2(pa / ra);
        std::cout << "\n";
        pd = rd; pa = ra;
    }
}

static std::vector<std::vector<std::vector<double>>> run(int n, int nSamples, int stepsPerSample, double &minA)
{
    ScalarField sf;
    sf.t = 0.;
    sf.dr = M_PI_2 / n;
    sf.dt = sf.dr / 2.;
    sf.r.resize(n + 1);
    sf.Phi.assign(n + 1, 0.);
    sf.Pi.resize(n + 1);
    for (int i = 0; i <= n; i++) { sf.r[i] = sf.dr * i; sf.Pi[i] = piInit(sf.r[i]); }

    std::vector<std::vector<std::vector<double>>> snaps;
    snaps.push_back({sf.Pi, sf.Phi});
    minA = 1.;
    std::vector<double> A, delta;
    for (int k = 1; k <= nSamples; k++) {
        for (int s = 0; s < stepsPerSample; s++) sf = evolveRK4(sf, &einsteinRHS);
        solveGeometry(sf.Pi, sf.Phi, sf.dr, A, delta);
        for (double a : A) minA = std::min(minA, a);
        snaps.push_back({sf.Pi, sf.Phi});
    }
    return snaps;
}

static double l2(const std::vector<double> &coarse, const std::vector<double> &fine, double dr)
{
    double sum = 0.;
    for (size_t i = 0; i < coarse.size(); i++) {
        double d = coarse[i] - fine[2 * i];
        sum += d * d * dr;
    }
    return std::sqrt(sum);
}

int main(int argc, char **argv)
{
    if (argc > 1) amp = std::atof(argv[1]);
    if (argc > 2) sigma = std::atof(argv[2]);
    std::cout << std::setprecision(5);
    std::cout << "amplitude = " << amp << ", sigma = " << sigma << "\n";

    constraintCheck();

    const std::vector<int> ns = {250, 500, 1000, 2000};
    const int nSamples = 20; // t_k = k*pi/20, dt = pi/(4n) -> n/5 steps per sample
    std::vector<std::vector<std::vector<std::vector<double>>>> data;
    std::vector<double> minA(ns.size());
    for (size_t j = 0; j < ns.size(); j++) {
        std::cerr << "running n = " << ns[j] << std::endl;
        data.push_back(run(ns[j], nSamples, ns[j] / 5, minA[j]));
    }

    mkdir("convergence", 0755);
    std::ofstream ord("convergence/einstein_order.dat");
    ord << std::setprecision(10);
    ord << "# amplitude " << amp << " sigma " << sigma << "\n# t  order(250/500/1000)  order(500/1000/2000)\n";
    std::cout << "\n== self-convergence of Pi (L2), min A over the run (n=2000): " << minA[3] << "\n";
    std::cout << "t\torder(250/500/1000)\torder(500/1000/2000)\n";
    for (int k = 1; k <= nSamples; k++) {
        double e[3];
        for (int j = 0; j < 3; j++)
            e[j] = l2(data[j][k][0], data[j + 1][k][0], M_PI_2 / ns[j]);
        double o1 = std::log2(e[0] / e[1]), o2 = std::log2(e[1] / e[2]);
        ord << k * M_PI / 20. << "\t" << o1 << "\t" << o2 << "\n";
        std::cout << k * M_PI / 20. << "\t" << o1 << "\t\t" << o2 << "\n";
    }
    return 0;
}
