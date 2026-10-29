#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "derivative.hpp"
#include "structure.hpp"
#include "evolution.hpp"
#include "einstein.hpp"
#include "progress.hpp"
#include <cstring>
#include <cstdlib>
#include <map>

using namespace std;


ScalarField initialize(int n, double sigma, double amplitude, double (*initialPi)(double, double)) {
    double rMax = M_PI_2;

    vector<double> r;
    ScalarField sf;
    sf.t = 0.;
    sf.dr = rMax / double(n);
    sf.dt = sf.dr / 2.;

    sf.r.resize(n+1);
    sf.Phi.resize(n+1);
    sf.Pi.resize(n+1);

    for (int i = 0; i < sf.Phi.size(); i++){
        sf.r[i] = sf.dr * double(i);
        sf.Phi[i] = 0.;

        double x = sf.r[i];
        sf.Pi[i] = amplitude * initialPi(x, sigma);
    }
    
    return sf;
}

void dPi2File(string filename, ScalarField &isf)
{
    int n = isf.Pi.size();
    ofstream file(filename);
    file << std::setprecision(15);

    vector<double> dPi;
    dPi.resize(n);
    vector<double> dPhi;
    dPhi.resize(n);
    for (int i = 0; i < n; i++)
    {
        // AdS parities: Pi even/odd and Phi odd/even at x = 0 / x = pi/2
        dPi.at(i) = derP(isf.Pi, i, isf.dr, 1., -1.);
        dPhi.at(i) = derP(isf.Phi, i, isf.dr, -1., 1.);
        double x = isf.r[i];
        file << x << "\t" << isf.Pi[i] << "\t" << dPi[i] << "\t" << isf.Phi[i] << "\t" << dPhi[i] << endl;
    }

    file.close();
}

// usage: ./build/tests/test_convergence_point [ads|einstein] [amplitude]      defaults: ads 1
//   ads      : fixed AdS background, A = 1, delta = 0 (linear wave equation)
//   einstein : A and delta from the constraints at every step (backreaction).
//              The run is stopped (exit code 3) as soon as min A < 1/n or A is not finite,
//              i.e. the fields are so intense that a horizon is forming.
// See einstein.cpp for the equations. Both use the same RHS code (computeRHS).
int main(int argc, char **argv)
{
    std::string mode = argc > 1 ? argv[1] : "ads";
    double amplitude = argc > 2 ? std::atof(argv[2]) : 1.;
    RHSFunction rhs;
    if (mode == "ads") rhs = &adsRHS;
    else if (mode == "einstein") rhs = &einsteinRHS;
    else { cerr << "mode must be 'ads' or 'einstein'" << endl; return 1; }
    cout << "mode = " << mode << ", amplitude = " << amplitude << endl;

    std::map<double, double> my_map = {};
    double minA = 1.;

    for(int n = 250; n <= 1000; n*=2){
        ScalarField isf = initialize(n, 1./16., amplitude, &gaussTan);

        dPi2File("is" + to_string(n) + ".txt", isf);

        // evolve every resolution to the same final time T = 4n * dt = pi
        // (dt = dr/2 = pi/(4n)), otherwise the runs are compared at different times
        ScalarField iter1 = isf;
        if (mode == "einstein"){
            double a0 = minimumA(isf.Pi, isf.Phi, isf.dr);
            if (!(a0 >= 1. / n)){
                cout << "COLLAPSE: the initial data already has min A = " << a0 << " < 1/n = " << 1. / n
                     << " (n = " << n << "); lower the amplitude" << endl;
                return 3;
            }
        }
        for (int step = 0; step < 4 * n; step++){
            iter1 = evolveRK4(iter1, rhs);
            {
                char text[96];
                std::snprintf(text, sizeof(text), "n=%d  t = %.3f / %.3f", n, iter1.t, 4. * n * iter1.dt);
                progressUpdate(double(step + 1) / double(4 * n), text);
            }
            if (mode == "einstein"){
                int where = 0;
                double a = minimumA(iter1.Pi, iter1.Phi, iter1.dr, &where);
                if (!(a >= 1. / n)){
                    progressEnd();
                    cout << "COLLAPSE: " << (std::isnan(a) ? "A is not finite (numerical blow-up: collapse or dt too large for these fields)"
                                                           : "min A = " + to_string(a) + " < 1/n = " + to_string(1. / n))
                         << " at x = " << iter1.r[where] << " (i = " << where << "), t = " << iter1.t
                         << ", step " << step + 1 << " of " << 4 * n << ", n = " << n << endl;
                    dPi2File("iter" + to_string(n) + ".txt", iter1);
                    cout << "state written to iter" << n << ".txt; no convergence ratio computed" << endl;
                    return 3;
                }
                minA = std::min(minA, a);
            }
        }

        progressEnd();
        dPi2File("iter" + to_string(n) + ".txt", iter1);
        my_map[n] = iter1.Pi[n/10];
    }

    std::cout << std::fixed;
    std::cout << std::setprecision(15);
    cout << my_map[250] << "\t" << my_map[500] << "\t" << my_map[1000] << endl;
    cout << "ratio = " << log(abs((my_map[250] - my_map[500]) / (my_map[500] - my_map[1000]))) / log(2)
         << endl;

    if (mode == "einstein")
        cout << "min A over all runs (checked every step) = " << minA
             << "   (close to 0 means collapse)" << endl;

    return 0;
}
