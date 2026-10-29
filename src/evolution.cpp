#include "evolution.hpp"
#include "derivative.hpp"

#include <iostream>


ScalarField evolve(
    ScalarField sf,
    double (*dPi)(double, double, double, double, double, double),
    double (*dPhi)(double, double, double, double, double, double),
    PhiOperator phiOp,
    PiOperator piOp
){
    // value passed to the callbacks as "dPhidx": plain d/dx by default,
    // or a user supplied operator acting on the whole Phi vector (e.g. a divergence)
    std::vector<double> dPhik;
    auto computeDPhi = [&](std::vector<double> &v){
        if (phiOp){
            phiOp(v, sf.dr, dPhik);
        } else {
            dPhik.resize(v.size());
            for (size_t j = 0; j < v.size(); j++) dPhik[j] = der(v, j, sf.dr, PhiField);
        }
    };

    std::vector<double> dPik;
    auto computeDPi = [&](std::vector<double> &v){
        if (piOp){
            piOp(v, sf.dr, dPik);
        } else {
            dPik.resize(v.size());
            for (size_t j = 0; j < v.size(); j++) dPik[j] = der(v, j, sf.dr, PiField);
        }
    };

    std::vector<double> k1Pi, k2Pi, k3Pi, k4Pi;
    std::vector<double> k1Phi, k2Phi, k3Phi, k4Phi;
    std::vector<double> Pik, Phik;
    int s = sf.r.size();
    k1Pi.resize(s);
    k2Pi.resize(s);
    k3Pi.resize(s);
    k4Pi.resize(s);
    k1Phi.resize(s);
    k2Phi.resize(s);
    k3Phi.resize(s);
    k4Phi.resize(s);
    Pik.resize(s);
    Phik.resize(s);

    // k1
    for (int i = 0; i < s; i++){
        Pik[i] = sf.Pi[i];
        Phik[i] = sf.Phi[i];
    }
    computeDPi(Pik);
    computeDPhi(Phik);

    #pragma omp parallel for
    for (int i = 0; i < s; i++){
        k1Pi[i] = dPi(
            sf.t,
            sf.r[i],
            Pik[i],
            Phik[i],
            dPik[i],
            dPhik[i]
        );
        k1Phi[i] = dPhi(
            sf.t,
            sf.r[i],
            Pik[i],
            Phik[i],
            dPik[i],
            dPhik[i]
        );
    }

    // k2
    #pragma omp parallel for
    for (int i = 0; i < s; i++){
        Pik[i] = sf.Pi[i] + k1Pi[i] * sf.dt / 2.;
        Phik[i] = sf.Phi[i] + k1Phi[i] * sf.dt / 2.;
    }
    computeDPi(Pik);
    computeDPhi(Phik);
    #pragma omp parallel for
    for (int i = 0; i < s; i++){
        k2Pi[i] = dPi(
            sf.t + sf.dt/2.,
            sf.r[i],
            Pik[i],
            Phik[i],
            dPik[i],
            dPhik[i]
        );
        k2Phi[i] = dPhi(
            sf.t + sf.dt/2.,
            sf.r[i],
            Pik[i],
            Phik[i],
            dPik[i],
            dPhik[i]
        );
    }

    // k3
    #pragma omp parallel for
    for (int i = 0; i < s; i++){
        Pik[i] = sf.Pi[i] + k2Pi[i] * sf.dt / 2.;
        Phik[i] = sf.Phi[i] + k2Phi[i] * sf.dt / 2.;
    }
    computeDPi(Pik);
    computeDPhi(Phik);
    #pragma omp parallel for
    for (int i = 0; i < s; i++){
        k3Pi[i] = dPi(
            sf.t + sf.dt/2.,
            sf.r[i],
            Pik[i],
            Phik[i],
            dPik[i],
            dPhik[i]
        );
        k3Phi[i] = dPhi(
            sf.t + sf.dt/2.,
            sf.r[i],
            Pik[i],
            Phik[i],
            dPik[i],
            dPhik[i]
        );
    }

    // k4
    #pragma omp parallel for
    for (int i = 0; i < s; i++){
        Pik[i] = sf.Pi[i] + k3Pi[i] * sf.dt;
        Phik[i] = sf.Phi[i] + k3Phi[i] * sf.dt;
    }
    computeDPi(Pik);
    computeDPhi(Phik);
    #pragma omp parallel for
    for (int i = 0; i < s; i++){
        k4Pi[i] = dPi(
            sf.t + sf.dt,
            sf.r[i],
            Pik[i],
            Phik[i],
            dPik[i],
            dPhik[i]
        );
        k4Phi[i] = dPhi(
            sf.t + sf.dt,
            sf.r[i],
            Pik[i],
            Phik[i],
            dPik[i],
            dPhik[i]
        );
    }

    // evolution
    ScalarField nextIter = sf;
    nextIter.t += nextIter.dt;

    std::vector<double> factors {1./6., 2./6., 2./6., 1./6.};
    for (int i = 0; i < s; i++){
        std::vector<double> kPhi {k1Phi[i], k2Phi[i], k3Phi[i], k4Phi[i]};

        nextIter.Phi[i] += sf.dt * std::inner_product(
            kPhi.begin(), kPhi.end(), factors.begin(), 0.
        );

        std::vector<double> kPi {k1Pi[i], k2Pi[i], k3Pi[i], k4Pi[i]};
        nextIter.Pi[i] += sf.dt * std::inner_product(
            kPi.begin(), kPi.end(), factors.begin(), 0.
        );

        // if(i == 0){
        //     std::cout << "kPi" << endl;
        //     for (double i : kPi)
        //         std::cout << i << '\t';
        //     std::cout << endl;
        // }
    }

    return nextIter;
}

ScalarField evolveRK4(ScalarField sf, RHSFunction rhs){
    const size_t s = sf.r.size();
    const double dt = sf.dt;
    std::vector<double> k1Pi, k2Pi, k3Pi, k4Pi, k1Phi, k2Phi, k3Phi, k4Phi;
    std::vector<double> Pik(s), Phik(s);

    rhs(sf.Pi, sf.Phi, sf.dr, k1Pi, k1Phi);
    #pragma omp parallel for
    for (size_t i = 0; i < s; i++){
        Pik[i] = sf.Pi[i] + 0.5 * dt * k1Pi[i];
        Phik[i] = sf.Phi[i] + 0.5 * dt * k1Phi[i];
    }
    rhs(Pik, Phik, sf.dr, k2Pi, k2Phi);
    #pragma omp parallel for
    for (size_t i = 0; i < s; i++){
        Pik[i] = sf.Pi[i] + 0.5 * dt * k2Pi[i];
        Phik[i] = sf.Phi[i] + 0.5 * dt * k2Phi[i];
    }
    rhs(Pik, Phik, sf.dr, k3Pi, k3Phi);
    #pragma omp parallel for
    for (size_t i = 0; i < s; i++){
        Pik[i] = sf.Pi[i] + dt * k3Pi[i];
        Phik[i] = sf.Phi[i] + dt * k3Phi[i];
    }
    rhs(Pik, Phik, sf.dr, k4Pi, k4Phi);

    ScalarField next = sf;
    next.t += dt;
    #pragma omp parallel for
    for (size_t i = 0; i < s; i++){
        next.Pi[i] += dt / 6. * (k1Pi[i] + 2. * k2Pi[i] + 2. * k3Pi[i] + k4Pi[i]);
        next.Phi[i] += dt / 6. * (k1Phi[i] + 2. * k2Phi[i] + 2. * k3Phi[i] + k4Phi[i]);
    }
    return next;
}
