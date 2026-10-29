#ifndef STRUCTURE_HPP
#define STRUCTURE_HPP

#include <vector>
#include <cmath>

struct ScalarField
{
    double t;
    double dt;
    std::vector<double> r;
    double dr;
    std::vector<double> Phi;
    std::vector<double> Pi;
};

double gaussTan(double x, double sigma);

#endif