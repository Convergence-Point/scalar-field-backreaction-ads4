#ifndef DERIVATIVE_HPP
#define DERIVATIVE_HPP

#include <vector>
#include <array>
#include <iostream>
#include <cstdlib>
#include <numeric>
#include <cmath>

using namespace std;

enum Field {PiField, PhiField};
double der(const vector<double> &v, int p, double dx, Field field);

double derP(const vector<double> &v, int p, double dx, double leftParity, double rightParity);

// AdS4 operators (x in [0, pi/2], phi = f rho^3 + O(rho^5) at pi/2); see derivative.cpp
void adsDPi(const vector<double> &Pi, double dr, vector<double> &out);
void adsDivergence(const vector<double> &Phi, double dr, vector<double> &out);
// The same system in "expanded" form, symmetrized with Theta = tan x Pi, Psi = tan x Phi so that it is
// stable (see derivative.cpp). Use the two together: Pi_t from adsDivergenceExpanded, Phi_t from adsDPiExpanded.
void adsDivergenceExpanded(const vector<double> &Phi, double dr, vector<double> &out);
void adsDPiExpanded(const vector<double> &Pi, double dr, vector<double> &out);

#endif
