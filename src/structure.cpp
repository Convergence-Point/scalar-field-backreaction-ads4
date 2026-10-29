#include "structure.hpp"

double gaussTan(double x, double sigma){
    return pow(M_E, -tan(x) * tan(x) / sigma / sigma / M_PI / M_PI);
}
