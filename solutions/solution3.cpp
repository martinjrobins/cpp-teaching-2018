
#include "OdeInt.h"
#include <cassert>
#include <cmath>
#include <iostream>

int main(void) {
    for (double t = 0; t < 10; t += 1) {
        const double x0 = 1;
        const double dt = 1e-3;
        const double result = integrate(x0,t,dt);
        const double truth = std::exp(-t);
        const double error = std::abs((result - truth)/truth);
        assert(error < 1e-2);
    }

    std::cout << "All done, success!" << std::endl;
    
    return 0;
}
