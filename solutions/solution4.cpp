
#include "OdeIntTemplate.h"
#include <cassert>
#include <cmath>
#include <vector>
#include <iostream>

#include <Eigen/Dense>

struct negx {
    void operator()(const double &x, double &dxdt) { dxdt = -x; }
};

typedef std::vector<double> state_type;
typedef Eigen::Matrix<double, 2, 1> eigen_type;

template<typename T>
struct spring {
    void operator()(const T& x, T& dxdt) { 
        dxdt[0] = x[1]; 
        dxdt[1] = -x[0]; 
    }
};

int main(void) {

    for (double t = 0; t < 10; t += 1) {
        const double x0 = 1;
        const double dt = 1e-3;
        const double result = integrate(x0,t,dt,negx());
        const double truth = std::exp(-t);
        const double error = std::abs((result - truth)/truth);
        assert(error < 1e-2);
    }

    for (double t = 0; t < 10; t += 1) {
        state_type x0(2);
        x0[0] = 1;
        x0[1] = 0;
        const double dt = 1e-3;
        state_type result = integrate(x0,t,dt,spring<state_type>());
        const double truth = std::cos(t);
        const double error = std::abs((result[0] - truth)/truth);
        assert(error < 1e-2);
    }

    for (double t = 0; t < 10; t += 1) {
        eigen_type x0;
        x0[0] = 1;
        x0[1] = 0;
        const double dt = 1e-3;
        eigen_type result = integrate(x0,t,dt,spring<eigen_type>());
        const double truth = std::cos(t);
        const double error = std::abs((result[0] - truth)/truth);
        assert(error < 1e-2);
    }

    std::cout << "All done, success!" << std::endl;
    
    return 0;
}
