#ifndef ODEINT_H
#define ODEINT_H

#include <vector>
#include <iostream>

class OdeInt {
    public:
        void step();
        double getx() { return x; };
        void setx(const double new_x) { x = new_x; };
        double getdt() { return dt; };
        void setdt(const double new_dt) { dt = new_dt; };
    private:
        void dxdt(const double x, double& dxdt);
        double x;
        double dt;
};

double integrate(const double x0, const double t1, const double dt);

#endif
