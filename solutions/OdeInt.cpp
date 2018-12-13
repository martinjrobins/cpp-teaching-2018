#include "OdeInt.h"

double integrate(const double x0, const double t1, const double dt) {
    OdeInt stepper;
    stepper.setdt(dt);
    stepper.setx(x0);
    double t;
    for (t = 0; t < t1; t += dt) {
        stepper.step();
    }
    stepper.setdt(t1-t);
    stepper.step();

    return stepper.getx();
}

void OdeInt::dxdt(const double x, double& dxdt) {
    dxdt = -x;
}

void OdeInt::step() {
    double dxdt_at_x;
    dxdt(x,dxdt_at_x);
        
    x += dt*dxdt_at_x;
}



