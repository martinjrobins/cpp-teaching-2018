#ifndef ODEINT_H
#define ODEINT_H

#include <vector>
#include <iostream>

template<typename F, typename T>
class OdeInt {
    public:
        OdeInt(const F dxdt):
            dxdt(dxdt) {};
        void step();
        const T& getx() { return x; };
        void setx(const T& new_x) { x = new_x; };
        double getdt() { return dt; };
        void setdt(const double new_dt) { dt = new_dt; };
    private:
        T x;
        double dt;
        F dxdt;
};

template<typename F, typename T>
const T integrate(const T& x0, const double t1, const double dt, const F& dxdt) {
    OdeInt<F,T> stepper(dxdt);
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

template<typename F, typename T>
void OdeInt<F,T>::step() {
    T dxdt_at_x;
    dxdt(x,dxdt_at_x);
        
    x += dt*dxdt_at_x;
}

/*
 * std::vector specialisation
 */
template<typename F, typename T>
class OdeInt<F,std::vector<T> > {
    public:
        OdeInt(const F dxdt):
            dxdt(dxdt) {};
        void step();
        const std::vector<T>& getx() { return x; };
        void setx(const std::vector<T>& new_x) { 
            x.resize(new_x.size());
            std::copy(new_x.begin(), new_x.end(), x.begin()); 
        };
        double getdt() { return dt; };
        void setdt(const double new_dt) { dt = new_dt; };
    private:
        std::vector<T> x;
        double dt;
        F dxdt;
};


template<typename F, typename T>
void OdeInt<F,std::vector<T> >::step() {
    std::vector<T> dxdt_at_x(x.size());
    dxdt(x,dxdt_at_x);
        
    for (int i = 0; i < x.size(); ++i) {
        x[i] += dt*dxdt_at_x[i];
    }
}

#endif
