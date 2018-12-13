#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <random>
#include <cassert>
#include <array>
#include <vector>

void question1() {
    const double x = 1.0;
    const double y = 1.0;

    const double r = std::sqrt(std::pow(x,2) + std::pow(y,2));

    std::cout << "r = "<< r << std::endl;
}

void question2() {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniform(-1.0,1.0);
    const int N = 1e6;
    int count = 0;
    for (int i = 0; i < N; ++i) {
        const double x = uniform(generator);
        const double y = uniform(generator);
        const double r2 = std::pow(x,2) + std::pow(y,2);

        if (r2 < 1.0) {
            ++count;
        }
    }

    std::cout << "pi is about "<< 4.0*static_cast<double>(count)/static_cast<double>(N) << std::endl;
}

void question3() {
    const int N = 1000;
    double sum = 0.0;
    for (int i = 1; i < N; ++i) {
        sum += 1.0/static_cast<double>(std::pow(i,2));
    }

    std::cout << "pi is about "<< std::sqrt(6.0*sum) << std::endl;
}

void question4() {
    const int N = 100;
    double a = 1.0;
    double b = 1.0/std::sqrt(2);
    double t = 0.25;
    double p = 1.0;
    for (int i = 1; i < N; ++i) {
        const double an = a;
        const double bn = b;
        a = (an + bn)/2;
        b = std::sqrt(an*bn);
        t -= p*std::pow(a - an,2);
        p *= 2;
    }

    std::cout << "pi is about "<< std::pow(a + b,2)/(4*t) << std::endl;
}

void question5() {
    const double h = 1.0/100.0;
    double y = 1.0;
    std::ofstream output("xy.dat");
    assert(output.is_open());
    for (double x = 0; x < 1; x += h) {
        y /= 1.0 + h;
        std::cout << "outputting x = "<<x<<" y = "<< y << std::endl;
        output << x << " " << y << std::endl;
    }
    output.close();
}

void question6() {
    std::string line;
    std::ifstream input("xy.dat");
    assert(input.is_open());
    double max_error = std::numeric_limits<double>::max();
    double x,y;
    while (std::getline (input,line)) {
        std::istringstream s(line);
        s >> x;
        s >> y;
        std::cout << "reading x = "<<x<<" y = "<< y << std::endl;
        const double error = std::abs(y - std::exp(-x));
        if (max_error > error) {
            max_error = error;
        }
    }
    std::cout << "max error = "<< max_error << std::endl;
    input.close();
}


void question7() {
    std::array<double,3> x = {1.0, 2.0, 3.0};
    std::array<double,3> y = {1.0, 2.0, 3.0};

    double dot = 0.0;
    for (int i = 0; i < 3; ++i) {
        dot += x[i]*y[i];
    }

    std::cout << "dot with arrays = "<< dot << std::endl;
}

void question7b(const int n) {
    std::vector<double> x(n);
    std::vector<double> y(n);

    for (int i = 0; i < n; ++i) {
        x[i] = i;
        y[i] = i;
    }

    double dot = 0.0;
    for (int i = 0; i < n; ++i) {
        dot += x[i]*y[i];
    }

    std::cout << "dot with vectors = "<< dot << std::endl;
}

void question8() {
    std::array<std::array<double,3>,3> A = {{{5, 8, 2}, {8, 3, 1}, {5, 3, 9}}};
    std::array<std::array<double,3>,3> B = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    std::array<std::array<double,3>,3> C = {};

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                C[i][j] += A[i][k] * B[k][j]; 
            }
        }
    }

    std::cout << "C = " << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << "| ";
        for (int j = 0; j < 3; ++j) {
            std::cout << C[i][j];
            if (j == 2) {
                std::cout << " |" << std::endl;
            } else {
                std::cout << ", ";
            }
        }
    }
}


int main(void) {
    question1();
    question2();
    question3();
    question4();
    question5();
    question6();
    question7();
    question7b(100);
    question8();
        
    return 0;
}
