#include <cmath>
#include <iostream>
#include <cassert>
#include <array>
#include <vector>
#include <Eigen/Core>
#include <chrono>
#include <random>

void question1_f(int *in) {
    std::cout << "old value = "<< *in << std::endl;
    *in = 4;
}

void question1() {
    int x = 3;
    question1_f(&x);
    std::cout << "new value = "<< x << std::endl;
}

void question2_fp(float *x, float *y) {
    float tmp = *x;
    *x = *y;
    *y = tmp;
}

void question2_fr(float &x, float &y) {
    float tmp = x;
    x = y;
    y = tmp;
}

void question2() {
    float x = 1.0;
    float y = 2.0;
    std::cout << "before swap: x = "<< x <<" y = "<< y << std::endl;
    question2_fp(&x,&y);
    std::cout << "after swap: x = "<< x <<" y = "<< y << std::endl;
    question2_fr(x,y);
    std::cout << "after swap: x = "<< x <<" y = "<< y << std::endl;
}

double dot_product(const std::array<double,3>& x, const std::array<double,3>& y) {
    double dot = 0.0;
    for (unsigned int i = 0; i < x.size(); ++i) {
        dot += x[i]*y[i];
    }
    return dot;
}

void question3() {
    const int n = 100;
    std::array<double,3> x = {1, 1, 1};
    std::array<double,3> y = {1, 1, 1};
    
    std::cout << "dot product = "<< dot_product(x,y) << std::endl; 
}

template <long unsigned int N>
double p_norm(const std::array<double,N>& x, const unsigned int p = 2) {
    double result = 0.0;
    for (unsigned i = 0; i < x.size(); ++i) {
        result += std::pow(std::abs(x[i]),p); 
    }

    result = 0.0;
    for (double val: x) {
        result += std::pow(std::abs(val),p); 
    }

    result = 0.0;
    for (auto it = std::begin(x); it != std::end(x); ++it) {
        result += std::pow(std::abs(*it),p); 
    }

    result = std::accumulate(std::begin(x),std::end(x),0.0,
                [&](double accum, double val) { 
                    return accum + std::pow(std::abs(val),p); 
                });

    return result;
}

void question4() {
    std::array<double,5> x = {1, 2, 3, 4, 5};
    std::cout << "p_norm = "<< p_norm(x) << std::endl;
}

void question5() {
    const size_t N = 3;
    auto lp_norm = [](const std::array<double,N>& x, const unsigned int p = 2) {
        double result = 0.0;
        for (double val: x) {
            result += std::pow(std::abs(val),p); 
        }
        return result;
    };
    const double p = 2;
    auto lp_norm_cap = [p](const std::array<double,N>& x) {
        double result = 0.0;
        for (double val: x) {
            result += std::pow(std::abs(val),p); 
        }
        return result;
    };
    std::array<double,N> x = {1, 2, 3};

    std::cout << "p_norm = "<< lp_norm(x) << std::endl;
    std::cout << "p_norm = "<< lp_norm_cap(x) << std::endl;
}

double p_norm(const std::vector<double>& x, const unsigned int p = 2) {
    double result = 0.0;
    for (unsigned i = 0; i < x.size(); ++i) {
        result += std::pow(std::abs(x[i]),p); 
    }

    result = 0.0;
    for (double val: x) {
        result += std::pow(std::abs(val),p); 
    }

    result = 0.0;
    for (auto it = std::begin(x); it != std::end(x); ++it) {
        result += std::pow(std::abs(*it),p); 
    }

    result = std::accumulate(std::begin(x),std::end(x),0.0,
                [&](double accum, double val) { 
                    return accum + std::pow(std::abs(val),p); 
                });

    return result;
}


void question6() {
    const size_t n = 100;
    std::vector<double> x(n,1);
    std::cout << "p_norm = "<< p_norm(x) << std::endl;
}


std::vector<double> multiply(std::vector<double>& A, 
                             const std::array<size_t,2> Asize,
                             std::vector<double>& B,
                             const std::array<size_t,2> Bsize) {
    assert(Asize[1] == Bsize[0]);
    assert(Asize[0]*Asize[1] == A.size());
    assert(Bsize[0]*Bsize[1] == A.size());
    std::vector<double> C(Asize[0]*Bsize[1],0);

    for (size_t i = 0; i < Asize[0]; ++i) {
        for (size_t j = 0; j < Bsize[1]; ++j) {
            for (size_t k = 0; k < Asize[1]; ++k) {
                C[i*Bsize[1]+j] += A[i*Asize[1]+k] * B[k*Bsize[1]+j]; 
            }
        }
    }
    return C;
}


void question7() {
    std::vector<double> A = {5, 8, 2, 8, 3, 1, 5, 3, 9};
    std::vector<double> B = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    std::vector<double> C = multiply(A, {3,3}, B, {3,3});

    std::cout << "C = " << std::endl;
    for (int i = 0; i < 3; ++i) {
        std::cout << "| ";
        for (int j = 0; j < 3; ++j) {
            std::cout << C[i*3 + j];
            if (j == 2) {
                std::cout << " |" << std::endl;
            } else {
                std::cout << ", ";
            }
        }
    }
}

// implements C += A*B
void multiply_and_add(const std::vector<double>& A, 
                      const std::array<size_t,2> Asize,
                      const std::vector<double>& B,
                      const std::array<size_t,2> Bsize,
                      std::vector<double>& C) {
    assert(Asize[1] == Bsize[0]);
    assert(Asize[0]*Asize[1] == A.size());
    assert(Bsize[0]*Bsize[1] == A.size());
    assert(C.size() == Asize[0]*Bsize[1]);

    for (size_t i = 0; i < Asize[0]; ++i) {
        for (size_t j = 0; j < Bsize[1]; ++j) {
            for (size_t k = 0; k < Asize[1]; ++k) {
                C[i*Bsize[1]+j] += A[i*Asize[1]+k] * B[k*Bsize[1]+j]; 
            }
        }
    }
}

// block = A[block_row_index:block_row_index+block_size[0],
//           block_col_index:block_col_index+block_size[1]]
void copy_block(const std::vector<double>& A, 
            const std::array<size_t,2> Asize,
            const std::array<size_t,2> block_size,
            const size_t block_row_index, 
            const size_t block_col_index,
            std::vector<double>& block) {
 
    const size_t real_row_index = block_row_index*block_size[0];
    const size_t real_col_index = block_col_index*block_size[1];
    for (size_t i = 0; i < block_size[0]; ++i) {
        for (size_t j = 0; j < block_size[1]; ++j) {
            block[i*block_size[1] + j] = A[real_row_index*Asize[1]+real_col_index];
        }
    }
}

// block = 0
void zero_block(std::vector<double>& block, 
            const std::array<size_t,2> block_size) {
 
    for (size_t i = 0; i < block_size[0]; ++i) {
        for (size_t j = 0; j < block_size[1]; ++j) {
            block[i*block_size[1] + j] = 0;
        }
    }
}

// A[block_row_index:block_row_index+block_size[0],
//   block_col_index:block_col_index+block_size[1]] = block
void write_block(std::vector<double>& A, 
                 const std::array<size_t,2> Asize,
                 const std::vector<double>& block, 
                 const std::array<size_t,2> block_size,
                 const size_t block_row_index, 
                 const size_t block_col_index) {
 
    const size_t real_row_index = block_row_index*block_size[0];
    const size_t real_col_index = block_col_index*block_size[1];
    for (size_t i = 0; i < block_size[0]; ++i) {
        for (size_t j = 0; j < block_size[1]; ++j) {
            A[real_row_index*Asize[1]+real_col_index] = block[i*block_size[1] + j];
        }
    }
}

// returns A*B using blocking, or tiling technique
std::vector<double> multiply_blocking(const std::vector<double>& A, 
                             const std::array<size_t,2> Asize,
                             const std::vector<double>& B,
                             const std::array<size_t,2> Bsize,
                             const std::array<size_t,2> block_size) {
    assert(Asize[1] == Bsize[0]);
    assert(Asize[0]*Asize[1] == A.size());
    assert(Bsize[0]*Bsize[1] == A.size());
    std::array<size_t,2> Csize = {Asize[0],Bsize[1]};
    std::vector<double> C(Asize[0]*Bsize[1],0);

    const size_t total_block_size = block_size[0]*block_size[1];

    // pre-allocate blocks
    std::vector<double> Ablock(total_block_size);
    std::vector<double> Bblock(total_block_size);
    std::vector<double> Cblock(total_block_size);

    // loop through all blocks (i,j) calculating:
    // C(i,j) = A(i,k)*B(k,j)
    for (size_t i = 0; i < Asize[0]/block_size[0]; ++i) {
        for (size_t j = 0; j < Bsize[1]/block_size[1]; ++j) {
            zero_block(Cblock,block_size);
            for (size_t k = 0; k < Asize[1]/block_size[1]; ++k) {
                copy_block(A,Asize,block_size,i,k,Ablock);
                copy_block(B,Bsize,block_size,k,j,Bblock);
                multiply_and_add(Ablock,block_size,
                                 Bblock,block_size,
                                 Cblock);
            }
            write_block(C,Csize,Cblock,block_size,i,j);
        }
    }
    return C;
}

void question8() {
    const size_t n = 1000;
    std::vector<double> A(n*n);
    std::vector<double> B(n*n);
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniform(0.0,1.0);
    for (size_t i = 0; i < n*n; ++i) {
        A[i] = uniform(generator);
        B[i] = uniform(generator);
    }

    // Calculate C = A*B using naive method
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> C = multiply(A,{n,n},B,{n,n});
    auto t2 = std::chrono::high_resolution_clock::now();
    auto my_time = std::chrono::duration_cast<
                        std::chrono::microseconds>(t2 - t1).count();

    // Calculate C = A*B using blocking
    t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> C2 = multiply_blocking(A,{n,n},B,{n,n},{n/10,n/10});
    t2 = std::chrono::high_resolution_clock::now();
    auto my_blocking_time = std::chrono::duration_cast<
                        std::chrono::microseconds>(t2 - t1).count();

    Eigen::MatrixXd Aeigen = Eigen::MatrixXd::Random(n,n);
    Eigen::MatrixXd Beigen = Eigen::MatrixXd::Random(n,n);

    // Calculate C = A*B using Eigen
    t1 = std::chrono::high_resolution_clock::now();
    Eigen::MatrixXd Ceigen = Aeigen*Beigen;
    t2 = std::chrono::high_resolution_clock::now();
    auto eigen_time = std::chrono::duration_cast<
                        std::chrono::microseconds>(t2 - t1).count();

    // Output timings
    std::cout << "multiply without blocking = "<< my_time << 
                 " versus multiply with blocking = "<< my_blocking_time <<
                 " versus eigen time = "<< eigen_time << std::endl;
}

int main(void) {
    question1();
    question2();
    question3();
    question4();
    question5();
    question6();
    question7();
    question8();
        
    return 0;
}
