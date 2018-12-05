#ifndef VECTORDEF
#define VECTORDEF

#include <vector>



//  **********************
//  *  Class of vectors  *
//  **********************


//  Class written in such a way that code similar to Matlab
//  code may be written


#include <cmath>
#include "Exception.hpp"//  This class throws errors using the class "error"

class Vector;

// function prototypes
double norm(Vector& v, int p=2);
// Prototype signature of length() friend function
size_t length(const Vector& v);


class Vector
{
    private:
        // member variables
        std::vector<double> mData;   // data stored in vector

    public:
        // constructors
        // No default constructor
        // overridden copy constructor
        Vector(const Vector& v1);
        // construct vector of given size
        Vector(size_t sizeVal);

        // binary operators
        friend Vector operator+(const Vector& v1, 
                const Vector& v2);
        friend Vector operator-(const Vector& v1, 
                const Vector& v2);
        friend double operator*(const Vector& v1, 
                const Vector& v2);
        friend Vector operator*(const Vector& v, const double& a);
        friend Vector operator*(const double& a, const Vector& v);
        friend Vector operator/(const Vector& v, const double& a);

        // unary operators
        friend Vector operator+(const Vector& v);
        friend Vector operator-(const Vector& v);

        //other operators
        //assignment
        Vector& operator=(const Vector& v);
        //indexing
        double& operator()(size_t i);
        //output
        friend std::ostream& operator<<(std::ostream& output, const Vector& v);

        //norm (as a member method)
        double norm(int p=2) const;
        // functions that are friends
        friend double norm(Vector& v, int p);
        friend size_t length(const Vector& v);
};


#endif
