#include <new>          // operator new[]
#include <stdlib.h>     // functions malloc and free
#include <vector>	// std::vector class


void create_array(int n) {
	/*
	 * compiler error here
	 */
	double x0[n];

	/*
	 * inefficient way
         */
	double x1[1000];
	assert(n<1000);

	/*
     * the c way
     */
	double *x2 = (double *) malloc (n*sizeof(double));
	free(x2);

	/*
	 * the c++ way
	 */
	double *x3 = new double[n];
	delete[] x3;

	/*
	 * the c++ STL way
	 */
	std::vector<double> x4(n); //automatic memory management
}

	

	
	
	
