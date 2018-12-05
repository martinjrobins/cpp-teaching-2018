#include "mex.h"
#include <iostream>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[])
{
	mxArray *v = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *data = mxGetPr(v);
	*data = 3.142;
	std::cout<<"Num args = "<<nrhs<<" \n";
	plhs[0] = v;
}

