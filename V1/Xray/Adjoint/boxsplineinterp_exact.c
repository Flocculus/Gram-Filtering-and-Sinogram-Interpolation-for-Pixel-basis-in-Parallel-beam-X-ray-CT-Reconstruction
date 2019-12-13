#include "mex.h"
#include "math.h"
#include <stdio.h>
const double PI  =3.141592653589793238463;
///////////246

double mypow
(
	double x,
	int n
	)
{
	double tmp = x * x;
	if (n == 5)
	{
		return tmp * tmp * x; 
	}
	if (n == 3)
	{
		return tmp * x;
	}
}

double interp20
(
	double point,
	double* shift,
	double* coe
	)
{
	double V = 0;
	if (point <= shift[0])
	{
		V = 0;
	}
	else if (point <= shift[1])
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0;
	}
	else if (point <= shift[2])
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0 + coe[1] * mypow((point-shift[1]),5) / 120.0;
	}
	else if (point <= shift[3])
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0 + coe[1] * mypow((point-shift[1]),5) / 120.0 + coe[2] * mypow((point-shift[2]),5) / 120.0;
	}
	else if (point <= shift[4])
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0 + coe[1] * mypow((point-shift[1]),5) / 120.0 + coe[2] * mypow((point-shift[2]),5) / 120.0 + coe[3] * mypow((point-shift[3]),5) / 120.0;
	}
	else if (point <= shift[5])
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0 + coe[1] * mypow((point-shift[1]),5) / 120.0 + coe[2] * mypow((point-shift[2]),5) / 120.0 + coe[3] * mypow((point-shift[3]),5) / 120.0 + coe[4] * mypow((point-shift[4]),5) / 120.0;
	}
	else if (point <= shift[6])
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0 + coe[1] * mypow((point-shift[1]),5) / 120.0 + coe[2] * mypow((point-shift[2]),5) / 120.0 + coe[3] * mypow((point-shift[3]),5) / 120.0 + coe[4] * mypow((point-shift[4]),5) / 120.0 + coe[5] * mypow((point-shift[5]),5) / 120.0;
	}	
	else if (point <= shift[7])
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0 + coe[1] * mypow((point-shift[1]),5) / 120.0 + coe[2] * mypow((point-shift[2]),5) / 120.0 + coe[3] * mypow((point-shift[3]),5) / 120.0 + coe[4] * mypow((point-shift[4]),5) / 120.0 + coe[5] * mypow((point-shift[5]),5) / 120.0 + coe[6] * mypow((point-shift[6]),5) / 120.0;
	}
	else if (point <= shift[8])
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0 + coe[1] * mypow((point-shift[1]),5) / 120.0 + coe[2] * mypow((point-shift[2]),5) / 120.0 + coe[3] * mypow((point-shift[3]),5) / 120.0 + coe[4] * mypow((point-shift[4]),5) / 120.0 + coe[5] * mypow((point-shift[5]),5) / 120.0 + coe[6] * mypow((point-shift[6]),5) / 120.0 + coe[7] * mypow((point-shift[7]),5) / 120.0;
	}
	else if (point <= shift[9])
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0 + coe[1] * mypow((point-shift[1]),5) / 120.0 + coe[2] * mypow((point-shift[2]),5) / 120.0 + coe[3] * mypow((point-shift[3]),5) / 120.0 + coe[4] * mypow((point-shift[4]),5) / 120.0 + coe[5] * mypow((point-shift[5]),5) / 120.0 + coe[6] * mypow((point-shift[6]),5) / 120.0 + coe[7] * mypow((point-shift[7]),5) / 120.0 + coe[8] * mypow((point-shift[8]),5) / 120.0;
	}
	else
	{
		V = coe[0] * mypow((point-shift[0]),5) / 120.0 + coe[1] * mypow((point-shift[1]),5) / 120.0 + coe[2] * mypow((point-shift[2]),5) / 120.0 + coe[3] * mypow((point-shift[3]),5) / 120.0 + coe[4] * mypow((point-shift[4]),5) / 120.0 + coe[5] * mypow((point-shift[5]),5) / 120.0 + coe[6] * mypow((point-shift[6]),5) / 120.0 + coe[7] * mypow((point-shift[7]),5) / 120.0 + coe[8] * mypow((point-shift[8]),5) / 120.0 + coe[9] * mypow((point-shift[9]),5) / 120.0;
	}
	return V;
}


double interp12
(
	double point,
	double* shift,
	double* coe
	)
{
	double V = 0;
	if (point <= shift[0])
	{
		V = 0;
	}
	else if (point <= shift[1])
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0;
	}
	else if (point <= shift[2])
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0 + coe[1] * mypow((point-shift[1]),3) / 6.0;
	}
	else if (point <= shift[3])
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0 + coe[1] * mypow((point-shift[1]),3) / 6.0 + coe[2] * mypow((point-shift[2]),3) / 6.0;
	}
	else if (point <= shift[4])
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0 + coe[1] * mypow((point-shift[1]),3) / 6.0 + coe[2] * mypow((point-shift[2]),3) / 6.0 + coe[3] * mypow((point-shift[3]),3) / 6.0;
	}
	else if (point <= shift[5])
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0 + coe[1] * mypow((point-shift[1]),3) / 6.0 + coe[2] * mypow((point-shift[2]),3) / 6.0 + coe[3] * mypow((point-shift[3]),3) / 6.0 + coe[4] * mypow((point-shift[4]),3) / 6.0;
	}
	else
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0 + coe[1] * mypow((point-shift[1]),3) / 6.0 + coe[2] * mypow((point-shift[2]),3) / 6.0 + coe[3] * mypow((point-shift[3]),3) / 6.0 + coe[4] * mypow((point-shift[4]),3) / 6.0 + coe[5] * mypow((point-shift[5]),3) / 6.0;
	}	
	return V;
}


double interp8
(
	double point,
	double* shift,
	double* coe
	)
{
	double V = 0;
	if (point <= shift[0])
	{
		V = 0;
	}
	else if (point <= shift[1])
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0;
	}
	else if (point <= shift[2])
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0 + coe[1] * mypow((point-shift[1]),3) / 6.0;
	}
	else if (point <= shift[3])
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0 + coe[1] * mypow((point-shift[1]),3) / 6.0 + coe[2] * mypow((point-shift[2]),3) / 6.0;
	}
	else
	{
		V = coe[0] * mypow((point-shift[0]),3) / 6.0 + coe[1] * mypow((point-shift[1]),3) / 6.0 + coe[2] * mypow((point-shift[2]),3) / 6.0 + coe[3] * mypow((point-shift[3]),3) / 6.0;
	}
	return V;
}



double interp4
(
	double point,
	double* shift,
	double* coe
	)
{
	double V = 0;
	if (point <= shift[0])
	{
		V = 0;
	}
	else if (point <= shift[1])
	{
		V = coe[0] * (point - shift[0]);
	}
	else
	{
		V = coe[0] * (point - shift[0]) + coe[1] * (point - shift[1]);
	}
	return V;
}




double value
(
	int num_shift,
	double interp,
	double* coe,
	double* shift,
	double flag
	)
{
	double result = 0;
	if (flag != 0)
	{
		if (num_shift == 20)
		{
			result = interp20(interp,shift,coe);	
		}
		else
		{
			result = interp8(interp,shift,coe);
		}
	}
	else
	{
		if (num_shift == 12)
		{
			result = interp12(interp,shift,coe);
		}
		else
		{
			result = interp4(interp,shift,coe);
		}
	}
	return result;
}




void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){

    if (nrhs != 8)//size,1d角度数组,宽度，
        mexErrMsgTxt("Wrong number of input arguments.\n");
   // 检查输入变量数量是否正确，否则报错

    if (nlhs > 1)   
        mexErrMsgTxt("Too many output argumnents.\n");
    // 检查输出变量数量是否正确，否则报错


    const mxArray *shiftcell;//allshift
    shiftcell = prhs[0];
    double  *shift;

    int numtheta;//number of degree
	numtheta = mxGetDimensions(prhs[0])[0];
	//mexPrintf("%d sol\n", numtheta);
    
    int numshift;//20 or 12 or 8 or 4

    double *allcoe = mxGetPr(prhs[1]);//allcoe

 	const mxArray *coecell;//allshift
    coecell = prhs[1];
    double  *coe;

    double *g = mxGetPr(prhs[2]);
 	int ysize = mxGetM(prhs[2]);
 	double *offset = mxGetPr(prhs[3]);	
 	int numNN = mxGetM(prhs[3]);
 	double *ymgrid = mxGetPr(prhs[4]);
 	double *area = mxGetPr(prhs[5]);
 	double ystep = *mxGetPr(prhs[6]);
 	double detectorwidth = *mxGetPr(prhs[7]);

 	plhs[0] = mxCreateDoubleMatrix(sqrt(numNN),sqrt(numNN),mxREAL);
 	double *OUT = mxGetPr(plhs[0]);
 	
 	double width;
 	double interp;
	for (int I_theta = 0; I_theta < numtheta; I_theta++) 
	{
	    shift = mxGetPr(mxGetCell(shiftcell,I_theta));
	    numshift = mxGetM(mxGetCell(shiftcell,I_theta));
	    coe = mxGetPr(mxGetCell(coecell,I_theta));
	    width = area[I_theta];
	    for (int i = 0; i < numNN; i++) 
	    {
		    interp = offset[I_theta*numNN+i];
		    double xl = floor(((interp - width)-ymgrid[0])/ystep);//四舍五入取整，确保插值正确
			double xr = ceil(((interp + width)-ymgrid[0])/ystep);//四舍五入取整，确保插值正确
			
			if (xr > ysize -1)
			{
				xr = ysize -1;
			}
			if (xl <0)
			{
				xl = 0;
			}
			double result = 0;
			for (int j = xl; j<= xr; j++)
			{
				result += g[I_theta*ysize+j]*value(numshift, -fabs(ymgrid[j] - interp), coe, shift, detectorwidth);
			}
			OUT[i]+=result;
		}
	}	
}


