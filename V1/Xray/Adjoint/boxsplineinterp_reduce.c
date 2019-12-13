#include "mex.h"
#include "math.h"
#include <stdio.h>
const double PI  =3.141592653589793238463;

double interp8
(
	double point,
	double* shift,
	double* coe
	)
{
	double V = 0;
	if (point <=0)
	{
		V = 0;
	}
	else if (point <= shift[1])
	{
		V = coe[0] * point * point * point / 6.0;
	}
	else if (point <= shift[2])
	{
		V = coe[0] * point * point * point / 6.0 + coe[1] * (point - shift[1]) * (point - shift[1]) * (point - shift[1]) / 6.0;
	}
	else if (point <= shift[3])
	{
		V = (coe[0] * point * point * point + coe[1] * (point - shift[1]) * (point - shift[1]) * (point - shift[1]) + coe[2] * (point - shift[2]) * (point - shift[2]) * (point - shift[2])) / 6.0;
	}
	else if (point <= shift[4])
	{
		V = (coe[0] * point * point * point + coe[1] * (point - shift[1]) * (point - shift[1]) * (point - shift[1]) + coe[2] * (point - shift[2]) * (point - shift[2]) * (point - shift[2]) + coe[3] * (point - shift[3]) * (point - shift[3]) * (point - shift[3])) / 6.0;
	}
	else if (point <= shift[5])
	{
		V = (coe[0] * point * point * point + coe[1] * (point - shift[1]) * (point - shift[1]) * (point - shift[1]) + coe[2] * (point - shift[2]) * (point - shift[2]) * (point - shift[2]) + coe[3] * (point - shift[3]) * (point - shift[3]) * (point - shift[3]) + coe[4] * (point - shift[4]) * (point - shift[4]) * (point - shift[4])) / 6.0;
	}
	else if (point <= shift[6])
	{
		V = (coe[0] * point * point * point + coe[1] * (point - shift[1]) * (point - shift[1]) * (point - shift[1]) + coe[2] * (point - shift[2]) * (point - shift[2]) * (point - shift[2]) + coe[3] * (point - shift[3]) * (point - shift[3]) * (point - shift[3]) + coe[4] * (point - shift[4]) * (point - shift[4]) * (point - shift[4]) + coe[5] * (point - shift[5]) * (point - shift[5]) * (point - shift[5])) / 6.0;
	}
	else if (point <= shift[7])
	{
		V = (coe[0] * point * point * point + coe[1] * (point - shift[1]) * (point - shift[1]) * (point - shift[1]) + coe[2] * (point - shift[2]) * (point - shift[2]) * (point - shift[2]) + coe[3] * (point - shift[3]) * (point - shift[3]) * (point - shift[3]) + coe[4] * (point - shift[4]) * (point - shift[4]) * (point - shift[4]) + coe[5] * (point - shift[5]) * (point - shift[5]) * (point - shift[5]) + coe[6] * (point - shift[6]) * (point - shift[6]) * (point - shift[6])) / 6.0;
	}
	else
	{
		V = (coe[0] * shift[7] * shift[7] * shift[7] + coe[1] * (shift[7] - shift[1]) * (shift[7] - shift[1]) * (shift[7] - shift[1]) + coe[2] * (shift[7] - shift[2]) * (shift[7] - shift[2]) * (shift[7] - shift[2]) + coe[3] * (shift[7] - shift[3]) * (shift[7] - shift[3]) * (shift[7] - shift[3]) + coe[4] * (shift[7] - shift[4]) * (shift[7] - shift[4]) * (shift[7] - shift[4]) + coe[5] * (shift[7] - shift[5]) * (shift[7] - shift[5]) * (shift[7] - shift[5]) + coe[6] * (shift[7] - shift[6]) * (shift[7] - shift[6]) * (shift[7] - shift[6])) / 6.0;
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
	if (point <=0)
	{
		V = 0;
	}
	else if (point <= shift[1])
	{
		V = coe[0] * point * point / 2.0;
	}
	else if (point <= shift[2])
	{
		V = coe[0] * shift[1] * shift[1] / 2.0 + (point - shift[1]) * coe[0] * shift[1];
	}
	else if (point <= shift[3])
	{
		V = coe[0] * shift[1] * shift[1] + (shift[2] - shift[1]) * coe[0] * shift[1] - coe[0] * (shift[3] - point) * (shift[3] - point) / 2.0;
	}
	else
	{
		V = coe[0] * shift[1] * shift[1] + (shift[2] - shift[1]) * coe[0] * shift[1];
	}
	return V;
}

double interp2
(
	double point,
	double* shift,
	double* coe
	)
{
	double V = 0;
	if (point <=0)
	{
		V = 0;
	}
	else if (point <= shift[1])
	{
		V = coe[0] * point;
	}
	else
	{
		V = coe[0] * shift[1];
	}
	return V;
}


double value
(
	int num_shift,
	double left,
	double right,
	double* coe,
	double* shift,
	double flag
	)
{
	double result = 0;
	double vL = 0;
	double vR = 0;
	if (flag != 0)
	{
		if (num_shift == 8)
		{
			vL = interp8(left,shift,coe);
			vR = interp8(right,shift,coe);
			result = vR - vL;	
		}
		else
		{
			vL = interp4(left,shift,coe);
			vR = interp4(right,shift,coe);
			result = vR - vL;
		}
	}
	else
	{
		if (num_shift == 4)
		{
			vL = interp4(left,shift,coe);
			vR = interp4(right,shift,coe);
			result = vR - vL;
		}
		else
		{
			vL = interp2(left,shift,coe);
			vR = interp2(right,shift,coe);
			result = vR - vL;
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
    int numtheta;//number of degree
	numtheta = mxGetDimensions(prhs[0])[0];
	//mexPrintf("%d sol\n", numtheta);
    double  *shift;
    int numshift;//8 or 4 or 2

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
	    //mexPrintf("%f shift\n", shift[0]);
	    //mexPrintf("%f coe\n", coe[0]);
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
			//double left = 0;
			//double right = 0;
			for (int j = xl; j<= xr; j++)
			{
				//left = xl;
				//right = xl + 1.0;

				//result += g[I_theta*ysize+j]*value(numshift,-fabs(ymgrid[j]-interp),coe,shift,detectorwidth);
				result += g[I_theta*ysize+j]*value(numshift, (ymgrid[j] - 0.5*ystep - interp + width), (ymgrid[j] + 0.5*ystep - interp + width), coe, shift, detectorwidth);
			}
			OUT[i]+=result;
		}
	}	
}


