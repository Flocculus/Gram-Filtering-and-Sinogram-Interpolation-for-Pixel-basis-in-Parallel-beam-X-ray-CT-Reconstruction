#include "mex.h"
#include "math.h"
#include <stdio.h>
const double PI  =3.141592653589793238463;

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
    if (n == 1)
    {
        return x;
    }
}

double value_1d6dbox
(
	double* shift,
	double* coe,
	double interp,
	double xstep
)
{
	double sol = 0;
	for (int i = 0; i < 13; i++)
	{
		if (interp/xstep > shift[i])
		{
			sol += coe[i]*mypow((interp/xstep-shift[i]),5)/120.0;
		}
	}
	if(sol <= 1E-14)	
	{
		sol = 0;
	}	
	return sol;
}


double value_1d4dbox
(
	double* shift,
	double* coe,
	double interp,
	double xstep
)
{
	double sol = 0;
	for (int i = 0; i < 4; i++)
	{
		//mexPrintf("%f realdis, %f shift\n", interp/xstep , shift[i]);
		if (interp/xstep > shift[i])
		{
			//mexPrintf("%f sol\n", mypow((interp/xstep-shift[i]),3)/6.0);
			sol += coe[i]*mypow((interp/xstep-shift[i]),3)/6.0;
		}
	}
	if(sol <= 1E-14)	
	{
		sol = 0;
	}	
	return sol;
}

double value_1d2dbox
(
    double* shift,
    double* coe,
    double interp,
    double xstep
)
{
    double sol = 0;
    for (int i = 0; i < 1; i++)
    {
        //mexPrintf("%f realdis, %f shift\n", interp/xstep , shift[i]);
        if (interp/xstep > shift[i])
        {
            //mexPrintf("%f sol\n", mypow((interp/xstep-shift[i]),3)/6.0);
            sol += coe[i]*mypow((interp/xstep-shift[i]),1);
        }
    }
    if(sol <= 1E-14)    
    {
        sol = 0;
    }   
    return sol;
}


double getvalue
(
double c,
double s,
double w,
double interp,
double xstep
)
{	
	double sol = 0;
    if ( w != 0)
    {
		if ((c > 1E-14)&&(s > 1E-14))
		{
			double para = 1.0/c/c/s/s/w/w;
			double shift[27] = {-fabs(-c-s-w), -fabs(-s-w), -fabs(-c-w), -fabs(-s-c), -fabs(c-s-w), -fabs(s-c-w), -fabs(w-s-c), -fabs(-w), -fabs(-s), -fabs(-c), -fabs(c-w), -fabs(c-s), -fabs(s-w), 0, -fabs(w-s), -fabs(s-c), -fabs(w-c), -fabs(c), -fabs(s), -fabs(w), -fabs(c+s-w), -fabs(c+w-s), -fabs(w+s-c), -fabs(c+s), -fabs(c+w), -fabs(s+w),-fabs(c+s+w)};
            double coe[27] = {1*para, -2*para, -2*para, -2*para, 1*para, 1*para, 1*para, 4*para, 4*para, 4*para, -2*para, -2*para, -2*para, -8*para, -2*para, -2*para, -2*para, 4*para, 4*para, 4*para, 1*para, 1*para, 1*para, -2*para, -2*para, -2*para, 1*para};
			


            sol = xstep*value_1d6dbox(shift, coe, -interp, xstep);


		}
		else
		{	
			double para = 1.0/w/w;
			double shift[9] = {-fabs(-w-1), -fabs(-w), -fabs(-w+1), -abs(-1), 0, -abs(1), -fabs(w-1), -fabs(w), -fabs(w+1)};
			double coe[9] = {1*para, -2*para, 1*para, -2*para, 4*para, -2*para, 1*para ,-2*para, 1*para};
			sol = xstep*value_1d4dbox(shift, coe, -interp, xstep);
		}
    }
    else
    {
        if ((c > 1E-14)&&(s > 1E-14))
        {
            double para = 1.0/c/c/s/s;
            double shift[9] = {-fabs(-c-s), -fabs(-s), -fabs(c-s), -fabs(-c), 0, -fabs(c), -fabs(s-c), -fabs(s), -fabs(c+s)};
            double coe[9] = {1*para, -2*para, 1*para, -2*para, 4*para, -2*para, 1*para, -2*para, 1*para};
            


            sol = xstep*value_1d4dbox(shift, coe, -interp, xstep);


        }
        else
        {   
            double shift[3] = {-abs(-1), 0, -abs(1)};
            double coe[3] = {1, -2, 1};
            sol = xstep*value_1d2dbox(shift, coe, -interp, xstep);
        }
    }

	return sol;
}




void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){

    if (nrhs != 5)//size,1d角度数组,宽度，
        mexErrMsgTxt("Wrong number of input arguments.\n");
   // 检查输入变量数量是否正确，否则报错

    if (nlhs > 1)   
        mexErrMsgTxt("Too many output argumnents.\n");
    // 检查输出变量数量是否正确，否则报错
    
    
	//int number = mxGetN(prhs[1]);//?   
 	//int n = (int)*mxGetPr(prhs[0]);
    //double *theta = mxGetPr(prhs[1]);
    //double width = *mxGetPr(prhs[2]);
    //mexPrintf("%d thetas, %d size, %f width,\n", number, n, width);
    int number = mxGetN(prhs[0]);//角度数目  
    double *theta = mxGetPr(prhs[0]);//角度数组
    int n = *mxGetPr(prhs[1]);//图片size
    double width = *mxGetPr(prhs[2]);//探测器宽度倍数
    double xstep = *mxGetPr(prhs[3]);//缩放倍数
    double ystep = *mxGetPr(prhs[4]);//ystep
    //double *g = mxGetPr(prhs[5]);//g,number col, ysize row
    //double ystep = *mxGetPr(prhs[6]);











    plhs[0] = mxCreateDoubleMatrix(1,2*n*n,mxREAL);
    double *OUT = mxGetPr(plhs[0]);
    double *center;
    double *distance;
    center = (double *)calloc(n*n, sizeof(double));
    distance = (double *)calloc(2*n*n, sizeof(double));
    double Vcos = 0;
    double Vsin = 0;
	double length = 0;
    for(int i = 0; i <number; i++)
    {
    	Vcos = cos(-theta[i]/180.0*PI);
    	Vsin = sin(-theta[i]/180.0*PI);
    	//length = (fabs(Vcos) + fabs(Vsin) + width)*xstep;
        length = (fabs(Vcos)*xstep + fabs(Vsin)*xstep + width*ystep);
    	//mexPrintf("%f cos, %f sin, %f length\n", Vcos, Vsin,length);
    	
    	for (int j = 0; j<n; j++)
    	{
    		for (int k = 0; k<n; k++)
    		{
    			center[j*n + k] = (((0.5 - n/2.0 +k)*Vcos + (0.5 - n/2.0 +j)*Vsin))*xstep;

    		}
    	}


    	for (int j = 0; j < n*n; j++)
    	{
    		distance[j] = fabs(center[0] - center[j]);
    		
    	}

    	for (int j = 0; j < n; j++)
    	{
    		for (int k = 0; k < n; k++)
    		{
    			distance[n*n + n*j + k] = fabs(center[k] - center[j*n]);
    		}
    	}

    	for(int j = 0; j < 2*n*n; j++)
    	{
    		if (distance[j] < length)
    		{
    			OUT[j] = OUT[j] + getvalue(fabs(Vcos), fabs(Vsin), width*ystep/xstep, distance[j], xstep);
    			
    		}
    	}
    }
    free(center);
    free(distance);
}


