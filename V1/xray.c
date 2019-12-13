#include "mex.h"
#include "math.h"
#include <stdio.h>
const double PI  =3.141592653589793238463;

double value_1d2dbox
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
        
        if (interp/xstep > shift[i])
        {
            //mexPrintf("%f shift, %f coe, %f interp\n", shift[i],coe[i],interp);
            sol += coe[i]*(interp/xstep-shift[i]);
        }
        //mexPrintf("%f sol\n", sol);
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
double interp,
double xstep
)
{   
    double sol = 0;
    if (interp <= (c+s)*xstep/2.0)
    {
        if ((c > 1E-14)&&(s > 1E-14))
        {
            double para = 1.0/c/s;
            //mexPrintf("flag %f",xstep);
            double shift[8] = {(-c-s)/2.0, (c-s)/2.0, (s-c)/2.0,(s+c)/2.0};
            double coe[8] = {1*para, -1*para, -1*para, 1*para};
            sol = xstep*value_1d2dbox(shift, coe, -interp,xstep);


        }
        else
        {
            double val = 0;
            //mexPrintf("flag %f \n",fabs(interp/xstep));   
            if (fabs(interp/xstep)<=0.50)
            {
                sol = xstep;
            }
        }
    }
    return sol;
}







void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]){

    if (nrhs != 6)//
        mexErrMsgTxt("Wrong number of input arguments.\n");
   // 检查输入变量数量是否正确，否则报错

    if (nlhs > 1)   
        mexErrMsgTxt("Too many output argumnents.\n");
    // 检查输出变量数量是否正确，否则报错

    int number = mxGetN(prhs[0]);//角度数目  
    double *theta = mxGetPr(prhs[0]);//角度数组
    int ysize = mxGetN(prhs[1]);//y采样数目
    double *ysamplelocation = mxGetPr(prhs[1]);//y采样地点
    int n = *mxGetPr(prhs[2]);//图片size
    //double width = *mxGetPr(prhs[3]);//探测器宽度
    double xstep = *mxGetPr(prhs[3]);//缩放倍数
    double *g = mxGetPr(prhs[4]);//image
    double ystep = *mxGetPr(prhs[5]);

    plhs[0] = mxCreateDoubleMatrix(ysize,number,mxREAL);
    double *OUT = mxGetPr(plhs[0]);

    double *center;
    double distance = 0;
    center = (double *)calloc(n*n, sizeof(double));
    


    double Vcos = 0;
    double Vsin = 0;
    double length = 0;
    int L = 0;
    int R = 0;
    for(int i = 0; i <number; i++)//degree
    {
        Vcos = cos(-theta[i]/180.0*PI-PI/2.0);
        Vsin = sin(-theta[i]/180.0*PI-PI/2.0);
        length = 0.6*(fabs(Vcos) + fabs(Vsin))*xstep;
        //mexPrintf("%f cos, %f sin,\n", Vcos, Vsin);
        
        for (int j = 0; j<n; j++)//row
        {
            for (int k = 0; k<n; k++)//col
            {
                center[k*n + j] = ((0.5 - n/2.0 +k)*Vcos + (0.5 - n/2.0 +j)*Vsin)*xstep;
            }
        }


        for (int j = 0; j < n*n; j++)
        {
            L = ceil((center[j] + 8 - 1E-10 - length)/ystep);
            if (L < 0)
            {
                L = 0;
            }

            R = floor((center[j] + 8 - 1E-10 + length)/ystep);
            if (R > ysize)
            {
                R = ysize;
            }
            //mexPrintf("%d j %d L %d R %f center %f g[j]\n", j,L,R,center[j],g[j]); 
            for (int k = L; k <= R; k++)
            {       
                distance = fabs(center[j] - ysamplelocation[k]);
                OUT[i*ysize+k] += getvalue(fabs(Vcos), fabs(Vsin),distance,xstep) * g[j];
            }
            
        }

       
        
    }
    free(center);



}