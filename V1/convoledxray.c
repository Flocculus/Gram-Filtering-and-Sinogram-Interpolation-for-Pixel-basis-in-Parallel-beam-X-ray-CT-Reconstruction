#include "mex.h"
#include "math.h"
#include <stdio.h>
const double PI  =3.141592653589793238463;

double value_1d3dbox
(
    double* shift,
    double* coe,
    double interp,
    double xstep
)
{

    double sol = 0;
    for (int i = 0; i < 8; i++)
    {
        if (interp/xstep > shift[i])
        {
            //mexPrintf("%f shift, %f coe, %f interp\n", shift[i],coe[i],interp);
            sol += coe[i]*pow((interp/xstep-shift[i]),2)/2.0;
        }
        //mexPrintf("%f sol\n", sol);
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
double w,
double interp,
double xstep
)
{   
    double sol = 0;
    if (interp < (c+s+w)*xstep/2.0)
    {
        if ((c > 1E-14)&&(s > 1E-14))
        {
            double para = 1.0/c/s/w;
            //mexPrintf("flag %f",xstep);
            double shift[8] = {(-c-s-w)/2.0, (c-s-w)/2.0, (s-c-w)/2.0, (w-s-c)/2.0, (c+s-w)/2.0, (c+w-s)/2.0, (s+w-c)/2.0, (s+c+w)/2.0};
            double coe[8] = {1*para, -1*para, -1*para, -1*para, 1*para, 1*para, 1*para, -1*para};
            sol = xstep*value_1d3dbox(shift, coe, -interp,xstep);


        }
        else
        {   double para = 1.0/w;
            double shift[4] = {-(w+1)/2.0, (w-1)/2.0, (1-w)/2.0, (w+1)/2.0};
            double coe[4] = {1*para, -1*para, -1*para, 1*para};
            sol = xstep*value_1d2dbox(shift, coe, -interp,xstep);
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
    double xstep = *mxGetPr(prhs[3]);//缩放倍数
    double *g = mxGetPr(prhs[4]);//image
    double ystep = *mxGetPr(prhs[5]);
    //double width = *mxGetPr(prhs[6]);//探测器宽度

    plhs[0] = mxCreateDoubleMatrix(ysize,number,mxREAL);
    double *OUT = mxGetPr(plhs[0]);

    double *center;
    double distance = 0;
    center = (double *)calloc(n*n, sizeof(double));
    


    double Vcos = 0;
    double Vsin = 0;
    double length = 0;
    double L = 0;
    double R = 0;
    for(int i = 0; i <number; i++)//degree
    {
        Vcos = cos(-theta[i]/180.0*PI-PI/2.0);
        Vsin = sin(-theta[i]/180.0*PI-PI/2.0);
        length = 0.5*(fabs(Vcos)*xstep + fabs(Vsin)*xstep + ystep);
        //mexPrintf("%f cos, %f sin,\n", Vcos, Vsin);
        
        for (int j = 0; j<n; j++)//row
        {
            for (int k = 0; k<n; k++)//col
            {
                center[k*n + j] = ((0.5 - n/2.0 +k)*Vcos + (0.5 - n/2.0 +j)*Vsin)*xstep;
                //mexPrintf("%f center\n", center[j*n + k]);
            }
        }


        for (int j = 0; j < n*n; j++)
        {
            L = floor((center[j] + 8.0 - 1E-10 - length)/ystep);
            if (L < 0)
            {
                L = 0;
            }

            R = ceil((center[j] + 8.0 - 1E-10 + length)/ystep);
            if (R > ysize)
            {
                R = ysize;
            }
            /*
            if ((i == 179)&&(j == 10))
            {   
                mexPrintf("%f cos, %f sin, %f length\n", fabs(Vcos),fabs(Vsin),length);
                mexPrintf("%f l, %f r,\n", (center[j] + 8.0 - 1E-10 - length), (center[j] + 8.0 - 1E-10 + length));
                mexPrintf("%f L, %f R,\n", L, R);
                mexPrintf("%f length, %f ystep\n", length, ystep);
            }
            */
            
            for (int k = L; k <= R; k++)
            {       
                distance = fabs(center[j] - ysamplelocation[k]);
                OUT[i*ysize+k] += getvalue(fabs(Vcos), fabs(Vsin), ystep/xstep, distance,xstep) * g[j];
            }
            
        }

       
        
    }
    free(center);



}