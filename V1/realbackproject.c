#include "mex.h"
#include "math.h"
#include <stdio.h>
const double PI  =3.141592653589793238463;


double value_1d6dboxccssDD
(
    double* shift,
    double* coe,
    double interp,
    double xstep
)
{
    double sol = 0;
    for (int i = 0; i < 27; i++)
    {   
        
        if (interp/xstep > shift[i])
        {
            //mexPrintf("%f shift, %f coe, %f interp\n", shift[i],coe[i],interp);
            sol += coe[i]*(interp/xstep-shift[i])*(interp/xstep-shift[i])*(interp/xstep-shift[i])*(interp/xstep-shift[i])*(interp/xstep-shift[i])/120.0;
        }
        //mexPrintf("%f sol\n", sol);
    }
    if(sol <= 1E-14)    
    {
        sol = 0;
    }   
    return sol;
}


double value_1d4dboxSSDD
(
    double* shift,
    double* coe,
    double interp,
    double xstep
)
{
    double sol = 0;
    for (int i = 0; i < 9; i++)
    {   
        
        if (interp/xstep > shift[i])
        {
            //mexPrintf("%f shift, %f coe, %f interp\n", shift[i],coe[i],interp);
            sol += coe[i]*(interp/xstep-shift[i])*(interp/xstep-shift[i])*(interp/xstep-shift[i])/6.0;
        }
        //mexPrintf("%f sol\n", sol);
    }
    if(sol <= 1E-14)    
    {
        sol = 0;
    }   
    return sol;
}

double value_1d4dboxccss
(
    double* shift,
    double* coe,
    double interp,
    double xstep
)
{
    double sol = 0;
    for (int i = 0; i < 9; i++)
    {   
        
        if (interp/xstep > shift[i])
        {
            //mexPrintf("%f shift, %f coe, %f interp\n", shift[i],coe[i],interp);
            sol += coe[i]*(interp/xstep-shift[i])*(interp/xstep-shift[i])*(interp/xstep-shift[i])/6.0;
        }
        //mexPrintf("%f sol\n", sol);
    }
    if(sol <= 1E-14)    
    {
        sol = 0;
    }   
    return sol;
}


double value_1d2dboxSS
(
    double* shift,
    double* coe,
    double interp,
    double xstep
)
{
    double sol = 0;
    for (int i = 0; i < 3; i++)
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
double xstep,
double length,
int flag,
double detector
)
{   
    double sol = 0;
    if (interp <= length)
    {
        if (flag == 0)
        {
            if ((c > 1E-14)&&(s > 1E-14))//4d ccss
            {
                double para = 1.0/c/s/c/s;
                double shift[9] = {(-c-s),(-s),(-c),(0.0),(c-s),(s-c),(c),(s),(c+s)};
                double coe[9] = {1*para, -2*para, -2*para, 4*para, 1*para, 1*para, -2*para, -2*para, 1*para};
                sol = xstep*value_1d4dboxccss(shift, coe, -interp,xstep);


            }
            else//2d step step
            {
                double para = 1.0;
                double shift[3] = {-1.0, 0, 1.0};
                double coe[3] = {1*para, -2.0*para, -1*para};
                sol = xstep*value_1d2dboxSS(shift, coe, -interp,xstep);
            }
        }
        else
        {
            if ((c > 1E-14)&&(s > 1E-14))//6d ccssDD
            {
                double para = 1.0*xstep*xstep/c/s/c/s/detector/detector;
                double shift[27] = {(-c-s-detector/xstep),(-s-detector/xstep),(-c-detector/xstep),(0.0-detector/xstep),(c-s-detector/xstep),(s-c-detector/xstep),(c-detector/xstep),(s-detector/xstep),(c+s-detector/xstep),   (-c-s),(-s),(-c),(0.0),(c-s),(s-c),(c),(s),(c+s),   (-c-s+detector/xstep),(-s+detector/xstep),(-c+detector/xstep),(0.0+detector/xstep),(c-s+detector/xstep),(s-c+detector/xstep),(c+detector/xstep),(s+detector/xstep),(c+s+detector/xstep)};
                double coe[27] = {1*para, -2*para, -2*para, 4*para, 1*para, 1*para, -2*para, -2*para, 1*para,       -2*para, 4*para, 4*para, -8*para, -2*para, -2*para, 4*para, 4*para, -2*para,     1*para, -2*para, -2*para, 4*para, 1*para, 1*para, -2*para, -2*para, 1*para};
                sol = xstep*value_1d6dboxccssDD(shift, coe, -interp,xstep);


            }
            else//4d SSDD
            {
                double para = 1.0*xstep*xstep/detector/detector;
                double shift[9] = {-detector/xstep-1, -detector/xstep, 1-detector/xstep, -1, 0, 1, detector/xstep-1, detector/xstep, detector/xstep+1};
                double coe[9] = {1*para, -2*para, 1*para, -2*para, 4*para, -2*para, 1*para, -2*para, 1*para};
                sol = xstep*value_1d4dboxSSDD(shift, coe, -interp,xstep);
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
    int n = *mxGetPr(prhs[1]);//图片size
    double xstep = *mxGetPr(prhs[2]);//缩放倍数
    double *g = mxGetPr(prhs[3]);//image
    int flag = *mxGetPr(prhs[4]);
    double detector = *mxGetPr(prhs[5]);

    plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);
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
        if (flag == 0)
        {
            length = 0.5*(fabs(Vcos) + fabs(Vsin) + fabs(Vcos) + fabs(Vsin))*xstep;
        }
        else
        {
            length = 0.5*(fabs(Vcos) + fabs(Vsin) + fabs(Vcos) + fabs(Vsin) + 2*detector/xstep)*xstep;
        }
        
        //mexPrintf("%d i, %f cos, %f sin,\n",i, Vcos, Vsin);
        
        for (int j = 0; j<n; j++)//row
        {
            for (int k = 0; k<n; k++)//col
            {
                center[k*n + j] = ((0.5 - n/2.0 +k)*Vcos + (0.5 - n/2.0 +j)*Vsin)*xstep;
            }
        }


        for (int j = 0; j < n*n; j++)
        {
            /*
            L = ceil((center[j] + 8 - 1E-10 - length)/xstep);
            if (L < 0)
            {
                L = 0;
            }

            R = floor((center[j] + 8 - 1E-10 + length)/xstep);
            if (R > ysize)
            {
                R = ysize;
            }
            */
            //mexPrintf("%d j %d L %d R %f center %f g[j]\n", j,L,R,center[j],g[j]); 
            for (int k = 0; k <= n*n; k++)
            {       
                distance = fabs(center[j] - center[k]);
                OUT[j] += getvalue(fabs(Vcos), fabs(Vsin),distance,xstep,length,flag,detector) * g[k];
            }
            
        }

       
        
    }
    free(center);



}