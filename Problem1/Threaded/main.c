#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


double gaussJacobian(double x[4], double y[4], double xsi, double eta);
double gaussIntegrate(double x[4], double y[4], double(*f)(double,double));
double gaussIntegrateRecursive(double x[4], double y[4], double(*f)(double,double), int n);

double fun(double x, double y)         { return cos(x) + y * y; }
double stupid(double x, double y)      { return 1.0; }

int main()
{
   double x[4] = { 2.0, -1.0, -1.0,  1.0} ;
   double y[4] = { 2.0,  1.0, -1.0, -1.0} ;
   int n;
   printf("Jacobian value in (0.0,0.0)  : %14.7e \n",gaussJacobian(x,y,0.0,0.0));
   printf("Jacobian value in (0.5,0.5)  : %14.7e \n",gaussJacobian(x,y,0.5,0.5));
   printf("Surface integration by Gauss-Legendre  : %14.7e \n",gaussIntegrate(x,y,stupid));
   printf("More funny integration         : %14.7e \n",gaussIntegrate(x,y,fun));
    clock_t start=clock();
    for (n=0;  n <= 12; n++){
        printf("Recursive integration (n = %2d) : %14.7e \n",n,gaussIntegrateRecursive(x,y,fun,n));}
    clock_t stop=clock();
    double time_spent = (double)(stop - start) / CLOCKS_PER_SEC;
    printf("Execution time for the recursion : %f\n", time_spent);
   return 0;
}