//
//  Created by Laurent Debersaques and Maxime De Mol on 1/02/13.
//  Copyright (c) 2013 DaKot. All rights reserved.
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double phi(int order, double xsi, double eta);
double phix(double x[4], double y[4], double xsi, double eta);
double phiy(double x[4], double y[4], double xsi, double eta);

double gaussJacobian(double x[4], double y[4], double xsi, double eta)
{
    double jacobian;
    double j[2][2];
    
    j[0][0]=(x[0]*(1+eta)*(1.0/4.0))-(x[1]*(1+eta)*(1.0/4.0))-(x[2]*(1-eta)*(1.0/4.0))+(x[3]*(1-eta)*(1.0/4.0));
    j[0][1]=(x[0]*(1+xsi)*(1.0/4.0))+(x[1]*(1-xsi)*(1.0/4.0))-(x[2]*(1-xsi)*(1.0/4.0))-(x[3]*(1+xsi)*(1.0/4.0));
    j[1][0]=(y[0]*(1+eta)*(1.0/4.0))-(y[1]*(1+eta)*(1.0/4.0))-(y[2]*(1-eta)*(1.0/4.0))+(y[3]*(1-eta)*(1.0/4.0));
    j[1][1]=(y[0]*(1+xsi)*(1.0/4.0))+(y[1]*(1-xsi)*(1.0/4.0))-(y[2]*(1-xsi)*(1.0/4.0))-(y[3]*(1+xsi)*(1.0/4.0));
    jacobian=(j[0][0]*j[1][1])-(j[0][1]*j[1][0]);
    
    return jacobian;
}


double gaussIntegrate(double x[4], double y[4], double (*f) (double, double))
{
    double I;
    
    double val=1/sqrt(3);
    double xsik[4]={ val, -val, -val, val};
    double etak[4]={ val, val, -val, -val};
    
    double i[4];
    int j;
    
    for (j=0; j<4; j++) {
        i[j]=f(phix(x,y,xsik[j],etak[j]),phiy(x,y,xsik[j],etak[j]))*gaussJacobian(x,y,xsik[j],etak[j]);
    }
    
    I=i[0]+i[1]+i[2]+i[3];
    
    return I;
}

double gaussIntegrateRecursive(double x[4], double y[4], double (*f)(double,double), int n)
{
     double I;
    
    if (n==0) {
        I=gaussIntegrate(x,y,f);
    }
    else {
        double midUpX=(x[0]+x[1])/2.0;
        double midUpY=(y[0]+y[1])/2.0;
        double midLeftX=(x[1]+x[2])/2.0;
        double midLeftY=(y[1]+y[2])/2.0;
        double midDownX=(x[2]+x[3])/2.0;
        double midDownY=(y[2]+y[3])/2.0;
        double midRightX=(x[3]+x[0])/2.0;
        double midRightY=(y[3]+y[0])/2.0;
        // Vu que les mmedianes se coupent toujours au milieu
        double midCenterX=(midUpX+midDownX)/2.0;
        double midCenterY=(midUpY+midDownY)/2.0;
        
        double xNorth[4]={x[0], midUpX, midCenterX, midRightX};
        double yNorth[4]={y[0], midUpY, midCenterY, midRightY};
        double xWest[4]={midUpX, x[1], midLeftX, midCenterX};
        double yWest[4]={midUpY, y[1], midLeftY, midCenterY};
        double xSouth[4]={midCenterX, midLeftX, x[2], midDownX};
        double ySouth[4]={midCenterY, midLeftY, y[2], midDownY};
        double xEst[4]={midRightX, midCenterX, midDownX, x[3]};
        double yEst[4]={midRightY, midCenterY, midDownY, y[3]};
        
        I=gaussIntegrateRecursive(xNorth,yNorth,f,n-1)+gaussIntegrateRecursive(xWest,yWest,f,n-1)+gaussIntegrateRecursive(xSouth,ySouth,f,n-1)+gaussIntegrateRecursive(xEst,yEst,f,n-1);
    }
    
     return I;
}

double phi(int order, double xsi, double eta)
{
    double p;
    switch (order) {
        case 1:
            p=(1+xsi)*(1+eta)*(1.0/4.0);
            break;
        case 2:
            p=(1-xsi)*(1+eta)*(1.0/4.0);
            break;
        case 3:
            p=(1-xsi)*(1-eta)*(1.0/4.0);
            break;
        case 4:
            p=(1+xsi)*(1-eta)*(1.0/4.0);
            break;
            
        default:
            printf("Syntax Error In Function Phi");
            exit(EXIT_FAILURE);
            break;
    }
    return p;
}

double phix(double x[4], double y[4], double xsi, double eta)
{
    double X=0;
    int i;
    
    for (i=1; i<=4; i++) {
        X=X+x[i-1]*phi(i,xsi,eta);
    }

    return X;
}

double phiy(double x[4], double y[4], double xsi, double eta)
{
    double Y=0;
    int i;
    
    for (i=1; i<=4; i++) {
        Y=Y+y[i-1]*phi(i,xsi,eta);
    }
    
    return Y;
}