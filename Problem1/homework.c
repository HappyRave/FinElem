//
//  Created by Laurent Debersaques and Maxime De Mol on 1/02/13.
//  Copyright (c) 2013 DaKot. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double phi(int order, double xsi, double eta);
double varTransform(double z[4], double xsi, double eta);

double gaussJacobian(double x[4], double y[4], double xsi, double eta) {

    double jacobian;
    double j[2][2];
    
    j[0][0]= (x[0]-x[1])*(1.0+eta) + (x[3]-x[2])*(1.0-eta);
    j[0][1]= (x[0]-x[3])*(1.0+xsi) + (x[1]-x[2])*(1.0-xsi);
    j[1][0]= (y[0]-y[1])*(1.0+eta) + (y[3]-y[2])*(1.0-eta);
    j[1][1]= (y[0]-y[3])*(1.0+xsi) + (y[1]-y[2])*(1.0-xsi);
    jacobian=( j[0][0]*j[1][1] - j[0][1]*j[1][0] ) / 16.0;
    
    return jacobian;
}

double gaussIntegrate(double x[4], double y[4], double (*f) (double, double)){

    double I = 0;
    
    double val = 1/sqrt(3);
    double xsik[4]= {val, -val, -val, val};
    double etak[4]= {val, val, -val, -val};
    
    int j; double xsi, eta;
    for(j=0; j<4; j++) {
		xsi = xsik[j]; eta = etak[j];
		I += (*f)(varTransform(x,xsi,eta), varTransform(y,xsi,eta)) * gaussJacobian(x,y,xsi,eta);
    }
    
    return I;
}

double gaussIntegrateRecursive(double x[4], double y[4], double (*f)(double,double), int n){

     double I;
    
    if (n==0) {

        I=gaussIntegrate(x,y,(*f));

    } else {

		// Edge-center coordinates
        double midUpX=(x[0]+x[1])/2.0, midUpY=(y[0]+y[1])/2.0;
        double midLeftX=(x[1]+x[2])/2.0, midLeftY=(y[1]+y[2])/2.0;
        double midDownX=(x[2]+x[3])/2.0, midDownY=(y[2]+y[3])/2.0;
        double midRightX=(x[3]+x[0])/2.0, midRightY=(y[3]+y[0])/2.0;

        // Center coordinate
        double midCenterX=(midUpX+midDownX)/2.0;
        double midCenterY=(midUpY+midDownY)/2.0;
        
		// u = upper, b = bottom, R = right, L = left
        double uRX[4]={x[0], midUpX, midCenterX, midRightX};
        double uRY[4]={y[0], midUpY, midCenterY, midRightY};
        double uLX[4]={midUpX, x[1], midLeftX, midCenterX};
        double uLY[4]={midUpY, y[1], midLeftY, midCenterY};
        double bLX[4]={midCenterX, midLeftX, x[2], midDownX};
        double bLY[4]={midCenterY, midLeftY, y[2], midDownY};
        double bRX[4]={midRightX, midCenterX, midDownX, x[3]};
        double bRY[4]={midRightY, midCenterY, midDownY, y[3]};
        
        I = gaussIntegrateRecursive(uRX,uRY,f,n-1)
			+ gaussIntegrateRecursive(uLX,uLY,f,n-1)
			+ gaussIntegrateRecursive(bLX,bLY,f,n-1)
			+ gaussIntegrateRecursive(bRX,bRY,f,n-1);

    }
    
    return I;

}

double phi(int order, double xsi, double eta){

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

double varTransform(double z[4], double xsi, double eta){

    double Z = 0;

    int i;
    for (i=1; i<=4; i++) {
        Z = Z + z[i-1]*phi(i,xsi,eta);
    }

    return Z;

}
