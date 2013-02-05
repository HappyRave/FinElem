//
//  Created by Laurent Debersaques and Maxime De Mol on 1/02/13.
//  Copyright (c) 2013 DaKot. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

double gaussJacobian(double x[4], double y[4], double xsi, double eta);
double gaussIntegrate(double x[4], double y[4], double (*f) (double, double));
double gaussIntegrateRecursive(double x[4], double y[4], double (*f)(double,double), int n);
void *gaussIntegrateRecursiveParallel(void *arguments);
double gaussIntegrateRecursiveSimple(double x[4], double y[4], double (*f)(double,double), int n);
double phi(int order, double xsi, double eta);
double varTransform(double z[4], double xsi, double eta);

typedef struct GaussArgs GaussArgs;
struct GaussArgs {
    double *I;
    double x[4];
    double y[4];
    double (*f) (double, double);
    int n;
};

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
    
    if (n==0) {
        
        return gaussIntegrate(x,y,f);
        
    } else {

    
    pthread_t thread1;
    pthread_t thread2;
    pthread_t thread3;
    pthread_t thread4;

    double I[4]={0,0,0,0};

    // Edge-center coordinates
    double midUpX=(x[0]+x[1])/2.0, midUpY=(y[0]+y[1])/2.0;
    double midLeftX=(x[1]+x[2])/2.0, midLeftY=(y[1]+y[2])/2.0;
    double midDownX=(x[2]+x[3])/2.0, midDownY=(y[2]+y[3])/2.0;
    double midRightX=(x[3]+x[0])/2.0, midRightY=(y[3]+y[0])/2.0;
    
    // Center coordinate
    double midCenterX=(midUpX+midDownX)/2.0;
    double midCenterY=(midUpY+midDownY)/2.0;
    
    GaussArgs uRQuadri = {&I[0],{x[0], midUpX, midCenterX, midRightX},{y[0], midUpY, midCenterY, midRightY},f,n-1};
    GaussArgs uLQuadri = {&I[1],{midUpX, x[1], midLeftX, midCenterX},{midUpY, y[1], midLeftY, midCenterY},f,n-1};
    GaussArgs bLQuadri = {&I[2],{midCenterX, midLeftX, x[2], midDownX},{midCenterY, midLeftY, y[2], midDownY},f,n-1};
    GaussArgs bRQuadri = {&I[3],{midRightX, midCenterX, midDownX, x[3]},{midRightY, midCenterY, midDownY, y[3]},f,n-1};
    
    pthread_create( &thread1, NULL, gaussIntegrateRecursiveParallel, (void*) &uRQuadri);
    pthread_create( &thread2, NULL, gaussIntegrateRecursiveParallel, (void*) &uLQuadri);
    pthread_create( &thread3, NULL, gaussIntegrateRecursiveParallel, (void*) &bLQuadri);
    pthread_create( &thread4, NULL, gaussIntegrateRecursiveParallel, (void*) &bRQuadri);
    
    pthread_join( thread1, NULL);
    pthread_join( thread2, NULL);
    pthread_join( thread3, NULL);
    pthread_join( thread4, NULL);
    
    return I[0]+I[1]+I[2]+I[3];
    }

}

void *gaussIntegrateRecursiveParallel(void *arguments)
{
    GaussArgs *pointer;
    pointer = (GaussArgs*) arguments;
    GaussArgs args=*pointer;
    double x[4]={args.x[0],args.x[1],args.x[2],args.x[3]};
    double y[4]={args.y[0],args.y[1],args.y[2],args.y[3]};
    double (*f)(double,double)=args.f;
    int n=args.n;
    
    *args.I=gaussIntegrateRecursiveSimple(x,y,f,n);
    return 0;
    
}

double gaussIntegrateRecursiveSimple(double x[4], double y[4], double (*f)(double,double), int n){
    
    double I;
    
    if (n==0) {
        
        I=gaussIntegrate(x,y,f);
        
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
        
        I = gaussIntegrateRecursiveSimple(uRX,uRY,f,n-1)
        + gaussIntegrateRecursiveSimple(uLX,uLY,f,n-1)
        + gaussIntegrateRecursiveSimple(bLX,bLY,f,n-1)
        + gaussIntegrateRecursiveSimple(bRX,bRY,f,n-1);
        
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
