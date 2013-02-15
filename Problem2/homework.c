//
//  Created by Laurent Debersaques and Maxime De Mol on 14/02/13.
//  Copyright (c) 2013 DaKot. All rights reserved.
//

#include "fem.h"

double phi(int order, double xsi, double eta);
double varTransform(double z[4], double xsi, double eta);

double earthRadius()
{
    return 6371000.0;
}

double earthJacobian(double x[4], double y[4], double xsi, double eta, double R)
{
	double toRad= M_PI/180.0;
	double dXdXsi= toRad*((x[0]-x[1])*(1.0+eta) + (x[3]-x[2])*(1.0-eta))/4.0;
    double dXdEta= toRad*((x[0]-x[3])*(1.0+xsi) + (x[1]-x[2])*(1.0-xsi))/4.0;
    double dYdXsi= toRad*((y[0]-y[1])*(1.0+eta) + (y[3]-y[2])*(1.0-eta))/4.0;
    double dYdEta= toRad*((y[0]-y[3])*(1.0+xsi) + (y[1]-y[2])*(1.0-xsi))/4.0;
	
	return R*R*cos(toRad*varTransform(y,xsi,eta))*(dXdXsi*dYdEta-dXdEta*dYdXsi);
}


double earthGridInterpolate(femGrid *theGrid, double x, double y)
{
	// get elements from theGrid
    float Ox=theGrid->Ox;
	float Oy=theGrid->Oy;
	float dx=theGrid->dx;
    int nx=theGrid->nx;
    int ny=theGrid->ny;
	
	// get the relative x and y position to Ox and Oy in the elem array
	float Rx=x-Ox;
	float Ry=y-Oy;
	int Ax=floor((Rx/dx));
	int Ay=floor((Ry/dx));
	
	// the points for the interpolation
	int pt1[2]={Ax+1, Ay+1};
	int pt2[2]={Ax, Ay+1};
	int pt3[2]={Ax, Ay};
	int pt4[2]={Ax+1, Ay};
	
	if(Ax >= 0 && Ax+1 <= nx-1 && Ay >= 0 && Ay+1 <= ny-1){
		// the bathymetry for the 4 points
		float bat1= theGrid->elem[pt1[0]][pt1[1]];
		float bat2= theGrid->elem[pt2[0]][pt2[1]];
		float bat3= theGrid->elem[pt3[0]][pt3[1]];
		float bat4= theGrid->elem[pt4[0]][pt4[1]];
		
		// value of x1, x2, y1, y2
		float x1=Ox+pt1[0]*dx;
		float x2=Ox+pt3[0]*dx;
		float y1=Oy+pt1[1]*dx;
		float y2=Oy+pt3[1]*dx;
		
		// bilinear interpolation
		double bat=(bat3*(x1-x)*(y1-y)+bat4*(x-x2)*(y1-y)+bat2*(x1-x)*(y-y2)+bat1*(x-x2)*(y-y2))/((x1-x2)*(y1-y2));
		
		return bat;
	}
	return 0;
}

double earthIntegrateBathymetry(femMesh *theMesh, femGrid *theGrid, femIntegration *theRule, double radius)
{
	double I = 0;
	
	/* Shortcut variables */
	// Mesh
	int *elem = theMesh->elem;;
    double *X = theMesh->X; // Longitude
    double *Y = theMesh->Y; // Lattitude
    int nElem = theMesh->nElem;
	// Integration rule
  	int nRule = theRule->n;
  	const double *xsik = theRule->xsi;
  	const double *etak = theRule->eta;
  	const double *weightk = theRule->weight;
	
	// Browsing each quadrilateral
	int quadNr = 0;
	for(quadNr = 0; quadNr<nElem ; quadNr++){
		
		// Gathering vertices composing the quad
		int* verticesID = &elem[quadNr*4];
		double xQuad[4] = {X[verticesID[0]], X[verticesID[1]], X[verticesID[2]], X[verticesID[3]]};
		double yQuad[4] = {Y[verticesID[0]], Y[verticesID[1]], Y[verticesID[2]], Y[verticesID[3]]};
		
		// Computing integration
		int iPID = 0; double xsi = 0; double eta = 0; double weight = 1;
		for(iPID = 0; iPID < nRule ; iPID++){
			xsi = xsik[iPID]; eta = etak[iPID]; weight = weightk[iPID];
			double interpolationResult = earthGridInterpolate(theGrid, varTransform(xQuad,xsi,eta), varTransform(yQuad,xsi,eta));
			if(interpolationResult < 0){
				I += weight * interpolationResult * earthJacobian(xQuad,yQuad,xsi,eta,earthRadius());
			}
		}
	}
	return -I;
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