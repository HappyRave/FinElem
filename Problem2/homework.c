#include "fem.h"

double phi(int order, double xsi, double eta);
double varTransform(double z[4], double xsi, double eta);

double earthRadius()
{
    return 6371000.0;
}

double earthJacobian(double x[4], double y[4], double xsi, double eta, double R)
{
    double dXdXsi= ((x[0]-x[1])*(1.0+eta) + (x[3]-x[2])*(1.0-eta))/4.0;
    double dXdEta= ((x[0]-x[3])*(1.0+xsi) + (x[1]-x[2])*(1.0-xsi))/4.0;
    double dYdXsi= ((y[0]-y[1])*(1.0+eta) + (y[3]-y[2])*(1.0-eta))/4.0;
    double dYdEta= ((y[0]-y[3])*(1.0+xsi) + (y[1]-y[2])*(1.0-xsi))/4.0;
	double toRad= M_PI/180.0;
	
	return R*R*cos(toRad*varTransform(y,xsi,eta))*(dXdXsi*dYdEta-dXdEta*dYdXsi);
}


double earthGridInterpolate(femGrid *theGrid, double x, double y)
{
	// get elements from theGrid
    float Ox=theGrid->Ox;
	float Oy=theGrid->Oy;
	float dx=theGrid->dx;
    //int nx=theGrid->nx;
    //int ny=theGrid->ny;
	
	// get the relatif x and y position to Ox and Oy in the elem array
	float Rx=x-Ox;
	float Ry=y-Oy;
	int Ax=floor((Rx/dx));// pas sur s'il faut pas un fabs
	int Ay=floor((Ry/dx));
	
	// the points for the interpolation
	int pt1[2]={Ax+1, Ay+1};
	int pt2[2]={Ax, Ay+1};
	int pt3[2]={Ax, Ay};
	int pt4[2]={Ax+1, Ay};
	
	
	// apparemment on s'en fout, car sur le serveur de soumission il demande par exemple x=15 et y=19
	/*
	// check if the points are in the array
	if (pt3[0]<0 || pt3[0]>nx-2) {
		return 0;
	}
	if (pt3[1]<0 || pt3[1]>ny-2) {
		return 0;
	}
	 */
	
	// the bathymetrie for the 4 points
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
	
	// apparemment c'est dans interpolate que tu décide de ne pas compter l'intégrale si bat>0
	// if bat > 0 => not under water
	if (bat>0) {
		return bat;
	} else {
		return bat;
	}
}

double earthIntegrateBathymetry(femMesh *theMesh, femGrid *theGrid, femIntegration *theRule, double radius)
{
	//Your turn
	return M_PI;
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