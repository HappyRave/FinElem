#include "fem.h"

double earthRadius()
{
    return 6371000.0;
}

double earthJacobian(double x[4], double y[4], double xsi, double eta, double R)
{
     double jacobian = 0.0;
     return jacobian;
}


double earthGridInterpolate(femGrid *theGrid, double x, double y)
{
	// get elements from theGrid
    float Ox=theGrid->Ox;
	float Oy=theGrid->Oy;
	float dx=theGrid->dx;
    int nx=theGrid->nx;
    int ny=theGrid->ny;
	
	// get the relatif x and y position to Ox and Oy in the elem array
	float Rx=x-Ox;
	float Ry=y-Oy;
	int Ax=floor(fabs(Rx/dx));
	int Ay=floor(fabs(Ry/dx));
	
	
	// the points for the interpolation
	int pt1[2]={Ax, Ay};
	int pt2[2]={pt1[0]+1, pt1[1]};
	int pt3[2]={pt1[0], pt1[1]+1};
	int pt4[2]={pt1[0]+1, pt1[1]+1};
	
	// check if the points are in the array
	if (pt1[0]<0 || pt1[0]>nx-2) {
		return 0;
	}
	if (pt1[1]<0 || pt1[1]>ny-2) {
		return 0;
	}
	
	// the bathymetrie for the 4 points
	float bat1= theGrid->elem[pt1[0]][pt1[1]];
	float bat2= theGrid->elem[pt2[0]][pt2[1]];
	float bat3= theGrid->elem[pt3[0]][pt3[1]];
	float bat4= theGrid->elem[pt4[0]][pt4[1]];
	
	// value of x1, x2, y1, y2
	float x1=Ox+Ax*dx;
	float x2=Ox+(Ax+1)*dx;
	float y1=Oy+Ay*dx;
	float y2=Oy+(Ay+1)*dx;
	
	// bilinear interpolation
	float bat=(bat1*(x2-x)*(y2-y)+bat2*(x-x1)*(y2-y)+bat3*(x2-x)*(y-y1)+bat4*(x-x1)*(y-y1))/((x2-x1)*(y2-y1));
	
	// if bat > 0 => not unde water
	if (bat>0) {
		return 0;
	} else {
		return bat;
	}
}

double earthIntegrateBathymetry(femMesh *theMesh, femGrid *theGrid, femIntegration *theRule, double radius)
{
    return 3.14;
}

