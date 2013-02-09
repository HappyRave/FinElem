#include "fem.h"

double phi(int order, double xsi, double eta);
double varTransform(double z[4], double xsi, double eta);
double gaussJacobian(double x[4], double y[4], double xsi, double eta);
double surface(double x[4], double y[4]);

double earthRadius()
{
    return 6371000.0;
}

double earthJacobian(double x[4], double y[4], double xsi, double eta, double R)
{
    double t[2][2];
    
    t[0][0]= (x[0]-x[1])*(1.0+eta) + (x[3]-x[2])*(1.0-eta);
    t[0][1]= (x[0]-x[3])*(1.0+xsi) + (x[1]-x[2])*(1.0-xsi);
    t[1][0]= (y[0]-y[1])*(1.0+eta) + (y[3]-y[2])*(1.0-eta);
    t[1][1]= (y[0]-y[3])*(1.0+xsi) + (y[1]-y[2])*(1.0-xsi);
    
	double jacobianSquare=( t[0][0]*t[1][1] - t[0][1]*t[1][0] ) / 16.0;
	
	double jacobianSphere = fabs(R*R*sin(2*M_PI/360*varTransform(y,xsi,eta)));
	//return jacobianSphere;
	
	/* FAIL */
	double j[3][3];
	j[0][0]=cos(2*M_PI/360*varTransform(x,xsi,eta))*cos(2*M_PI/360*varTransform(y,xsi,eta));
	j[0][1]=R*(t[0][0])/4.0*sin(2*M_PI/360*varTransform(x,xsi,eta))*(t[1][0])/4.0*sin(2*M_PI/360*varTransform(y,xsi,eta));
	j[0][2]=R*(t[0][1])/4.0*sin(2*M_PI/360*varTransform(x,xsi,eta))*(t[1][1])/4.0*sin(2*M_PI/360*varTransform(y,xsi,eta));
	j[1][0]=sin(2*M_PI/360*varTransform(x,xsi,eta))*cos(2*M_PI/360*varTransform(y,xsi,eta));
	j[1][1]=-R*(t[0][0])/4.0*cos(2*M_PI/360*varTransform(x,xsi,eta))*(t[1][0])/4.0*sin(2*M_PI/360*varTransform(y,xsi,eta));
	j[1][2]=-R*(t[0][1])/4.0*cos(2*M_PI/360*varTransform(x,xsi,eta))*(t[1][1])/4.0*sin(2*M_PI/360*varTransform(y,xsi,eta));
	j[2][0]=cos(2*M_PI/360*varTransform(y,xsi,eta));
	j[2][1]=-R*(t[1][0])/4.0*sin(2*M_PI/360*varTransform(y,xsi,eta));
	j[2][2]=-R*(t[1][1])/4.0*sin(2*M_PI/360*varTransform(y,xsi,eta));
	return (j[0][0]*j[1][1]*j[2][2]+j[1][0]*j[2][1]*j[0][2]+j[2][0]*j[0][1]*j[1][2])-(j[2][0]*j[1][1]*j[0][2]+j[0][0]*j[2][1]*j[1][2]+j[1][0]*j[0][1]*j[2][2]);
	 
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
	double bat=(bat1*(x2-x)*(y2-y)+bat2*(x-x1)*(y2-y)+bat3*(x2-x)*(y-y1)+bat4*(x-x1)*(y-y1))/((x2-x1)*(y2-y1));
	
	// if bat > 0 => not under water
	if (bat>0) {
		return 0;
	} else {
		return bat;
	}
}

double earthIntegrateBathymetry(femMesh *theMesh, femGrid *theGrid, femIntegration *theRule, double radius)
{
    int nElem=theMesh->nElem;
	double I=0.0;
	int i=0;
	
	for (i=0; i<nElem; i=i+4) {
		// get lat and lon for the 4 corners of the quadrilateral
		int quad[4]={theMesh->elem[i], theMesh->elem[i+1], theMesh->elem[i+2], theMesh->elem[i+3]};
		double theta[4]={theMesh->X[quad[0]], theMesh->X[quad[1]], theMesh->X[quad[2]], theMesh->X[quad[3]]};
		double phi[4]={theMesh->Y[quad[0]], theMesh->Y[quad[1]], theMesh->Y[quad[2]], theMesh->Y[quad[3]]};
		printf("theta: %f %2f %3f %4f\n", theta[0], theta[1], theta[2], theta[3]);
		printf("phi: %f %2f %3f %4f\n", phi[0], phi[1], phi[2], phi[3]);
		// convert to cartesian
		double x[4];
		double y[4];
		double z[4];
		int j;
		for (j=0; j<4; j++) {
			x[j]= radius*sin((2*M_PI/360)*theta[j])*cos((2*M_PI/360)*phi[j]);
			y[j]= radius*cos((2*M_PI/360)*theta[j])*sin((2*M_PI/360)*phi[j]);
			z[j]= radius*cos((2*M_PI/360)*phi[j]);
		}
		printf("x: %f %2f %3f %4f\n", x[0], x[1], x[2], x[3]);
		printf("y: %f %2f %3f %4f\n", y[0], y[1], y[2], y[3]);
		//printf("surface == %f\n",
		
		// integrate
		int k; double xsi,eta,w,bat,Ielem=0.0;
		for (k=0; k<4; k++) {
			xsi=theRule->xsi[k];
			eta=theRule->eta[k];
			w=theRule->weight[k];
			bat=-earthGridInterpolate(theGrid, varTransform(theta,xsi,eta), varTransform(phi,xsi,eta));
			Ielem+=w*bat*earthJacobian(theta,phi,xsi,eta,radius);
			printf("Jacob %f\n", earthJacobian(theta,phi,xsi,eta,radius));
		}
		I+=Ielem;
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

double surface(double x[4], double y[4]){
	
    double S = 0;
    
    double val = 1/sqrt(3);
    double xsik[4]= {val, -val, -val, val};
    double etak[4]= {val, val, -val, -val};
    
    int j; double xsi, eta;
    for(j=0; j<4; j++) {
		xsi = xsik[j]; eta = etak[j];
		S += gaussJacobian(x,y,xsi,eta);
    }
    
    return S;
}

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
