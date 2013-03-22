#include"fem.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef struct sparseEl sparseEl;
struct sparseEl{
	int nodeID;
	int nodeCloserToMin;
	int nodeCloserToMax;
	sparseEl* next;
};

int searchForIt(sparseEl** sparseMatrix, int minNode, int maxNode, char doYouWantCloserToMin);
double xP(femDiffusionProblem *theProblem, int mapID);
double yP(femDiffusionProblem *theProblem, int mapID);

void q3c0_x(double *xsi, double *eta)
{
    xsi[0] =  1.0;       eta[0] =  1.0;
    xsi[1] = -1.0;       eta[1] =  1.0;
    xsi[2] = -1.0;       eta[2] = -1.0;
    xsi[3] =  1.0;       eta[3] = -1.0;
    xsi[4] =  2.0/3.0;   eta[4] =  1.0;
    xsi[5] = -2.0/3.0;   eta[5] =  1.0;
    xsi[6] = -1.0;       eta[6] =  2.0/3.0;
    xsi[7] = -1.0;       eta[7] = -2.0/3.0;
    xsi[8] = -2.0/3.0;   eta[8] = -1.0;
    xsi[9] =  2.0/3.0;   eta[9] = -1.0;
    xsi[10] =  1.0;      eta[10] = -2.0/3.0;
    xsi[11] =  1.0;      eta[11] =  2.0/3.0;
    xsi[12] =  2.0/3.0;  eta[12] =  2.0/3.0;
    xsi[13] = -2.0/3.0;  eta[13] =  2.0/3.0;
    xsi[14] = -2.0/3.0;  eta[14] = -2.0/3.0;
    xsi[15] =  2.0/3.0;  eta[15] = -2.0/3.0;
}

void q3c0_phi(double xsi, double eta, double *phi)
{
	phi[0] = (1.0+xsi) * (2.0/3.0+xsi) * (-2.0/3.0+xsi) * (1.0+eta) * (2.0/3.0+eta) * (-2.0/3.0+eta)*81.0/100.0;
	phi[1] = (1.0-xsi) * (2.0/3.0-xsi) * (-2.0/3.0-xsi) * (1.0+eta) * (2.0/3.0+eta) * (-2.0/3.0+eta)*81.0/100.0;
	phi[2] = (1.0-xsi) * (2.0/3.0-xsi) * (-2.0/3.0-xsi) * (1.0-eta) * (2.0/3.0-eta) * (-2.0/3.0-eta)*81.0/100.0;
	phi[3] = (1.0+xsi) * (2.0/3.0+xsi) * (-2.0/3.0+xsi) * (1.0-eta) * (2.0/3.0-eta) * (-2.0/3.0-eta)*81.0/100.0;
	phi[4] = (1.0-xsi) * (2.0/3.0+xsi) * (1.0+xsi) * (1.0+eta) * (2.0/3.0+eta) * (-2.0/3.0+eta)*243.0/200.0;
	phi[5] = (1.0-xsi) * (2.0/3.0-xsi) * (1.0+xsi) * (1.0+eta) * (2.0/3.0+eta) * (-2.0/3.0+eta)*243.0/200.0;
	phi[6] = (1.0-xsi) * (2.0/3.0-xsi) * (-2.0/3.0-xsi) * (1.0-eta) * (2.0/3.0+eta) * (1.0+eta)*243.0/200.0;
	phi[7] = (1.0-xsi) * (2.0/3.0-xsi) * (-2.0/3.0-xsi) * (1.0-eta) * (2.0/3.0-eta) * (1.0+eta)*243.0/200.0;
	phi[8] = (1.0-xsi) * (2.0/3.0-xsi) * (1.0+xsi) * (1.0-eta) * (2.0/3.0-eta) * (-2.0/3.0-eta)*243.0/200.0;
	phi[9] = (1.0-xsi) * (2.0/3.0+xsi) * (1.0+xsi) * (1.0-eta) * (2.0/3.0-eta) * (-2.0/3.0-eta)*243.0/200.0;
	phi[10] = (1.0+xsi) * (2.0/3.0+xsi) * (-2.0/3.0+xsi) * (1.0-eta) * (2.0/3.0-eta) * (1.0+eta)*243.0/200.0;
	phi[11] = (1.0+xsi) * (2.0/3.0+xsi) * (-2.0/3.0+xsi) * (1.0-eta) * (2.0/3.0+eta) * (1.0+eta)*243.0/200.0;
	phi[12] = (1.0-xsi) * (2.0/3.0+xsi) * (1.0+xsi) * (1.0-eta) * (2.0/3.0+eta) * (1.0+eta)*729.0/400;
	phi[13] = (1.0-xsi) * (2.0/3.0-xsi) * (1.0+xsi) * (1.0-eta) * (2.0/3.0+eta) * (1.0+eta)*729.0/400;
	phi[14] = (1.0-xsi) * (2.0/3.0-xsi) * (1.0+xsi) * (1.0-eta) * (2.0/3.0-eta) * (1.0+eta)*729.0/400;
	phi[15] = (1.0-xsi) * (2.0/3.0+xsi) * (1.0+xsi) * (1.0-eta) * (2.0/3.0-eta) * (1.0+eta)*729.0/400;
	
}
void q3c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
	dphidxsi[0] = (81.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0))/100.0 + (81.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi + 2.0/3.0))/100.0 + (81.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0;
	dphidxsi[1] = - (81.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0))/100.0 - (81.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 2.0/3.0))/100.0 - (81.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0;
	dphidxsi[2] = (81.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0))/100.0 + (81.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 2.0/3.0))/100.0 + (81.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0;
	dphidxsi[3] = - (81.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0))/100.0 - (81.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi + 2.0/3.0))/100.0 - (81.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0;
	dphidxsi[4] = - (243.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1))/200.0 - (243.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 2.0/3.0))/200.0 - (243.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi + 2.0/3.0))/200.0;
	dphidxsi[5] = (243.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1))/200.0 + (243.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0))/200.0 + (243.0*(eta + 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0))/200.0;
	dphidxsi[6] = (243.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0))/200.0 + (243.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 2.0/3.0))/200.0 + (243.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0;
	dphidxsi[7] = - (243.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0))/200.0 - (243.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 2.0/3.0))/200.0 - (243.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0;
	dphidxsi[8] = - (243.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1))/200.0 - (243.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0))/200.0 - (243.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0))/200.0;
	dphidxsi[9] = (243.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1))/200.0 + (243.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 2.0/3.0))/200.0 + (243.0*(eta - 1)*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi + 2.0/3.0))/200.0;
	dphidxsi[10] = (243.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0))/200.0 + (243.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi + 1)*(xsi + 2.0/3.0))/200.0 + (243.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0;
	dphidxsi[11] = - (243.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0))/200.0 - (243.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi + 1)*(xsi + 2.0/3.0))/200.0 - (243.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0;
	dphidxsi[12] = (729.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1))/400.0 + (729.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 2.0/3.0))/400.0 + (729.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi + 1)*(xsi + 2.0/3.0))/400.0;
	dphidxsi[13] = - (729.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1))/400.0 - (729.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0))/400.0 - (729.0*(eta - 1)*(eta + 1)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0))/400.0;
	dphidxsi[14] = (729.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1))/400.0 + (729.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0))/400.0 + (729.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0))/400.0;
	dphidxsi[15] = - (729.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1))/400.0 - (729.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 2.0/3.0))/400.0 - (729.0*(eta - 1)*(eta + 1)*(eta - 2.0/3.0)*(xsi + 1)*(xsi + 2.0/3.0))/400.0;
	
	dphideta[0] = (81.0*(eta + 1)*(eta - 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0 + (81.0*(eta + 1)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0 + (81.0*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0;
	dphideta[1] = - (81.0*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0 - (81.0*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0 - (81.0*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0;
	dphideta[2] = (81.0*(eta - 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0 + (81.0*(eta - 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0 + (81.0*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0;
	dphideta[3] = - (81.0*(eta - 1)*(eta - 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0 - (81.0*(eta - 1)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0 - (81.0*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/100.0;
	dphideta[4] = - (243.0*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/200.0 - (243.0*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/200.0 - (243.0*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/200.0;
	dphideta[5] = (243.0*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/200.0 + (243.0*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/200.0 + (243.0*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/200.0;
	dphideta[6] = (243.0*(eta - 1)*(eta + 1)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0 + (243.0*(eta - 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0 + (243.0*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0;
	dphideta[7] = - (243.0*(eta - 1)*(eta + 1)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0 - (243.0*(eta - 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0 - (243.0*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0;
	dphideta[8] = - (243.0*(eta - 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/200.0 - (243.0*(eta - 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/200.0 - (243.0*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/200.0;
	dphideta[9] = (243.0*(eta - 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/200.0 + (243.0*(eta - 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/200.0 + (243.0*(eta - 2.0/3.0)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/200.0;
	dphideta[10] = (243.0*(eta - 1)*(eta + 1)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0 + (243.0*(eta - 1)*(eta - 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0 + (243.0*(eta + 1)*(eta - 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0;
	dphideta[11] = - (243.0*(eta - 1)*(eta + 1)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0 - (243.0*(eta - 1)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0 - (243.0*(eta + 1)*(eta + 2.0/3.0)*(xsi + 1)*(xsi - 2.0/3.0)*(xsi + 2.0/3.0))/200.0;
	dphideta[12] = (729.0*(eta - 1)*(eta + 1)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/400.0 + (729.0*(eta - 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/400.0 + (729.0*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/400.0;
	dphideta[13] = - (729.0*(eta - 1)*(eta + 1)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/400.0 - (729.0*(eta - 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/400.0 - (729.0*(eta + 1)*(eta + 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/400.0;
	dphideta[14] = (729.0*(eta - 1)*(eta + 1)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/400.0 + (729.0*(eta - 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/400.0 + (729.0*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi - 2.0/3.0))/400.0;
	dphideta[15] = - (729.0*(eta - 1)*(eta + 1)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/400.0 - (729.0*(eta - 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/400.0 - (729.0*(eta + 1)*(eta - 2.0/3.0)*(xsi - 1)*(xsi + 1)*(xsi + 2.0/3.0))/400.0;
}

void p3c0_x(double *xsi, double *eta)
{
    xsi[0] = 0.0;     eta[0] = 0.0;
    xsi[1] = 1.0;     eta[1] = 0.0;
    xsi[2] = 0.0;     eta[2] = 1.0;
    xsi[3] = 1.0/3.0; eta[3] = 0.0;
    xsi[4] = 2.0/3.0; eta[4] = 0.0;
    xsi[5] = 2.0/3.0; eta[5] = 1.0/3.0;
    xsi[6] = 1.0/3.0; eta[6] = 2.0/3.0;
    xsi[7] = 0.0;     eta[7] = 2.0/3.0;
    xsi[8] = 0.0;     eta[8] = 1.0/3.0;
    xsi[9] = 1.0/3.0; eta[9] = 1.0/3.0;
}

void p3c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 - xsi - eta) * (1.0/3.0 - xsi - eta) * (2.0/3.0 - xsi - eta) *  9.0/2.0;
    phi[1] =               xsi * (xsi - 1.0/3.0)       * (xsi - 2.0/3.0)       *  9.0/2.0;
    phi[2] =               eta * (eta - 1.0/3.0)       * (eta - 2.0/3.0)       *  9.0/2.0;
    phi[3] = (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) * xsi                   * 27.0/2.0;
    phi[4] = (1.0 - xsi - eta) * (xsi - 1.0/3.0)       * xsi                   * 27.0/2.0;
    phi[5] =               xsi * eta                   * (xsi - 1.0/3.0)       * 27.0/2.0;
    phi[6] =               xsi * eta                   * (eta - 1.0/3.0)       * 27.0/2.0;
    phi[7] = (1.0 - xsi - eta) * (eta - 1.0/3.0)       * eta                   * 27.0/2.0;
    phi[8] = (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) * eta                   * 27.0/2.0;
    phi[9] = (1.0 - xsi - eta) * xsi                   * eta                   * 27.0;
}

void p3c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta) {
	dphidxsi[0] = - (9.0*(eta + xsi - 1)*(eta + xsi - 1.0/3.0))/2.0 - (9.0*(eta + xsi - 1)*(eta + xsi - 2.0/3.0))/2.0 - (9.0*(eta + xsi - 1.0/3.0)*(eta + xsi - 2.0/3.0))/2.0;
	dphidxsi[1] = (9.0*(xsi - 1.0/3.0)*(xsi - 2.0/3.0))/2.0 + (9.0*xsi*(xsi - 1.0/3.0))/2.0 + (9.0*xsi*(xsi - 2.0/3.0))/2.0;
	dphidxsi[2] = 0;
	dphidxsi[3] = (27.0*xsi*(eta + xsi - 1))/2.0 + (27.0*xsi*(eta + xsi - 2.0/3.0))/2.0 + (27.0*(eta + xsi - 1)*(eta + xsi - 2.0/3.0))/2.0;
	dphidxsi[4] = - (27.0*xsi*(eta + xsi - 1))/2.0 - (27.0*(xsi - 1.0/3.0)*(eta + xsi - 1))/2.0 - (27.0*xsi*(xsi - 1.0/3.0))/2.0;
	dphidxsi[5] = (27.0*eta*xsi)/2.0 + (27.0*eta*(xsi - 1.0/3.0))/2.0;
	dphidxsi[6] = (27.0*eta*(eta - 1.0/3.0))/2.0;
	dphidxsi[7] = -(27.0*eta*(eta - 1.0/3.0))/2.0;
	dphidxsi[8] = (27.0*eta*(eta + xsi - 1))/2.0 + (27.0*eta*(eta + xsi - 2.0/3.0))/2.0;
	dphidxsi[9] = - 27.0*eta*(eta + xsi - 1) - 27.0*eta*xsi;
	dphideta[0] = - (9.0*(eta + xsi - 1)*(eta + xsi - 1.0/3.0))/2.0 - (9.0*(eta + xsi - 1)*(eta + xsi - 2.0/3.0))/2.0 - (9.0*(eta + xsi - 1.0/3.0)*(eta + xsi - 2.0/3.0))/2.0;
	dphideta[1] = 0;
	dphideta[2] = (9.0*eta*(eta - 1.0/3.0))/2.0 + (9.0*eta*(eta - 2.0/3.0))/2.0 + (9.0*(eta - 1.0/3.0)*(eta - 2.0/3.0))/2.0;
	dphideta[3] = (27.0*xsi*(eta + xsi - 1))/2.0 + (27.0*xsi*(eta + xsi - 2.0/3.0))/2.0;
	dphideta[4] = -(27.0*xsi*(xsi - 1.0/3.0))/2.0;
	dphideta[5] = (27.0*xsi*(xsi - 1.0/3.0))/2.0;
	dphideta[6] = (27.0*eta*xsi)/2.0 + (27.0*xsi*(eta - 1.0/3.0))/2.0;
	dphideta[7] = - (27.0*eta*(eta + xsi - 1))/2.0 - (27.0*(eta - 1.0/3.0)*(eta + xsi - 1))/2.0 - (27.0*eta*(eta - 1.0/3.0))/2.0;
	dphideta[8] = (27.0*eta*(eta + xsi - 1))/2.0 + (27.0*eta*(eta + xsi - 2.0/3.0))/2.0 + (27.0*(eta + xsi - 1)*(eta + xsi - 2.0/3.0))/2.0;
	dphideta[9] = - 27.0*xsi*(eta + xsi - 1) - 27.0*eta*xsi;
}




void femDiffusionComputeMap(femDiffusionProblem *theProblem, int orderType)
{
    femMesh  *theMesh = theProblem->mesh;
	
    int i,j;
    int nElem = theMesh->nElem;
    int nNode = theMesh->nNode;
    int nLocal = theProblem->space->n; // Number of nodes by element (3, 6, 10, 4, 9 or 16)
    int nLocalNode = theProblem->mesh->nLocalNode; // Number of vertices by element (3 for triangles, 4 for quads)
	
    if (orderType >= 1) {
        for (i=0; i < nNode; i++) {
            theProblem->X[i] = theMesh->X[i];
            theProblem->Y[i] = theMesh->Y[i];
		}
        for (i=0; i < nElem; i++){
			//printf("---------------------------------%d\n", theProblem->edges->nEdge);
            for (j=0; j < nLocalNode; j++){
				//printf("Element i = %d, node %d assigned to unknown %d\n",i,i*nLocalNode+j,theProblem->mesh->elem[i*nLocalNode+j]);
                theProblem->map[i*nLocal+j] = theProblem->mesh->elem[i*nLocalNode+j];
			}
		}
	}
	
	//
	// A modifier :-)
	//
	
	if(orderType > 1){
		
		/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		 * 1 - Support for creating node(s) on the edges of the mesh + update of X and Y for those nodes
		 * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		 */
		
		int variableNr = nNode;
		
		/*
		 * 1.1 Setting up a sparse matrix to be able to recover variables from the edges' endpoints
		 */
		
		// Shortcut variables
		femEdge* edges = theProblem->edges->edges;
		int nEdge = theProblem->edges->nEdge;
		
		// Decl and init of the sparse matrix
		sparseEl** sparseMatrix = malloc(sizeof(sparseEl*)*nNode);
		for(i = 0; i < nNode; i++){
			sparseMatrix[i] = NULL;
		}
		
		// Filling the sparse matrix
		sparseEl* current;
		for(i = 0; i < nEdge; i++){
			
			femEdge currentEdge = edges[i];
			
			// Normalized representation : edge described by (minNode,maxNode)
			int minNode = MIN(currentEdge.node[0],currentEdge.node[1]);
			int maxNode = MAX(currentEdge.node[0],currentEdge.node[1]);
			
			// Creating the edge to add (represented by a sparseEl) to the sparse matrix
			sparseEl* newSparseEl = malloc(sizeof(sparseEl));
			newSparseEl->nodeID = maxNode;
			newSparseEl->nodeCloserToMin = variableNr++;
			
			// In quadratic order, there is one node per edge. In cubic order, there are two.
			if(orderType == 2){
				
				newSparseEl->nodeCloserToMax = newSparseEl->nodeCloserToMin;
				
				// We can already update the coordinates for that node :)
				theProblem->X[newSparseEl->nodeCloserToMin] = (theProblem->X[currentEdge.node[0]] + theProblem->X[currentEdge.node[1]])/2.0;
				theProblem->Y[newSparseEl->nodeCloserToMin] = (theProblem->Y[currentEdge.node[0]] + theProblem->Y[currentEdge.node[1]])/2.0;
				
			} else {
				
				newSparseEl->nodeCloserToMax = variableNr++;
				
				// We can already update the coordinates for those nodes :)
				theProblem->X[newSparseEl->nodeCloserToMin] = (2.0*theProblem->X[minNode] + theProblem->X[maxNode])/3.0;
				theProblem->Y[newSparseEl->nodeCloserToMin] = (2.0*theProblem->Y[minNode] + theProblem->Y[maxNode])/3.0;
				theProblem->X[newSparseEl->nodeCloserToMax] = (theProblem->X[minNode] + 2.0*theProblem->X[maxNode])/3.0;
				theProblem->Y[newSparseEl->nodeCloserToMax] = (theProblem->Y[minNode] + 2.0*theProblem->Y[maxNode])/3.0;
				
			}
			newSparseEl->next = NULL;
			
			// Adding the sparseEl at the end of the list associated to its min node
			if(sparseMatrix[minNode] == NULL){
				sparseMatrix[minNode] = newSparseEl;
			} else {
				current = sparseMatrix[minNode];
				while(current->next != NULL){
					current = current->next;
				}
				current->next = newSparseEl;
			}
		}
		
		// Note : we just created the structure but map still needs to be updated !
		
		/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		 * 2 - Support for creating node(s) on the elements of the mesh
		 * ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		 */
		
		// So we know when to start counting to number nodes inside an element
		int totalNodesOnEdges = 0; // shall contain one of the following : 0, 3, 4, 6 or 8
		int startNrForElements = 0;
		
		if (nLocalNode == 4) {
		    switch (orderType) {
		        case 2 : totalNodesOnEdges = 4; startNrForElements = nNode+nEdge; break;
		        case 3 : totalNodesOnEdges = 8; startNrForElements = nNode+2*nEdge; break;
			}
		} else if (nLocalNode == 3) {
		    switch (orderType) {
		        case 2 : totalNodesOnEdges = 3; startNrForElements = nNode+nEdge; break;
		        case 3 : totalNodesOnEdges = 6; startNrForElements = nNode+2*nEdge; break;
			}
		}
		
		/* +++++++++++++++++++++++++++++++++++++++++++++++++++++
		 * 3 - Assessing values to map + last updates on X and Y
		 * +++++++++++++++++++++++++++++++++++++++++++++++++++++
		 */
		
		int nrToSimplyIncrement = startNrForElements; // For unknows corresponding to nodes inside elements
		char doYouWantCloserToMin = 0; // To select which node on the edge we have to select (used only in the cubic case)
		
		
		for (i=0; i < nElem; i++) {
			
			
			// -------------- FOR DEBUGGING PURPOSES ONLY --------------------
			// -- (although there is no bug --> insert Brook's laugh here) --
			//printf("--------------------------------------------------");
			//printf("ELEMENT %d\n",i);
			//for(j = 0; j<nLocalNode ; j++){
			//	printf("value of map at %d : %d\n",j,theProblem->map[i*nLocal+j]);
			//}
			// ---------------------------------------------------------------
			
			int minNode = 0, maxNode = 0;
		    for (j=nLocalNode; j < nLocal; j++){
				
				// Concerning nodes on the edges
				if(j < nLocalNode + totalNodesOnEdges){
					
					// If quadratic, it's easy : select the node of the corresponding edge by sending a request to the sparse matrix
					if(orderType == 2){
						
						// Selecting the edge by getting its two endpoints
						minNode = MIN(theProblem->mesh->elem[(i-1)*nLocalNode+j], theProblem->mesh->elem[i*nLocalNode+(j-nLocalNode+1)%nLocalNode]);
						maxNode = MAX(theProblem->mesh->elem[(i-1)*nLocalNode+j], theProblem->mesh->elem[i*nLocalNode+(j-nLocalNode+1)%nLocalNode]);
						
						theProblem->map[i*nLocal+j] = searchForIt(sparseMatrix, minNode, maxNode, doYouWantCloserToMin);
						
						// If cubic, we need to pay attention at the direction we're browsing the edge before requesting a node !
					} else {
						
						// Little trick to determine what edge to consider given the variable j
						int originJ = j-nLocalNode;
						if((originJ%2) == 0){
							
							// Selecting the edge by getting its two endpoints
							int firstNode = originJ/2;
							minNode = MIN(theProblem->mesh->elem[i*nLocalNode+firstNode], theProblem->mesh->elem[i*nLocalNode+(firstNode+1)%nLocalNode]);
							maxNode = MAX(theProblem->mesh->elem[i*nLocalNode+firstNode], theProblem->mesh->elem[i*nLocalNode+(firstNode+1)%nLocalNode]);
							
							// This is where we decide which node to request
							if(minNode == theProblem->mesh->elem[i*nLocalNode+firstNode]){
								doYouWantCloserToMin = 1;
							} else {
								doYouWantCloserToMin = 0;
							}
							
						} else {
							doYouWantCloserToMin = 1-doYouWantCloserToMin;
						}
						
						theProblem->map[i*nLocal+j] = searchForIt(sparseMatrix, minNode, maxNode, doYouWantCloserToMin);
						
					}
					
					
					// Concerning nodes inside the element (uber easy)
				} else {
					
					theProblem->map[i*nLocal+j] = nrToSimplyIncrement++;
					femDiffusionProblem* tP = theProblem; // Simple copy
					
					// Assigning coordinates to that node
					if( (nLocalNode == 3 && orderType == 3) || (nLocalNode == 4 && orderType == 2) ){
						theProblem->X[theProblem->map[i*nLocal+j]] = 0;
						theProblem->Y[theProblem->map[i*nLocal+j]] = 0;
						int k = 0;
						for(k = 0; k < nLocalNode; k++){
							theProblem->X[theProblem->map[i*nLocal+j]] += xP(tP,i*nLocal+k);
							theProblem->Y[theProblem->map[i*nLocal+j]] += yP(tP,i*nLocal+k);
						}
						theProblem->X[theProblem->map[i*nLocal+j]] /= nLocalNode;
						theProblem->Y[theProblem->map[i*nLocal+j]] /= nLocalNode;
						
					} else if(nLocalNode == 4 && orderType == 3) {
						
						
						// Ugly code ? Yes ! Who cares ?
						int o = i*nLocal;
						
						float x4 = xP(tP,o+4), y4 = yP(tP,o+4), x9 = xP(tP,o+9), y9 = yP(tP,o+9);
						float a49 = (y4-y9)/(x4-x9), b49 = y9-a49*x9;
						
						float x5 = xP(tP,o+5), y5 = yP(tP,o+5), x8 = xP(tP,o+8), y8 = yP(tP,o+8);
						float a58 = (y5-y8)/(x5-x8), b58 = y8-a58*x8;
						
						float x6 = xP(tP,o+6), y6 = yP(tP,o+6), x11 = xP(tP,o+11), y11 = yP(tP,o+11);
						float a611 = (y6-y11)/(x6-x11), b611 = y11-a611*x11;
						
						float x7 = xP(tP,o+7), y7 = yP(tP,o+7), x10 = xP(tP,o+10), y10 = yP(tP,o+10);
						float a710 = (y7-y10)/(x7-x10), b710 = y10-a710*x10;
						
						int originJ = j - (nLocalNode + totalNodesOnEdges);
						if(originJ == 0){
							theProblem->X[theProblem->map[i*nLocal+j]] = (b611-b49)/(a49-a611);
							theProblem->Y[theProblem->map[i*nLocal+j]] = a49*theProblem->X[theProblem->map[i*nLocal+j]]+b49;
						} else if(originJ == 1){
							theProblem->X[theProblem->map[i*nLocal+j]] = (b611-b58)/(a58-a611);
							theProblem->Y[theProblem->map[i*nLocal+j]] = a58*theProblem->X[theProblem->map[i*nLocal+j]]+b58;
						} else if(originJ == 2){
							theProblem->X[theProblem->map[i*nLocal+j]] = (b710-b58)/(a58-a710);
							theProblem->Y[theProblem->map[i*nLocal+j]] = a58*theProblem->X[theProblem->map[i*nLocal+j]]+b58;
						} else if(originJ == 3){
							theProblem->X[theProblem->map[i*nLocal+j]] = (b710-b49)/(a49-a710);
							theProblem->Y[theProblem->map[i*nLocal+j]] = a49*theProblem->X[theProblem->map[i*nLocal+j]]+b49;
						} else {
							printf("ERROR !!! CAN'T SET POSITION OF CENTER ELEMENT FOR QUAD AT CUBIC !\n");
							theProblem->X[theProblem->map[i*nLocal+j]] = 0;
							theProblem->Y[theProblem->map[i*nLocal+j]] = 0;
						}
						
					}
				}
				
				
				// -------------- FOR DEBUGGING PURPOSES ONLY --------------------
				//printf("value of map at %d : %d\n",j,theProblem->map[i*nLocal+j]);
				// ---------------------------------------------------------------
				
				
			}
		}
		
		// Setting coordinates --> this is Legatus' todo-code, but I have already done the job in -1- and -3- :)
		//for (i=nNode; i < theProblem->size; i++) {
		//        theProblem->X[i] = 0.0;
		//        theProblem->Y[i] = 0.0;
		//}
		
    }
	
	//
	// A decommenter pour imprimer votre tableau d'appartenance des noeuds
	// Attention : bien respecter la numÃ©rotation de l'enonce
	// Il faut Ãªtre attentif Ã  la numÃ©rotation des aretes
	//
	// femDiffusionPrintMap(theProblem);
}


void femDiffusionPrintMap(femDiffusionProblem *theProblem)
{
    femMesh  *theMesh = theProblem->mesh;
    femEdges *theEdges = theProblem->edges;
	
    int i,j;
    int nElem = theMesh->nElem;
    int nNode = theMesh->nNode;
    int nEdge = theEdges->nEdge;
    int nLocal = theProblem->space->n;
    int order = theProblem->space->order;
    
    printf( "Number of Elements %d \n", nElem);
    for (i = 0; i < nElem; ++i) {
        int *elem = &(theProblem->map[i*nLocal]);
        printf("%6d : ", i);
        for (j=0; j < nLocal; j++)
            printf("%6d ", elem[j]);
        printf("\n");  }
	
    printf( "Number of segments %d \n", nEdge);
    for (i = 0; i < nEdge;i++) {
        printf("%6d : ", i);
        printf("%6d ", theEdges->edges[i].node[0]);
        for (j=1; j < order; j++)
            printf("%6d ", i*(order-1)+nNode+j-1);
        printf("%6d ", theEdges->edges[i].node[1]);
        printf("\n");  }
    fflush(stdout);
}


int searchForIt(sparseEl** sparseMatrix, int minNode, int maxNode, char doYouWantCloserToMin){
	
	sparseEl* element = sparseMatrix[minNode];
	if(element == NULL){
		printf("ERROR ! SPARSE MATRIX EMPTY FOR NODE %d !!\n",minNode);
		return 0;
	} else {
		
		if(element->nodeID == maxNode){
			if(doYouWantCloserToMin == 1){
				return element->nodeCloserToMin;
			}
			return element->nodeCloserToMax;
		}
		
		while(element->nodeID != maxNode && element->next != NULL){
			element = element->next;
			if(element->nodeID == maxNode){
				if(doYouWantCloserToMin == 1){
					return element->nodeCloserToMin;
				}
				return element->nodeCloserToMax;
			}
		}
		
	}
	
	printf("ERROR ! COULDN'T FIND THE PAIR (%d,%d) IN SPARSE MATRIX !!\n",minNode,maxNode);
	return 0;
	
}

double xP(femDiffusionProblem *theProblem, int mapID){
	return theProblem->X[theProblem->map[mapID]];
}

double yP(femDiffusionProblem *theProblem, int mapID){
	return theProblem->Y[theProblem->map[mapID]];
}

void freeSparseMatrix(sparseEl** sparseMatrix, int nNode){
	
}
