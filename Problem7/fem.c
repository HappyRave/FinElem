/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2013 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

static const double _gaussEdge2Xsi[2]     = { 0.577350269189626,-0.577350269189626 };
static const double _gaussEdge2Weight[2]  = { 1.000000000000000, 1.000000000000000 };

static const double _gaussQuad4Xsi[4]     = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Eta[4]     = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Weight[4]  = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000 };
static const double _gaussQuad9Xsi[9]     = { 0.774596669241483, 0.000000000000000,-0.774596669241483, 
                                              0.774596669241483, 0.000000000000000,-0.774596669241483, 
                                              0.774596669241483, 0.000000000000000,-0.774596669241483 };
static const double _gaussQuad9Eta[9]     = { 0.774596669241483, 0.774596669241483, 0.774596669241483, 
                                              0.000000000000000, 0.000000000000000, 0.000000000000000, 
                                             -0.774596669241483,-0.774596669241483,-0.774596669241483 };
static const double _gaussQuad9Weight[9]  = { 0.308641975308642, 0.493827160493827, 0.308641975308642, 
                                              0.493827160493827, 0.790123456790123, 0.493827160493827, 
                                              0.308641975308642, 0.493827160493827, 0.308641975308642 };

static const double _gaussTri3Xsi[3]      = { 0.166666666666667, 0.666666666666667, 0.166666666666667 };
static const double _gaussTri3Eta[3]      = { 0.166666666666667, 0.166666666666667, 0.666666666666667 };
static const double _gaussTri3Weight[3]   = { 0.166666666666667, 0.166666666666667, 0.166666666666667 };
static const double _gaussTri12Xsi[12]    = { 0.249286745170910, 0.249286745170910, 0.501426509658179, 
                                              0.063089014491502, 0.063089014491502, 0.873821971016996, 
                                              0.310352451033785, 0.636502499121399, 0.053145049844816, 
                                              0.310352451033785, 0.636502499121399, 0.053145049844816 };
static const double _gaussTri12Eta[12]    = { 0.249286745170910, 0.501426509658179, 0.249286745170910,
                                              0.063089014491502, 0.873821971016996, 0.063089014491502,
                                              0.636502499121399, 0.053145049844816, 0.310352451033785,
                                              0.053145049844816, 0.310352451033785, 0.636502499121399 };
static const double _gaussTri12Weight[12] = { 0.058393137863189, 0.058393137863189, 0.058393137863189,
                                              0.025422453185104, 0.025422453185104, 0.025422453185104,
                                              0.041425537809187, 0.041425537809187, 0.041425537809187,
                                              0.041425537809187, 0.041425537809187, 0.041425537809187 };

                                              
                                            
                                                        


femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight; }
    else if (type == FEM_QUAD && n == 9) {
        theRule->n      = 9;
        theRule->xsi    = _gaussQuad9Xsi;
        theRule->eta    = _gaussQuad9Eta;
        theRule->weight = _gaussQuad9Weight; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else if (type == FEM_TRIANGLE && n == 12) {
        theRule->n      = 12;
        theRule->xsi    = _gaussTri12Xsi;
        theRule->eta    = _gaussTri12Eta;
        theRule->weight = _gaussTri12Weight; }
    else if (type == FEM_EDGE && n == 2) {
        theRule->n      = 2;
        theRule->xsi    = _gaussEdge2Xsi;
        theRule->eta    = NULL;
        theRule->weight = _gaussEdge2Weight; }

    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}


void _1c0_x(double *xsi) 
{
    xsi[0] = -1.0; 
    xsi[1] =  1.0; 
}


void _1c0_phi(double xsi, double *phi) 
{   
    phi[0] = (1.0 - xsi)/2.0;
    phi[1] = (1.0 + xsi)/2.0;   
    
}

void _1c0_dphidx(double xsi, double *dphidxsi) 
{   
    dphidxsi[0] =  -1.0/2.0;
    dphidxsi[1] =   1.0/2.0;  
}

void _3c0_x(double *xsi) 
{
    xsi[0] = -1.0; 
    xsi[1] = -1.0/3.0; 
    xsi[1] =  1.0/3.0; 
    xsi[1] =  1.0; 
}

void _3c0_phi(double xsi, double *phi)
{
  	phi[0] =   9./16 * (-1./3 - xsi) * ( 1./3 - xsi) * (1.   - xsi);
  	phi[1] = -27./16 * (-1.   - xsi) * ( 1./3 - xsi) * (1.   - xsi);
  	phi[2] =  27./16 * (-1.   - xsi) * (-1./3 - xsi) * (1.   - xsi);
  	phi[3] =  -9./16 * (-1.   - xsi) * (-1./3 - xsi) * (1./3 - xsi);
}

void _3c0_dphidx(double xsi, double *dphidxsi)
{
  	dphidxsi[0] =   9./16 * ( - ( 1./3 - xsi) * (1.   - xsi) - (-1./3 - xsi) * (1.   - xsi) - (-1./3 - xsi) * ( 1./3 - xsi) ) ;
  	dphidxsi[1] = -27./16 * ( - ( 1./3 - xsi) * (1.   - xsi) - (-1.   - xsi) * (1.   - xsi) - (-1.   - xsi) * ( 1./3 - xsi) ) ;
  	dphidxsi[2] =  27./16 * ( - (-1./3 - xsi) * (1.   - xsi) - (-1.   - xsi) * (1.   - xsi) - (-1.   - xsi) * (-1./3 - xsi) ) ;
  	dphidxsi[3] =  -9./16 * ( - (-1./3 - xsi) * (1./3 - xsi) - (-1.   - xsi) * (1./3 - xsi) - (-1.   - xsi) * (-1./3 - xsi) ) ;
}



void _q1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;  
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =   (1.0 + eta) / 4.0;  
    dphidxsi[1] = - (1.0 + eta) / 4.0;
    dphidxsi[2] = - (1.0 - eta) / 4.0;
    dphidxsi[3] =   (1.0 - eta) / 4.0;
    dphideta[0] =   (1.0 + xsi) / 4.0;  
    dphideta[1] =   (1.0 - xsi) / 4.0;
    dphideta[2] = - (1.0 - xsi) / 4.0;
    dphideta[3] = - (1.0 + xsi) / 4.0;

}

void _q2c0_x(double *xsi, double *eta) 
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
    xsi[4] =  0.0;  eta[4] =  1.0;
    xsi[5] = -1.0;  eta[5] =  0.0;
    xsi[6] =  0.0;  eta[6] = -1.0;
    xsi[7] =  1.0;  eta[7] =  0.0;
    xsi[8] =  0.0;  eta[8] =  0.0;
}

void _q2c0_phi(double xsi, double eta, double *phi)
{
    phi[0] =  xsi*(1.0+xsi)*eta*(1.0+eta)/4.0;
    phi[1] = -xsi*(1.0-xsi)*eta*(1.0+eta)/4.0;
    phi[2] =  xsi*(1.0-xsi)*eta*(1.0-eta)/4.0;
    phi[3] = -xsi*(1.0+xsi)*eta*(1.0-eta)/4.0;
    phi[4] =  (1.0-xsi*xsi)*eta*(1.0+eta)/2.0;
    phi[5] = -xsi*(1.0-xsi)*(1.0-eta*eta)/2.0;
    phi[6] = -(1.0-xsi*xsi)*eta*(1.0-eta)/2.0;
    phi[7] =  xsi*(1.0+xsi)*(1.0-eta*eta)/2.0;
    phi[8] =  (1.0-xsi*xsi)*(1.0-eta*eta);
}

void _q2c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =  (1.0+2.0*xsi)*eta*(1.0+eta)/4.0;
    dphidxsi[1] = (-1.0+2.0*xsi)*eta*(1.0+eta)/4.0;
    dphidxsi[2] =  (1.0-2.0*xsi)*eta*(1.0-eta)/4.0;
    dphidxsi[3] = -(1.0+2.0*xsi)*eta*(1.0-eta)/4.0;
    dphidxsi[4] =       -2.0*xsi*eta*(1.0+eta)/2.0;
    dphidxsi[5] = (-1.0+2.0*xsi)*(1.0-eta*eta)/2.0;
    dphidxsi[6] =        2.0*xsi*eta*(1.0-eta)/2.0;
    dphidxsi[7] =  (1.0+2.0*xsi)*(1.0-eta*eta)/2.0;
    dphidxsi[8] =       -2.0*xsi*(1.0-eta*eta);
    dphideta[0] =  xsi*(1.0+xsi)*(1.0+2.0*eta)/4.0;
    dphideta[1] = -xsi*(1.0-xsi)*(1.0+2.0*eta)/4.0;
    dphideta[2] =  xsi*(1.0-xsi)*(1.0-2.0*eta)/4.0;
    dphideta[3] = -xsi*(1.0+xsi)*(1.0-2.0*eta)/4.0;
    dphideta[4] =  (1.0-xsi*xsi)*(1.0+2.0*eta)/2.0;
    dphideta[5] =  xsi*(1.0-xsi)*2.0*eta/2.0;
    dphideta[6] = -(1.0-xsi*xsi)*(1.0-2.0*eta)/2.0;
    dphideta[7] =  xsi*(1.0+xsi)*(-2.0*eta)/2.0;
    dphideta[8] =  (1.0-xsi*xsi)*(-2.0*eta);
}

void _q3c0_x(double *xsi, double *eta) 
{
    xsi[0] =  1.0;       eta[0] =  1.0;
    xsi[1] = -1.0;       eta[1] =  1.0;
    xsi[2] = -1.0;       eta[2] = -1.0;
    xsi[3] =  1.0;       eta[3] = -1.0;
    xsi[4] =  1.0/3.0;   eta[4] =  1.0;
    xsi[5] = -1.0/3.0;   eta[5] =  1.0;
    xsi[6] = -1.0;       eta[6] =  1.0/3.0;
    xsi[7] = -1.0;       eta[7] = -1.0/3.0;
    xsi[8] = -1.0/3.0;   eta[8] = -1.0;
    xsi[9] =  1.0/3.0;   eta[9] = -1.0;
    xsi[10] =  1.0;      eta[10] = -1.0/3.0;
    xsi[11] =  1.0;      eta[11] =  1.0/3.0;
    xsi[12] =  1.0/3.0;  eta[12] =  1.0/3.0;
    xsi[13] = -1.0/3.0;  eta[13] =  1.0/3.0;
    xsi[14] = -1.0/3.0;  eta[14] = -1.0/3.0;
    xsi[15] =  1.0/3.0;  eta[15] = -1.0/3.0;
}


void _q3c0_phi(double xsi, double eta, double *phi) 
{
  	double fXi[4], fEta[4];
  	int map[16] = {2,8,9,3,7,14,15,10,6,13,12,11,1,5,4,0};
  	_3c0_phi(xsi, fXi);
  	_3c0_phi(eta, fEta);
  	int i, j;
  	for (i = 0; i < 4; ++i) 
    	for (j = 0; j < 4; ++j) 
      		phi[map[j*4+i]] = fXi[i] * fEta[j];
}

void _q3c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta) 
{

  	double fXi[4], fEta[4], dfXi[4], dfEta[4];
    int map[16] = {2,8,9,3,7,14,15,10,6,13,12,11,1,5,4,0};

  	_3c0_phi(xsi, fXi);
  	_3c0_phi(eta, fEta);
  	_3c0_dphidx(xsi, dfXi);
  	_3c0_dphidx(eta, dfEta);
  	int i, j;
  	for (i = 0; i < 4; ++i) {
    	for (j = 0; j < 4; ++j) {
      		dphidxsi[map[j*4+i]] = dfXi[i] * fEta[j];
      		dphideta[map[j*4+i]] = fXi[i] * dfEta[j]; }}
}


void _p1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;  
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;  
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;

}

void _p2c0_x(double *xsi, double *eta) 
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
    xsi[3] =  0.5;  eta[3] =  0.0;
    xsi[4] =  0.5;  eta[4] =  0.5;
    xsi[5] =  0.0;  eta[5] =  0.5;
}


void _p2c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1.0 - 3.0*(xsi+eta) + 2.0*(xsi+eta)*(xsi+eta);  
    phi[1] = xsi*(2.0*xsi-1.0);
    phi[2] = eta*(2.0*eta-1.0);
    phi[3] = 4.0*xsi*(1.0-xsi-eta);
    phi[4] = 4.0*xsi*eta;
    phi[5] = 4.0*eta*(1.0-xsi-eta);
}

void _p2c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = - 3.0 + 4.0*xsi + 4.0*eta;   
    dphidxsi[1] = - 1.0 + 4.0*xsi          ;
    dphidxsi[2] =   0.0                    ;
    dphidxsi[3] =   4.0 - 8.0*xsi - 4.0*eta;  
    dphidxsi[4] =                   4.0*eta;
    dphidxsi[5] =                 - 4.0*eta;
    dphideta[0] = - 3.0 + 4.0*xsi + 4.0*eta; 
    dphideta[1] =   0.0                    ;
    dphideta[2] =  -1.0           + 4.0*eta;
    dphideta[3] =       - 4.0*xsi          ;
    dphideta[4] =         4.0*xsi          ;
    dphideta[5] =   4.0 - 4.0*xsi - 8.0*eta;  
}

void _p3c0_x(double *xsi, double *eta) 
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

void _p3c0_phi(double xsi, double eta, double *phi) 
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

void _p3c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta) 
{   
    dphidxsi[0] = (- (1.0/3.0 - xsi - eta) * (2.0/3.0 - xsi - eta) - (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) - (1.0 - xsi - eta) * (1.0/3.0 - xsi - eta)) *  9.0/2.0; 
    dphidxsi[1] = (  (xsi - 1.0/3.0) * (xsi - 2.0/3.0) +  xsi * (xsi - 2.0/3.0) + xsi * (xsi - 1.0/3.0) ) *  9.0/2.0;  
    dphidxsi[2] = 0.0;  
    dphidxsi[3] = ( - (2.0/3.0 - xsi - eta) * xsi - (1.0 - xsi - eta) * xsi + (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) ) * 27.0/2.0;
    dphidxsi[4] = ( - (xsi - 1.0/3.0) * xsi + (1.0 - xsi - eta) * xsi + (1.0 - xsi - eta) * (xsi - 1.0/3.0) ) * 27.0/2.0;
    dphidxsi[5] = ( eta * (xsi - 1.0/3.0) + xsi * eta ) * 27.0/2.0;
    dphidxsi[6] =   eta * (eta - 1.0/3.0) * 27.0/2.0;  
    dphidxsi[7] = - (eta - 1.0/3.0) * eta * 27.0/2.0;
    dphidxsi[8] = ( - (2.0/3.0 - xsi - eta) * eta - (1.0 - xsi - eta) * eta  ) * 27.0/2.0;
    dphidxsi[9] = ( - xsi * eta + (1.0 - xsi - eta) * eta) * 27.0; 
    dphideta[0] = (- (1.0/3.0 - xsi - eta) * (2.0/3.0 - xsi - eta) - (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta) - (1.0 - xsi - eta) * (1.0/3.0 - xsi - eta)) *  9.0/2.0; 
    dphideta[1] = 0.0;
    dphideta[2] =  (  (eta - 1.0/3.0) * (eta - 2.0/3.0) +  eta * (eta - 2.0/3.0) + eta * (eta - 1.0/3.0) ) *  9.0/2.0; 
    dphideta[3] = ( - (2.0/3.0 - xsi - eta) * xsi - (1.0 - xsi - eta) * xsi ) * 27.0/2.0;
    dphideta[4] = - (xsi - 1.0/3.0)  * xsi * 27.0/2.0;
    dphideta[5] =   xsi * (xsi - 1.0/3.0) * 27.0/2.0;;
    dphideta[6] = ( xsi * (eta - 1.0/3.0) +  xsi * eta ) * 27.0/2.0;  
    dphideta[7] = (- (eta - 1.0/3.0) * eta + (1.0 - xsi - eta) * eta + (1.0 - xsi - eta) * (eta - 1.0/3.0) ) * 27.0/2.0;
    dphideta[8] = ( - (2.0/3.0 - xsi - eta) * eta - (1.0 - xsi - eta) * eta + (1.0 - xsi - eta) * (2.0/3.0 - xsi - eta)) * 27.0/2.0;
    dphideta[9] = ( - xsi * eta + (1.0 - xsi - eta) * xsi) * 27.0;              
             
}



femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_QUAD && n == 4) {
        theSpace->n       = 4;
        theSpace->order   = 1;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx; 
        theSpace->x1      = _1c0_x;
        theSpace->phi1    = _1c0_phi;
        theSpace->dphi1dx = _1c0_dphidx;}
    else if (type == FEM_QUAD && n == 9) {
        theSpace->n       = 9;
        theSpace->order   = 2;
        theSpace->x2      = _q2c0_x;
        theSpace->phi2    = _q2c0_phi;
        theSpace->dphi2dx = _q2c0_dphidx; }
    else if (type == FEM_QUAD && n == 16) {
        theSpace->n       = 16;
        theSpace->order   = 3;
        theSpace->x2      = _q3c0_x;
        theSpace->phi2    = _q3c0_phi;
        theSpace->dphi2dx = _q3c0_dphidx; 
        theSpace->x1      = _3c0_x;
        theSpace->phi1    = _3c0_phi;
        theSpace->dphi1dx = _3c0_dphidx; } 
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->order   = 1;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx;
        theSpace->x1      = _1c0_x;
        theSpace->phi1    = _1c0_phi;
        theSpace->dphi1dx = _1c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 6) {
        theSpace->n       = 6;
        theSpace->order   = 2;
        theSpace->x2      = _p2c0_x;
        theSpace->phi2    = _p2c0_phi;
        theSpace->dphi2dx = _p2c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 10) {
        theSpace->n       = 10;
        theSpace->order   = 3;
        theSpace->x2      = _p3c0_x;
        theSpace->phi2    = _p3c0_phi;
        theSpace->dphi2dx = _p3c0_dphidx; 
        theSpace->x1      = _3c0_x;
        theSpace->phi1    = _3c0_phi;
        theSpace->dphi1dx = _3c0_dphidx;}
    else Error("Cannot create such a discrete space !");
    return theSpace; 
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscreteXsi1(femDiscrete* mySpace, double *xsi)
{
    mySpace->x1(xsi);
}

void femDiscretePhi1(femDiscrete* mySpace, double xsi, double *phi)
{
    mySpace->phi1(xsi,phi);
}

void femDiscreteDphi1(femDiscrete* mySpace, double xsi, double *dphidxsi)
{
    mySpace->dphi1dx(xsi,dphidxsi);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[16], eta[16], phi[16], dphidxsi[16], dphideta[16];
    
    femDiscreteXsi2(mySpace,xsi,eta);
    for (i=0; i < n; i++) {
        
        femDiscretePhi2(mySpace,xsi[i],eta[i],phi);
        femDiscreteDphi2(mySpace,xsi[i],eta[i],dphidxsi,dphideta);

        for (j=0; j < n; j++)  {
            printf("(xsi=%+.1f,eta=%+.1f) : ",xsi[i],eta[i]);
            printf(" phi(%d)=%+.1f",j,phi[j]);  
            printf("   dphidxsi(%d)=%+.1f",j,dphidxsi[j]);  
            printf("   dphideta(%d)=%+.1f \n",j,dphideta[j]);  }
        printf(" \n"); }
}



femMesh *femMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash,*elem;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    fscanf(file, "Number of nodes %d \n", &theMesh->nNode);
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fscanf(file,"%d : %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i]); }
    
    char str[256]; fgets(str, sizeof(str), file);
    if (!strncmp(str,"Number of triangles",19))  { 
        sscanf(str,"Number of triangles %d \n", &theMesh->nElem);
        theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
        theMesh->nLocalNode = 3;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            fscanf(file,"%d : %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2]); }}
    else if (!strncmp(str,"Number of quads",15))  { 
        sscanf(str,"Number of quads %d \n", &theMesh->nElem);  
        theMesh->elem = malloc(sizeof(int)*4*theMesh->nElem);
        theMesh->nLocalNode = 4;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
        fscanf(file,"%d : %d %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2],&elem[3]); }}
  
    fclose(file);
     return theMesh;
}

void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->elem);
    free(theMesh);
}

void femMeshWrite(const femMesh *theMesh, const char *filename)
{
    int i,*elem;
    
    FILE* file = fopen(filename,"w");
    
    fprintf(file, "Number of nodes %d \n", theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fprintf(file,"%6d : %14.7e %14.7e \n",i,theMesh->X[i],theMesh->Y[i]); }
    
    if (theMesh->nLocalNode == 4) {
        fprintf(file, "Number of quads %d \n", theMesh->nElem);  
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
            fprintf(file,"%6d : %6d %6d %6d %6d \n", i,elem[0],elem[1],elem[2],elem[3]);   }}
    else if (theMesh->nLocalNode == 3) {
        fprintf(file, "Number of triangles %d \n", theMesh->nElem);  
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            fprintf(file,"%6d : %6d %6d %6d \n", i,elem[0],elem[1],elem[2]);   }}
    
    fclose(file);
}
                                                                                                          
femEdges *femEdgesCreate(femMesh *theMesh)
{
    femEdges *theEdges = malloc(sizeof(femEdges));
    int nLoc = theMesh->nLocalNode;
    int i,j,n = theMesh->nElem * nLoc;
    femEdge* edges = malloc(n * sizeof(femEdge));
    theEdges->mesh  = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;
    
    for (i = 0; i < theMesh->nElem; i++) {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++) {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc]; }}

    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), femEdgesCompare);

    int index = 0;          
    int nBoundary = 0;
    
    for (i=0; i < theEdges->nEdge; i++) {
      if (i == theEdges->nEdge - 1 || femEdgesCompare(&edges[i],&edges[i+1]) != 0) {
              edges[index] = edges[i];
              nBoundary++; }
      else {  edges[index] = edges[i];
              edges[index].elem[1] = edges[i+1].elem[0];
              i = i+1;}
      index++; }
      
    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
    return theEdges;
}

void femEdgesPrint(femEdges *theEdges)
{
    int i;    
    for (i = 0; i < theEdges->nEdge; ++i) {
        printf("%6d : %4d %4d : %4d %4d \n",i,
               theEdges->edges[i].node[0],theEdges->edges[i].node[1],
               theEdges->edges[i].elem[0],theEdges->edges[i].elem[1]); }
}

void femEdgesFree(femEdges *theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

int femEdgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*) edgeOne)->node;
    int *nodeTwo = ((femEdge*) edgeTwo)->node;  
    int  diffMin = fmin(nodeOne[0],nodeOne[1]) - fmin(nodeTwo[0],nodeTwo[1]);
    int  diffMax = fmax(nodeOne[0],nodeOne[1]) - fmax(nodeTwo[0],nodeTwo[1]);
    
    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1; 
                        return  0;
}

femSolver *femSolverFullCreate(int size)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_FULL;
    mySolver->solver = (femSolver *)femFullSystemCreate(size);
    return(mySolver);
}

femSolver *femSolverBandCreate(int size, int band)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_BAND;
    mySolver->solver = (femSolver *)femBandSystemCreate(size,band);
    return(mySolver);
}

femSolver *femSolverIterativeCreate(int size)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_ITER;
    mySolver->solver = (femSolver *)femIterativeSolverCreate(size);
    return(mySolver);
}

void femSolverFree(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemFree((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemFree((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverFree((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    free(mySolver);
}

void femSolverInit(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemInit((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemInit((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverInit((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}


double femSolverGet(femSolver *mySolver,int i,int j)
{
    double value = 0;
    switch (mySolver->type) {
        case FEM_FULL : value = femFullSystemGet((femFullSystem *)mySolver->solver,i,j); break;
        case FEM_BAND : value = femBandSystemGet((femBandSystem *)mySolver->solver,i,j); break;
        case FEM_ITER : value = (i==j); break;
        default : Error("Unexpected solver type"); }
    return(value);
}

void femSolverPrint(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrint((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrint((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrint((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

void femSolverPrintInfos(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrintInfos((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrintInfos((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrintInfos((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

  
void femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemAssemble((femFullSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_BAND : femBandSystemAssemble((femBandSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_ITER : femIterativeSolverAssemble((femIterativeSolver *)mySolver->solver,Aloc,Bloc,Uloc,map,nLoc); break;
        default : Error("Unexpected solver type"); }
}
    
void femSolverConstrain(femSolver *mySolver, int myNode, double myValue)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemConstrain((femFullSystem *)mySolver->solver,myNode,myValue); break;
        case FEM_BAND : femBandSystemConstrain((femBandSystem *)mySolver->solver,myNode,myValue); break;
        case FEM_ITER : femIterativeSolverConstrain((femIterativeSolver *)mySolver->solver,myNode,myValue); break;
        default : Error("Unexpected solver type"); }
}
  
double *femSolverEliminate(femSolver *mySolver)
{
    double *soluce;
    switch (mySolver->type) {
        case FEM_FULL : soluce = femFullSystemEliminate((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : soluce = femBandSystemEliminate((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : soluce = femIterativeSolverEliminate((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(soluce);
}



int femSolverConverged(femSolver *mySolver)
{
    int  testConvergence;
    switch (mySolver->type) {
        case FEM_FULL : testConvergence = 1; break;
        case FEM_BAND : testConvergence = 1; break;
        case FEM_ITER : testConvergence = femIterativeSolverConverged((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(testConvergence);
}


femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    theSystem->A = malloc(sizeof(double*) * size); 
    theSystem->B = malloc(sizeof(double) * size * (size+1));
    theSystem->A[0] = theSystem->B + size;  
    theSystem->size = size;
    int i;
    for (i=1 ; i < size ; i++) 
        theSystem->A[i] = theSystem->A[i-1] + size;
    femFullSystemInit(theSystem);

    return theSystem; 
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}


void femFullSystemInit(femFullSystem *mySystem)
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++) 
        mySystem->B[i] = 0;}


void femFullSystemPrint(femFullSystem *mySystem)
{
    double  **A, *B;
    int     i, j, size;
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        for (j=0; j < size; j++)
            if (A[i][j] == 0)  printf("         ");   
            else               printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}

void femFullSystemPrintInfos(femFullSystem *mySystem)
{
    int  size = mySystem->size;
    printf(" \n");
    printf("    Full Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(size+1));     
}

void  femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }
    
    for (i=0; i < size; i++) 
        A[myNode][i] = 0; 
    
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}

void femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        for(j = 0; j < nLoc; j++) {
            mySystem->A[map[i]][map[j]] += Aloc[i*nLoc+j]; }
    mySystem->B[map[i]] += Bloc[i]; }
}

double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-8 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    
    return(mySystem->B);    
}

femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);        
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    int i;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);
}
 
void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A); 
    free(myBandSystem);
}
 
void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++) 
        myBandSystem->B[i] = 0;        
}
 
void femBandSystemPrint(femBandSystem *myBand)
{
    double  **A, *B;
    int     i, j, band, size;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    for (i=0; i < size; i++) {
        for (j=i; j < i+band; j++)
            if (A[i][j] == 0) printf("         ");   
            else              printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}
  
void femBandSystemPrintInfos(femBandSystem *myBand)
{
    int size = myBand->size;
    int band = myBand->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(band+1));     
}


double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol)
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol]; 
    return(value);
}

double femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol)
{
    return(myFullSystem->A[myRow][myCol]); 
}

void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        int myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            if (myCol >= myRow)  myBandSystem->A[myRow][myCol] += Aloc[i*nLoc+j]; }
        myBandSystem->B[myRow] += Bloc[i]; }
}


 
void femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue) 
{
    double  **A, *B;
    int     i, size, band, ifirst, iend;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    ifirst = fmax(0,myNode - band + 1);
    iend   = myNode;
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }

    ifirst = myNode+1;
    iend = fmin(myNode + band,size);
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[myNode][i];
        A[myNode][i] = 0; }
        
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}
 
double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    /* Incomplete Cholesky factorization */ 

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    /* Back-substitution */

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(myBand->B);
}



femIterativeSolver *femIterativeSolverCreate(int size)
{
    femIterativeSolver *mySolver = malloc(sizeof(femIterativeSolver));
    mySolver->R = malloc(sizeof(double)*size*4);  
    mySolver->R0 = mySolver->R + size;      
    mySolver->D  = mySolver->R + size*2;       
    mySolver->D0 = mySolver->R + size*3;       
    mySolver->size = size;
    femIterativeSolverInit(mySolver);
    return(mySolver);
}

void femIterativeSolverFree(femIterativeSolver *mySolver)
{
    free(mySolver->R);
    free(mySolver);
}

void femIterativeSolverInit(femIterativeSolver *mySolver)
{
    int i;
    mySolver->iter = 0;
    mySolver->error = 10.0e+12;
    for (i=0 ; i < mySolver->size*4 ; i++) 
        mySolver->R[i] = 0;        
}
 
void femIterativeSolverPrint(femIterativeSolver *mySolver)
{
    double  *R;
    int     i, size;
    R    = mySolver->R;
    size = mySolver->size;

    for (i=0; i < size; i++) {
        printf("%d :  %+.1e \n",i,R[i]); }
}

void femIterativeSolverPrintInfos(femIterativeSolver *mySolver)
{
    if (mySolver->iter == 1)     printf("\n    Iterative solver \n");
    printf("    Iteration %4d : %14.7e\n",mySolver->iter,mySolver->error);
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    int i,j,myRow;
    if (mySolver->iter==0) {
        for (i = 0; i < nLoc; i++) { 
            myRow = map[i];
            mySolver->R[myRow] += Bloc[i];
            mySolver->D[myRow] += Bloc[i];
            for(j = 0; j < nLoc; j++) {
                mySolver->R[myRow] -= Aloc[i*nLoc+j]*Uloc[j];
                mySolver->D[myRow] -= Aloc[i*nLoc+j]*Uloc[j]; }}}
    else {
    for (i = 0; i < nLoc; i++) {
        myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            mySolver->D0[myRow] += Aloc[i*nLoc+j] * mySolver->D[myCol]; }}}
}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue) 
{
    mySolver->R[myNode] = myValue;
    mySolver->D0[myNode] = myValue;
    mySolver->D[myNode] = myValue;
}

// R
// R0  = dX
// D   = P
// D0  = AP

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;
    double error = 0.0; int i;
    double denAlpha = 0.0;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i])*(mySolver->R[i]);
        denAlpha += mySolver->D[i] * mySolver->D0[i]; }
    double alpha = error/denAlpha;
    
    if (mySolver->iter == 1) {
        for (i=0; i < mySolver->size; i++) {
            mySolver->R0[i] = 0.0; }}
    else {
        double numBeta = 0.0;
        for (i=0; i < mySolver->size; i++) {
            mySolver->R0[i] = alpha * mySolver->D[i];
            mySolver->R[i] = mySolver->R[i] - alpha * mySolver->D0[i];
            numBeta += mySolver->R[i] * mySolver->R[i]; }
        double beta = numBeta/error;
        for (i=0; i < mySolver->size; i++) {
            mySolver->D0[i] = 0.0;
            mySolver->D[i] = mySolver->R[i] + beta * mySolver->D[i]; }}
   
    mySolver->error = sqrt(error);
    return(mySolver->R0);
}



int femIterativeSolverConverged(femIterativeSolver *mySolver)
{
    int  testConvergence = 0;
    if (mySolver->iter  > 3000)     testConvergence = -1;
    if (mySolver->error < 10.0e-6)  testConvergence = 1;
    return(testConvergence);
}

femDiffusionProblem *femDiffusionCreate(const char *filename, femSolverType solverType, femRenumType renumType, int orderType)
{
    int i,band;
    int nodesOrderTriangle[] = {1,3,6,10};
    int integOrderTriangle[] = {1,3,12,12};
    int nodesOrderQuad[]     = {1,4,9,16};
    int integOrderQuad[]     = {1,4,9,9};
    
    femDiffusionProblem *theProblem = malloc(sizeof(femDiffusionProblem));
    theProblem->mesh  = femMeshRead(filename);           
    theProblem->edges = femEdgesCreate(theProblem->mesh); 
    int nNode = theProblem->mesh->nNode; 
    int nElem = theProblem->mesh->nElem; 
    int nEdge = theProblem->edges->nEdge; 
    if (theProblem->mesh->nLocalNode == 4) {
        switch (orderType) {
            case 0 : theProblem->size = nElem; break;
            case 1 : theProblem->size = nNode; break;
            case 2 : theProblem->size = nNode + nEdge + nElem; break;
            case 3 : theProblem->size = nNode + 2*nEdge + 4*nElem; break;  }
        theProblem->space = femDiscreteCreate(nodesOrderQuad[orderType],FEM_QUAD);
        theProblem->rule = femIntegrationCreate(integOrderQuad[orderType],FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        switch (orderType) {
            case 0 : theProblem->size = nElem; break;
            case 1 : theProblem->size = nNode; break;
            case 2 : theProblem->size = nNode + nEdge; break;
            case 3 : theProblem->size = nNode + 2*nEdge + nElem; break;  }
        theProblem->space = femDiscreteCreate(nodesOrderTriangle[orderType],FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(integOrderTriangle[orderType],FEM_TRIANGLE); }
 
    theProblem->number  = malloc(sizeof(int)*theProblem->size);
    theProblem->map  = malloc((theProblem->space->n)*sizeof(int)*theProblem->mesh->nElem);
    theProblem->soluce = malloc(sizeof(double)*theProblem->size);
    theProblem->X = malloc(sizeof(double)*theProblem->size);
    theProblem->Y = malloc(sizeof(double)*theProblem->size);
    for (i = 0; i < theProblem->size; i++)      
        theProblem->soluce[i] = 0.0;  
    femDiffusionComputeMap(theProblem,orderType);

    
    femDiffusionRenumber(theProblem,renumType);
         
    switch (solverType) {
        case FEM_FULL : 
                theProblem->solver = femSolverFullCreate(theProblem->size); break;
        case FEM_BAND : 
                band = femDiffusionComputeBand(theProblem);
                theProblem->solver = femSolverBandCreate(theProblem->size,band); break;
        case FEM_ITER : 
               theProblem->solver = femSolverIterativeCreate(theProblem->size); break;
        default : Error("Unexpected solver option"); }
        
    
    return theProblem;
}

void femDiffusionFree(femDiffusionProblem *theProblem)
{
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    femSolverFree(theProblem->solver);
    free(theProblem->number);
    free(theProblem->soluce);
    free(theProblem->map);
    free(theProblem->X);
    free(theProblem->Y);
    free(theProblem);
}
    

void femDiffusionMeshLocal(const femDiffusionProblem *theProblem, const int iElem, int *map, double *x, double *y, double *u)
{
    int j,nLocal = theProblem->space->n;
    
    for (j=0; j < nLocal; ++j) {
        map[j] = theProblem->map[iElem*nLocal+j];
        x[j]   = theProblem->X[map[j]];
        y[j]   = theProblem->Y[map[j]]; 
        u[j]   = theProblem->soluce[map[j]];
        map[j] = theProblem->number[map[j]];
        }   
}

void femDiffusionCompute(femDiffusionProblem *theProblem)
{
    femMesh *theMesh = theProblem->mesh;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
    femSolver *theSolver = theProblem->solver;
    femEdges *theEdges = theProblem->edges;
    int *number = theProblem->number;

    int n = theSpace->n;
     
    double Xloc[n],Yloc[n],phi[n],dphidxsi[n],dphideta[n],dphidx[n],dphidy[n],Aloc[n*n],Bloc[n],Uloc[n];
    int iEdge,iElem,iInteg,i,j,map[n];
   
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (i = 0; i < theSpace->n; i++)      Bloc[i] = 0;
        for (i = 0; i < (theSpace->n)*(theSpace->n); i++) Aloc[i] = 0;
        femDiffusionMeshLocal(theProblem,iElem,map,Xloc,Yloc,Uloc);  
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0; 
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {    
                dxdxsi += Xloc[i]*dphidxsi[i];       
                dxdeta += Xloc[i]*dphideta[i];   
                dydxsi += Yloc[i]*dphidxsi[i];   
                dydeta += Yloc[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    Aloc[i*(theSpace->n)+j] += (dphidx[i] * dphidx[j] 
                                            + dphidy[i] * dphidy[j]) * jac * weight; }}                                                                                            
            for (i = 0; i < theSpace->n; i++) {
                Bloc[i] += phi[i] * jac *weight; }}
        femSolverAssemble(theSolver,Aloc,Bloc,Uloc,map,theSpace->n); } 

    for (iEdge= 0; iEdge < theEdges->nEdge; iEdge++) {      
        if (theEdges->edges[iEdge].elem[1] < 0) {       
            femSolverConstrain(theSolver,number[theEdges->edges[iEdge].node[0]],0.0); 
            femSolverConstrain(theSolver,number[theEdges->edges[iEdge].node[1]],0.0); 
            for (j=1; j < theSpace->order; j++)
                femSolverConstrain(theSolver,number[iEdge*(theSpace->order-1)+theMesh->nNode+j-1],0.0);
            }}
  
    double *soluce = femSolverEliminate(theSolver);
    for (i = 0; i < theProblem->size; i++)
        theProblem->soluce[i] += soluce[number[i]];
  
}

static double *theGlobal;

int compare(const void *nodeOne, const void *nodeTwo) 
{
    int *iOne = (int *)nodeOne;
    int *iTwo = (int *)nodeTwo;
    double diff = theGlobal[*iOne] - theGlobal[*iTwo];
    if (diff < 0)    return  1;
    if (diff > 0)    return -1;
                     return  0;  
}

void femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType)
{
    int i, *inverse;
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theProblem->size; i++) 
                theProblem->number[i] = i;
            break;
        case FEM_XNUM : 
        case FEM_YNUM :
            inverse = malloc(sizeof(int)*theProblem->size);
            for (i = 0; i < theProblem->size; i++) 
                inverse[i] = i; 
            if (renumType == FEM_XNUM) theGlobal = theProblem->X;
            if (renumType == FEM_YNUM) theGlobal = theProblem->Y;
            qsort(inverse, theProblem->size, sizeof(int), compare);
            for (i = 0; i < theProblem->size; i++)
                theProblem->number[inverse[i]] = i;
            free(inverse);  
            break;

        default : Error("Unexpected renumbering option"); }
}

int femDiffusionComputeBand(femDiffusionProblem *theProblem)
{
    femMesh *theMesh = theProblem->mesh;
    int iElem,j,myMax,myMin,myBand;
    int nLocal = theProblem->space->n;
    int map[nLocal];
    myBand = 0;
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j) 
            map[j] = theProblem->number[theProblem->map[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; }         
    return(++myBand);
}




void femDiffusionComputeMap(femDiffusionProblem *theProblem, int orderType)
{
    femMesh  *theMesh = theProblem->mesh;
    femEdges *theEdges = theProblem->edges;
       
    int i,j,k, iElem;
    int nElem = theMesh->nElem;
    int nNode = theMesh->nNode;
    int nEdge = theEdges->nEdge;
    int nLocal = theProblem->space->n;
    int nLocalNode = theProblem->mesh->nLocalNode;
    double xsi[16],eta[16],phi[4];
    femDiscreteXsi2(theProblem->space,xsi,eta);
    femDiscrete *theLinearSpace;
    if (nLocalNode == 3) theLinearSpace = femDiscreteCreate(3,FEM_TRIANGLE);
    else		theLinearSpace = femDiscreteCreate(4,FEM_QUAD);	
    
    if (orderType >= 1) {
        for (i=0; i < theMesh->nNode; i++) {
            theProblem->X[i] = theMesh->X[i];
            theProblem->Y[i] = theMesh->Y[i]; }
        for (i=0; i < nElem; i++) 
            for (j=0; j < nLocalNode; j++)
                theProblem->map[i*nLocal+j] = theProblem->mesh->elem[i*nLocalNode+j]; }
                
    if (orderType == 2) {
        for (i=0; i < theEdges->nEdge; i++) {
            int node0 = theEdges->edges[i].node[0];
            int node1 = theEdges->edges[i].node[1];
            theProblem->X[i+nNode] = (theMesh->X[theEdges->edges[i].node[0]]
                     									+ theMesh->X[theEdges->edges[i].node[1]])/2.0;
            theProblem->Y[i+nNode] = (theMesh->Y[theEdges->edges[i].node[0]]
                                                        + theMesh->Y[theEdges->edges[i].node[1]])/2.0;
            for (j=0; j < 2; j++) {
                iElem = theEdges->edges[i].elem[j];
                if (iElem >= 0) {
                    int *elem = &theProblem->mesh->elem[iElem*nLocalNode];
                    for (k=0; k < nLocalNode; k++) {
                        int k1 = k+1;
                        if (k1 == nLocalNode) k1 = 0;
                        if (node0 == elem[k] && node1 == elem[k1]) theProblem->map[iElem*nLocal+nLocalNode+k] = i+nNode;
                        if (node1 == elem[k] && node0 == elem[k1]) theProblem->map[iElem*nLocal+nLocalNode+k] = i+nNode; }}}}}
                        
    if (orderType == 2 && nLocalNode ==  4) {
            for (i=0; i < nElem; i++) {
            	theProblem->map[i*nLocal+8] = i+nNode+nEdge;
                for (j=0; j < nLocalNode; j++) {
                	theProblem->X[i+nNode+nEdge] = 0.0;
                	theProblem->Y[i+nNode+nEdge] = 0.0; }
             	for (j=0; j < nLocalNode; j++) {
					theProblem->X[i+nNode+nEdge] += theMesh->X[theMesh->elem[i*nLocalNode+j]]/4.0;
 					theProblem->Y[i+nNode+nEdge] += theMesh->Y[theMesh->elem[i*nLocalNode+j]]/4.0; }}}
               
                    
    if (orderType == 3) {
        for (i=0; i < theEdges->nEdge; i++) {
            int node0 = theEdges->edges[i].node[0];
            int node1 = theEdges->edges[i].node[1];
            theProblem->X[i*2+nNode] = (2.0 * theMesh->X[theEdges->edges[i].node[0]]
                     									+ theMesh->X[theEdges->edges[i].node[1]])/3.0;
            theProblem->Y[i*2+nNode] = (2.0 * theMesh->Y[theEdges->edges[i].node[0]]
                                                        + theMesh->Y[theEdges->edges[i].node[1]])/3.0;
            theProblem->X[i*2+nNode+1] = (theMesh->X[theEdges->edges[i].node[0]]
                     									+ 2.0 * theMesh->X[theEdges->edges[i].node[1]])/3.0;
            theProblem->Y[i*2+nNode+1] = (theMesh->Y[theEdges->edges[i].node[0]]
                                                        + 2.0 * theMesh->Y[theEdges->edges[i].node[1]])/3.0;

            for (j=0; j < 2; j++) {
                iElem = theEdges->edges[i].elem[j];
                if (iElem >= 0) {
                    int *elem = &theProblem->mesh->elem[iElem*nLocalNode];
                    for (k=0; k < nLocalNode; k++) {
                        int k1 = k+1;
                        if (k1 == nLocalNode) k1 = 0;
                        if (node0 == elem[k] && node1 == elem[k1]) {
                        		theProblem->map[iElem*nLocal+nLocalNode+k*2] = i*2+nNode;
                                theProblem->map[iElem*nLocal+nLocalNode+k*2+1] = i*2+nNode + 1; }
                        if (node1 == elem[k] && node0 == elem[k1]) {
                        		theProblem->map[iElem*nLocal+nLocalNode+k*2] = i*2+nNode + 1;
                                theProblem->map[iElem*nLocal+nLocalNode+k*2+1] = i*2+nNode; }}}}}}
                                
     if (orderType == 3 && nLocalNode ==  3) {
            for (i=0; i < nElem; i++) {
            	theProblem->map[i*nLocal+9] = i+nNode+2*nEdge;
                for (j=0; j < nLocalNode; j++) {
                	theProblem->X[i+nNode+nEdge*2] = 0.0;
                	theProblem->Y[i+nNode+nEdge*2] = 0.0; }
             	for (j=0; j < nLocalNode; j++) {
					theProblem->X[i+nNode+nEdge*2] += theMesh->X[theMesh->elem[i*nLocalNode+j]]/3.0;
 					theProblem->Y[i+nNode+nEdge*2] += theMesh->Y[theMesh->elem[i*nLocalNode+j]]/3.0; }}}
    
    if (orderType == 3 && nLocalNode ==  4) {
            for (i=0; i < nElem; i++) {
            	for (k=0; k < 4; k++) {
                    femDiscretePhi2(theLinearSpace, xsi[12+k], eta[12+k], phi);
            		theProblem->map[i*nLocal+12+k] = i*4+nNode+2*nEdge+k;
                		for (j=0; j < nLocalNode; j++) {
                			theProblem->X[i*4+nNode+nEdge*2+k] = 0.0;
                			theProblem->Y[i*4+nNode+nEdge*2+k] = 0.0; }
             			for (j=0; j < nLocalNode; j++) {
							theProblem->X[i*4+nNode+nEdge*2+k] += phi[j] * theMesh->X[theMesh->elem[i*nLocalNode+j]];
 							theProblem->Y[i*4+nNode+nEdge*2+k] += phi[j] * theMesh->Y[theMesh->elem[i*nLocalNode+j]]; }
                    
                    
                    }}}
                            
                                
                femDiscreteFree(theLinearSpace);
                
    femDiffusionPrintMap(theProblem);
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


void femDiffusionPrintInfos(femDiffusionProblem *theProblem)
{
    int  size = theProblem->size;
    printf(" \n");
    printf("    Diffusion problem \n");
    printf("    Number of nodal unknowns      : %8d\n",size);
}


femAdvProblem *femAdvCreate(const char *meshFileName)
{
    femAdvProblem *myProblem = malloc(sizeof(femAdvProblem));
            
    myProblem->mesh = femMeshRead(meshFileName);
    myProblem->edges = femEdgesCreate(myProblem->mesh);
    myProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
    myProblem->rule1d = femIntegrationCreate(2,FEM_EDGE);
    myProblem->rule2d = femIntegrationCreate(3,FEM_TRIANGLE);
    
    int sizeLoc = myProblem->space->n;
    int sizeGlo = myProblem->mesh->nElem * sizeLoc + 1; 
    
    myProblem->size = sizeGlo;
    myProblem->C = malloc(sizeof(double)*sizeGlo);
    myProblem->U = malloc(sizeof(double)*sizeGlo);
    myProblem->V = malloc(sizeof(double)*sizeGlo);
    myProblem->F = malloc(sizeof(double)*sizeGlo);
    
	return myProblem;
}
    
void femAdvFree(femAdvProblem *myProblem)
{
    free(myProblem->F);
    free(myProblem->U);
    free(myProblem->V);
    free(myProblem->C);
    femIntegrationFree(myProblem->rule1d);
    femIntegrationFree(myProblem->rule2d);
    femDiscreteFree(myProblem->space);
    femEdgesFree(myProblem->edges);
    femMeshFree(myProblem->mesh);
    free(myProblem);
}


void femStommel(double x, double y, double *u, double *v, double *eta)
{
    //
    // Solution analytique de Stommel dans un carre [0,1]x[0,1]
    // Modelisation de l'elevation de l'ocean Atlantique dans un carre adimensionnel
    // Ce modele que l'on attribue generalement au grand oceanographe Henry M.
    // Stommel (1920-1992), est considere comme le premier modele qualitativement correct du Gulf Stream
    //
    
    const double tau0 = 0.1;
    const double L = 1e6;
    const double gamm = 1e-6;
    const double rho = 1000;
    const double delta = 1;
    const double g = 9.81;
    const double h = 1000;
    const double f0 = 1e-4;
    const double beta = 0.5e-11;
    
    double Y = y - 0.5;
    double epsilon = gamm / (L * beta);
    double Z1 = (-1 + sqrt(1 + (2 * M_PI * delta * epsilon) * (2 * M_PI * delta * epsilon))) / (2 * epsilon);
    double Z2 = (-1 - sqrt(1 + (2 * M_PI * delta * epsilon) * (2 * M_PI * delta * epsilon))) / (2 * epsilon);
    double D = ((exp(Z2) - 1) * Z1 + (1 - exp(Z1)) * Z2) / (exp(Z1) - exp(Z2));
    double f1 = M_PI / D * (1 + ((exp(Z2) - 1) * exp(x * Z1) + (1 - exp(Z1)) * exp(x * Z2)) / (exp(Z1) - exp(Z2)));
    double f2 = 1 / D* (((exp(Z2) - 1) * Z1 * exp(x * Z1) + (1 - exp(Z1)) * Z2 * exp(x * Z2)) / (exp(Z1) - exp(Z2)));
    
    eta[0] = D * tau0 * f0 * L / (M_PI * gamm * rho * delta * g * h) *
               ( - gamm / (f0 * delta * M_PI) * f2 * sin(M_PI * Y) 
                 + 1 / M_PI * f1 * (cos(M_PI * Y) * (1 + beta * Y) 
                 - beta / M_PI * sin(M_PI * Y) ) );
    u[0] = D * tau0 / (M_PI * gamm * rho * h) * f1 * sin(M_PI * Y);
    v[0] = D * tau0 / (M_PI * gamm * rho * delta * h) * f2 * cos(M_PI * Y);
}


double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n)
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}

void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}


