#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <GL/glfw.h>
#define GLFW_INCLUDE_GLU

#include <pthread.h>

#define Error(a)   femError(a,__LINE__,__FILE__)
#define Warning(a) femWarning(a,  __LINE__, __FILE__)
#define NR_THREADS 8

static double R = 6371220;

/*
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 Struct
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 */

typedef enum {FEM_TRIANGLE,FEM_QUAD,FEM_EDGE} femElementType;
typedef enum {FEM_FULL,FEM_BAND,FEM_ITER} femSolverType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;

typedef struct {
    int *elem;
    double *X;
    double *Y;
	double *Z;
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;

typedef struct {
    int elem[2];
    int node[2];
} femEdge;

typedef struct {
    femMesh *mesh;
    femEdge *edges;
    int nEdge;
    int nBoundary;
} femEdges;

typedef struct {
    int n;
    int order;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x1)(double *xsi);
    void (*phi1)(double xsi, double *phi);
    void (*dphi1dx)(double xsi, double *dphidxsi);
	
} femDiscrete;

typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    femSolverType type;
    void *solver;
} femSolver;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

typedef struct
{
    double *B;
    double **A;
    int size;
    int band;
} femBandSystem;

typedef struct
{
    double *R;
    double *R0;
    double *D;
    double *D0;
    double error;
    int size;
    int iter;
} femIterativeSolver;

typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femIntegration *rule1d;
    femIntegration *rule2d;
    femSolver *solver;
    int size;
    double *E;
    double *U;
    double *V;
    double *FE;
    double *FU;
    double *FV;
    double omega;
    double gravity;
    double gamma;
    double timeStep;
} femShallowProblem;


typedef struct threadInfo {
	
	int threadNr;
	femShallowProblem *theProblem;
	
} threadInfo;

static femShallowProblem *theProblem;

/*
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 Protoypes
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 */

femShallowProblem *femShallowCreate(const char *meshFileName);
femMesh *femMeshRead(const char *filename);
femEdges *femEdgesCreate(femMesh *theMesh);
femIntegration *femIntegrationCreate(int n, femElementType type);
femDiscrete *femDiscreteCreate(int n, femElementType type);
femSolver *femSolverFullCreate(int size);
femSolver *femSolverBandCreate(int size, int band);
void femError(char *text, int line, char *file);
void femWarning(char *text, int line, char *file);
int femEdgesCompare(const void *edgeOne, const void *edgeTwo);
femFullSystem *femFullSystemCreate(int size);
femBandSystem *femBandSystemCreate(int size, int band);
void femFullSystemInit(femFullSystem *mySystem);
void femBandSystemInit(femBandSystem *myBandSystem);
void femShallowFree(femShallowProblem *myProblem);
void femSolverFree(femSolver *mySolver);
void femIntegrationFree(femIntegration *theRule);
void femDiscreteFree(femDiscrete *theSpace);
void femEdgesFree(femEdges *theEdges);
void femMeshFree(femMesh *theMesh);
void femFullSystemFree(femFullSystem *theSystem);
void femBandSystemFree(femBandSystem *myBandSystem);
void femIterativeSolverFree(femIterativeSolver *mySolver);
double interpolate(double *phi, double *U, int *map, int n);
void* femShallowAddIntegralsElements(void* args);
void* femShallowAddIntegralsEdges(void* args);
void femShallowMultiplyInverseMatrix(femShallowProblem *myProblem);
//void femShallowComputeEuler(femShallowProblem *theProblem);
void femShallowComputeHeun(femShallowProblem *theProblem);
void femDiscretePhi1(femDiscrete* mySpace, double xsi, double *phi);
void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
void femEdgesMap(femEdges *theEdges, int index, int map[2][2]);
void femSolverInit(femSolver *mySolver);
void femFullSystemInit(femFullSystem *mySystem);
void femBandSystemInit(femBandSystem *myBandSystem);
void femIterativeSolverInit(femIterativeSolver *mySolver);
void femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc);
void femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc);
void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double *femSolverEliminate(femSolver *mySolver);
double* femFullSystemEliminate(femFullSystem *mySystem);
double  *femBandSystemEliminate(femBandSystem *myBand);
double *femIterativeSolverEliminate(femIterativeSolver *mySolver);
void tsunamiApplyInitialCondition(femShallowProblem *theProblem);
void tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub);
int tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem);
void* femComputeThread(void* args);


/*
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 FEM functions
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 */

femShallowProblem *femShallowCreate(const char *meshFileName)
{
    femShallowProblem *myProblem = malloc(sizeof(femShallowProblem));
	
	
	// More realistic parameters, if you want to try
	
	myProblem->omega   = 2*M_PI/86400.0;
	myProblem->gravity = 9.81;
	myProblem->gamma   = 1e-7;
    
    myProblem->mesh = femMeshRead(meshFileName);
    myProblem->edges = femEdgesCreate(myProblem->mesh);
	
    
    myProblem->rule1d = femIntegrationCreate(2,FEM_EDGE);
    if (myProblem->mesh->nLocalNode == 4) {
        myProblem->space = femDiscreteCreate(4,FEM_QUAD);
        myProblem->rule2d = femIntegrationCreate(9,FEM_QUAD);  }
    else if (myProblem->mesh->nLocalNode == 3) {
		//printf("BITCH");
        myProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        myProblem->rule2d = femIntegrationCreate(3,FEM_TRIANGLE); }
	
	int i;
    int sizeLoc = myProblem->space->n;
    int sizeGlo = myProblem->mesh->nElem * sizeLoc + 1;
	//printf("sizeGlo %d\n", sizeGlo);
    myProblem->solver = femSolverFullCreate(3*sizeLoc);
    
    myProblem->size = sizeGlo;
    myProblem->E = malloc(sizeof(double)*sizeGlo);
    myProblem->U = malloc(sizeof(double)*sizeGlo);
    myProblem->V = malloc(sizeof(double)*sizeGlo);
    myProblem->FE = malloc(sizeof(double)*sizeGlo);
    myProblem->FU = malloc(sizeof(double)*sizeGlo);
    myProblem->FV = malloc(sizeof(double)*sizeGlo);
    for (i=0; i < myProblem->size; i++) {
        myProblem->E[i] = 0.0;
        myProblem->U[i] = 0.0;
        myProblem->V[i] = 0.0;  }
	
	
    return myProblem;
}

femMesh *femMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));
	
    int i,trash,*elem;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");
	
    fscanf(file, "Number of nodes %d \n", &theMesh->nNode);
	//printf("Number of nodes %d \n", theMesh->nNode);
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
	theMesh->Z = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fscanf(file,"%d : %le %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i],&theMesh->Z[i]); }
    
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
    //femMeshWrite(theMesh,"toto.txt");
	//printf("theMesh->nLocalNode %d\n",theMesh->nLocalNode);
	return theMesh;
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
		//printf("theSpace-> %d\n",theSpace->n);
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

void femDiscretePhi1(femDiscrete* mySpace, double xsi, double *phi)
{
    mySpace->phi1(xsi,phi);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
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

femFullSystem *femFullSystemCreate(int size)
{
	//printf("size full system %d\n", size);
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

void femEdgesMap(femEdges *theEdges, int index, int map[2][2])
{
    int i,j,k;
    int n = theEdges->mesh->nLocalNode;
    
    for (j=0; j < 2; ++j) {
        int node = theEdges->edges[index].node[j];
        for (k=0; k < 2; k++) {
            int elem = theEdges->edges[index].elem[k];
            map[k][j] = (theEdges->mesh->nElem)*n;
            if (elem >= 0) {
                for (i=0; i < n; i++) {
                    if (theEdges->mesh->elem[elem*n + i] == node) {
                        map[k][j] = elem*n + i;  }}}}}
}

void femSolverInit(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemInit((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemInit((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverInit((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

void femFullSystemInit(femFullSystem *mySystem)
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++)
        mySystem->B[i] = 0;
}

void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++)
        myBandSystem->B[i] = 0;
}

void femIterativeSolverInit(femIterativeSolver *mySolver)
{
    int i;
    mySolver->iter = 0;
    mySolver->error = 10.0e+12;
    for (i=0 ; i < mySolver->size*4 ; i++)
        mySolver->R[i] = 0;
}

void femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemAssemble((femFullSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_BAND : femBandSystemAssemble((femBandSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_ITER : femIterativeSolverAssemble((femIterativeSolver *)mySolver->solver,Aloc,Bloc,Uloc,map,nLoc); break;
        default : Error("Unexpected solver type"); }
}

void femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) {
        for(j = 0; j < nLoc; j++) {
            mySystem->A[map[i]][map[j]] += Aloc[i*nLoc+j]; }
		mySystem->B[map[i]] += Bloc[i]; }
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


/*
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 TSUNAMI functions
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 */

void tsunamiCompute(double dt, int nmax, int sub, int order, const char *meshFileName, const char *baseResultName){
	
	theProblem = femShallowCreate(meshFileName);
	tsunamiApplyInitialCondition(theProblem);
	theProblem->timeStep = dt;
	tsunamiWriteFile(baseResultName, 0, theProblem->U, theProblem->V, theProblem->E, theProblem->mesh->nElem, theProblem->mesh->nLocalNode);
	
	int i;
	for (i=1; i<=nmax; i++) {
		femShallowComputeHeun(theProblem);
		if (i%sub == 0) {
			tsunamiWriteFile(baseResultName, i, theProblem->U, theProblem->V, theProblem->E, theProblem->mesh->nElem, theProblem->mesh->nLocalNode);
		}
	}
	
	femShallowFree(theProblem);
	
}

void tsunamiAnimate(double dt, int nmax, int sub, int order, const char *meshFileName, const char *baseResultName)
{
	int nElem,nNode,i,j,index,trash,*elem;
    double *X,*Y,*H,*E,*U,*V;
    int width,height,mouse;
    
    int nout = sub;
    double t;
    double R = 6371220;
    double BathMax = 9368;
    GLfloat colors[9], coord[9];
	
    FILE* file = fopen(meshFileName,"r");
    if (file == NULL) {
    	printf("Error : cannot open mesh file :-) \n");
        exit(0); }
 	fscanf(file, "Number of nodes %d \n", &nNode);
  	X = malloc(sizeof(double)*nNode);
  	Y = malloc(sizeof(double)*nNode);
  	H = malloc(sizeof(double)*nNode);
	for (i = 0; i < nNode; i++)
    	fscanf(file,"%d : %le %le %le  \n",&trash,&X[i],&Y[i],&H[i]);
    fscanf(file, "Number of triangles %d \n", &nElem);
  	elem = malloc(sizeof(int)*3*nElem);
 	U = malloc(sizeof(double)*3*nElem);
 	V = malloc(sizeof(double)*3*nElem);
 	E = malloc(sizeof(double)*3*nElem);
  	for (i = 0; i < nElem; i++)
    	fscanf(file,"%d : %d %d %d \n", &trash,&elem[i*3],&elem[i*3+1],&elem[i*3+2]);
  	fclose(file);
    
	
   	glfwInit();
   	glfwOpenWindow(640,480,0,0,0,0,1,0,GLFW_WINDOW );
	glfwSetWindowTitle( "MECA1120 Tsunami" );
    
    glfwEnable( GLFW_STICKY_KEYS );
    glfwSwapInterval( 1 );
	
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 0.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 8.0, 8.0, 8.0, 0.0 };
    
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    GLfloat light_radiance[] = {1., 1., 1., 1.};
	
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_radiance);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_radiance);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
	
    double t0 = 0;
    int frame0,frame = -1;
	
    do {
        t = glfwGetTime();
        frame0 = frame;
    	frame = (int) ((t-t0) * 2);
		
        if (frame0 != frame) {
			
            glfwGetMousePos( &mouse, NULL );  mouse = 389;
			
            char filename[256];
            const char *basename = "data/tsunami-%08d.txt";
            sprintf(filename, basename, frame * nout);
			// A decommenter lorsque vous avez des fichiers de resultats
			
			if (!access(filename, F_OK)) {
				glfwTerminate(); exit( EXIT_SUCCESS ); }
			
            printf("===  Reading local file %s %d %f \n",filename,frame,t);
            tsunamiReadFile(baseResultName,frame*nout,U,V,E,nElem);
			
            glfwGetWindowSize( &width, &height );
            height = height > 0 ? height : 1;
            glViewport( 0, 0, width, height );
			
            glClearColor( 0.9f, 0.9f, 0.8f, 0.0f );
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluPerspective(65.0f,(GLfloat)width/(GLfloat)height,1.0f,100.0f);
			
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            gluLookAt(0.0f,1.0f,0.0f,0.0f, 20.0f, 0.0f,0.0f,0.0f,1.0f);
            glTranslatef(0.0f,14.0f,0.0f);
            double tt = 0;
            glRotatef(0.3f*(GLfloat)mouse + (GLfloat)tt*10.0f,0.0f,0.0f,1.0f);
            
            GLUquadricObj *quadratic = gluNewQuadric();
            gluQuadricNormals(quadratic, GLU_SMOOTH);
            glColor3f(1.0,1.0,1.0);
            gluSphere(quadratic,5.95,400,200);
            // A commenter pour supprimer la sphere interieure
            // Conseille pour les maillages grossiers :-)
			
            for (i=0; i < nElem; ++i) {
                for (j=0; j < 3; ++j) {
                    index = elem[3*i+j];
                    double value = H[index]/BathMax;
                    value = E[3*i+j]*10;
                    if (value < 0) value = 0;
                    if (value > 1) value = 1;
                    colors[j*3+0] = 3.5*(value)*(value);
                    colors[j*3+1] = (1-value)*(value)*3.5;
                    colors[j*3+2] = (1-value)*(1-value);
                    double x = X[index];
                    double y = Y[index];
                    double Factor = (4*R*R + x*x + y*y)*(R/6);
                    coord[j*3+0] = 4*R*R * x / Factor;
                    coord[j*3+1] = 4*R*R * y / Factor;
                    coord[j*3+2] = (4*R*R - x*x - y*y)*R / Factor;  }
				
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glNormalPointer(GL_FLOAT, 0, coord);
                glColorPointer(3, GL_FLOAT, 0, colors);
                glDrawArrays(GL_TRIANGLES, 0, 3);
                glDisableClientState(GL_NORMAL_ARRAY);
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_VERTEX_ARRAY);
				
                glColor3f(0.0, 0.0, 0.0);
                glEnableClientState(GL_VERTEX_ARRAY);
                for (j=0; j < 9; ++j)
					coord[j] = coord[j] * 1.001;
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glDrawArrays(GL_LINE_LOOP, 0, 3);
                glDisableClientState(GL_VERTEX_ARRAY); }
            glfwSwapBuffers();
			
		}} while( glfwGetKey( GLFW_KEY_ESC ) != GLFW_PRESS && glfwGetWindowParam( GLFW_OPENED ) );
    
    glfwTerminate();
    exit( EXIT_SUCCESS );
	
}

double tsunamiInitialConditionOkada(double x, double y)
{
	double x3d = 4*R*R*x / (4*R*R + x*x + y*y);
    double y3d = 4*R*R*y / (4*R*R + x*x + y*y);
    double z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
    double lat = asin(z3d/R)*180/M_PI;
    double lon = atan2(y3d,x3d)*180/M_PI;
    double lonMin = 142;
    double lonMax = 143.75;
    double latMin = 35.9;
    double latMax = 39.5;
    double olon = (lonMin+lonMax)/2.0;
    double olat = (latMin+latMax)/2.0;
    double angle = -12.95*M_PI/180;
    double lon2 = olon + (lon-olon)*cos(angle) + (lat-olat)*sin(angle);
    double lat2 = olat - (lon-olon)*sin(angle) + (lat-olat)*cos(angle);
	if ( lon2 <= lonMax && lon2 >= lonMin && lat2 >= latMin && lat2 <= latMax )
		return 1.0;
    else    return 0.0;
}

void tsunamiApplyInitialCondition(femShallowProblem *theProblem){
	//Interpolation? Qu'est ce qu'il veut? A voir plus tard!
	
	double *E = theProblem->E;
	int n = theProblem->space->n;
	
	int elem,j,mapElem[n];
	double Xloc[n],Yloc[n],Zloc[n];
	
	for (elem=0; elem < theProblem->mesh->nElem; elem++) {
        int *mapCoord = &(theProblem->mesh->elem[elem*n]);
        for (j=0; j < n; ++j) {
            mapElem[j] = elem*n + j;
        	Xloc[j] = theProblem->mesh->X[mapCoord[j]];
        	Yloc[j] = theProblem->mesh->Y[mapCoord[j]];
			Zloc[j] = theProblem->mesh->Z[mapCoord[j]];
			E[mapElem[j]] = tsunamiInitialConditionOkada(Xloc[j],Yloc[j]);
		}
	}
}

void tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub) {
    int i,j;
    const char *basename = "%s-%08d.txt";
    char filename[256];
    sprintf(filename,basename,baseResultName,iter);
    FILE* file = fopen(filename,"w");
    fprintf(file, "Number of elem %d \n", nelem);
    fprintf(file, "Number of local values per element %d \n", nsub);
    for (i = 0; i < nelem; ++i) {
    	for (j = 0; j < nsub; ++j) {
        	int index = i*nsub+j;
        	fprintf(file,"%d;%d;%le;%le;%le;\n",i,j,U[index],V[index],E[index]); }}
    fclose(file);
}

int tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem)
{
    int i,j,trash,nelemFile,nsub;
    const char *basename = "%s-%08d.txt";
    char filename[256];
    sprintf(filename,basename,baseResultName,iter);
    FILE* file = fopen(filename,"r");
    fscanf(file, "Number of elem %d \n", &nelemFile);
    fscanf(file, "Number of local values per element %d \n", &nsub);
    if (nelem != nelemFile) {
        printf("Error : wrong data file %d %d:-) \n",nelem,nelemFile);
        exit(0); }
    for (i = 0; i < nelem; ++i) {
    	for (j = 0; j < nsub; ++j) {
        	int index = i*nsub+j;
        	fscanf(file,"%d;%d;%le;%le;%le;\n",&trash,&trash,&U[index],&V[index],&E[index]); }}
	
    fclose(file);
    return nsub;
}

/*
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 FEM computation
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 */

double interpolate(double *phi, double *U, int *map, int n)
{
    double u = 0.0; int i;
    for (i=0; i <n; i++)
        u += phi[i]*U[map[i]];
    return u;
}

void* femShallowAddIntegralsElements(void* args)
{
	
	threadInfo* tI = (threadInfo*)args;
	femShallowProblem *myProblem = tI->theProblem;
	int threadNr = tI->threadNr;
	
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
    double *Y = myProblem->mesh->Y;
	double *X = myProblem->mesh->X;
	double *Z = myProblem->mesh->Z;
    femIntegration *theRule = myProblem->rule2d;
    femDiscrete *theSpace = myProblem->space;
    int n = theSpace->n;
    double Xloc[n],Yloc[n],phi[n],dphidxsi[n],dphideta[n],dphidx[n],dphidy[n];
    double  xsi,eta,weight,jac;
    double  x,y,h,e,u,v,t,ge;
    int     i,j,k,elem,mapElem[n];
	
    double  g     = myProblem->gravity;
    double  gamma = myProblem->gamma;
    double  omega    = myProblem->omega;
	
	int nElem = myProblem->mesh->nElem;
	
	//printf("nElem : %d   ----> %d to %d  \n", nElem,(threadNr*nElem)/NR_THREADS ,((threadNr+1)*nElem)/NR_THREADS);
    for (elem=(threadNr*nElem)/NR_THREADS; elem < ((threadNr+1)*nElem)/NR_THREADS; elem++) {
        int *mapCoord = &(myProblem->mesh->elem[elem*n]);
        for (j=0; j < n; ++j) {
            mapElem[j] = elem*n + j;
        	Xloc[j] = myProblem->mesh->X[mapCoord[j]];
        	Yloc[j] = myProblem->mesh->Y[mapCoord[j]]; }
        for (k=0; k < theRule->n; k++) {
            xsi = theRule->xsi[k];
            eta = theRule->eta[k];
            weight = theRule->weight[k];
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            for (i = 0; i < n; i++) {
                dxdxsi += Xloc[i]*dphidxsi[i];
                dxdeta += Xloc[i]*dphideta[i];
                dydxsi += Yloc[i]*dphidxsi[i];
                dydeta += Yloc[i]*dphideta[i]; }
            jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }
            
			x = interpolate(phi,X,mapCoord,n);
            y = interpolate(phi,Y,mapCoord,n);
			h = interpolate(phi,Z,mapCoord,n);
			
			double z3d = R*(4.0*R*R - x*x - y*y) / (4.0*R*R + x*x + y*y);
			double lat = asin(z3d/R);
			
            const double coriolis = 2.0*omega*sin(lat);
            e = interpolate(phi,E,mapElem,n);
            u = interpolate(phi,U,mapElem,n);
            v = interpolate(phi,V,mapElem,n);
			
            t = (4.0 * R*R + x*x + y*y) / (4.0*R*R);
			ge = (g*e)/(2.0*R*R);
			
            for (i=0; i < n; i++) {
                BE[mapElem[i]] += ((h*u*dphidx[i] + h*v*dphidy[i]) * t)*jac*weight;
                BU[mapElem[i]] += ((coriolis*v - gamma*u)*phi[i] + g*e*dphidx[i]*t)*jac*weight;
                BV[mapElem[i]] += ((-coriolis*u - gamma*v)*phi[i] + g*e*dphidy[i]*t)*jac*weight;
				
				BE[mapElem[i]] += (phi[i]*(h*(x*u+y*v))/(R*R))*jac*weight;
                BU[mapElem[i]] += (phi[i]*x*ge)*jac*weight;
                BV[mapElem[i]] += (phi[i]*y*ge)*jac*weight;
			}
		}
	}
	return NULL;
	
}

void* femShallowAddIntegralsEdges(void* args)
{
    
	threadInfo* tI = (threadInfo*)args;
	femShallowProblem *myProblem = tI->theProblem;
	int threadNr = tI->threadNr;
	
    double  xEdge[2],yEdge[2],hEdge[2],phiEdge[2];
    double  xsi,weight,jac;
    double  eL,eR,uL,uR,vL,vR,unL,unR,x,y,h,t;
    double  qe,qu,qv;
    int     i,j,k,edge,mapEdge[2][2];
	
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
    femIntegration *theRule = myProblem->rule1d;
    femDiscrete *theSpace = myProblem->space;
    int 	sizeLoc = theSpace->n;
    int     sizeGlo = myProblem->mesh->nElem * sizeLoc + 1;
	
    double  g = myProblem->gravity;
	
	int nEdge = myProblem->edges->nEdge;
	//printf("nEdge : %d   ----> %d to %d  \n", nEdge,(threadNr*nEdge)/NR_THREADS ,((threadNr+1)*nEdge)/NR_THREADS);

    for (edge=(threadNr*nEdge)/NR_THREADS; edge < ((threadNr+1)*nEdge)/NR_THREADS; edge++) {
        femEdgesMap(myProblem->edges,edge,mapEdge);
        for (j=0; j < 2; ++j) {
        	int node = myProblem->edges->edges[edge].node[j];
        	xEdge[j] = myProblem->mesh->X[node];
        	yEdge[j] = myProblem->mesh->Y[node];
			hEdge[j] = myProblem->mesh->Z[node];}
        
        int boundary = (mapEdge[1][0] == sizeGlo-1);
        
        double dxdxsi = (xEdge[1] - xEdge[0]);
        double dydxsi = (yEdge[1] - yEdge[0]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        jac = norm / 2.0;
        for (k=0; k < theRule->n; k++) {
            xsi = theRule->xsi[k];
            weight = theRule->weight[k];
            femDiscretePhi1(theSpace,xsi,phiEdge);
			
            eL = interpolate(phiEdge,E,mapEdge[0],2);
            eR = boundary ? eL : interpolate(phiEdge,E,mapEdge[1],2);
            uL = interpolate(phiEdge,U,mapEdge[0],2);
            uR = interpolate(phiEdge,U,mapEdge[1],2);
            vL = interpolate(phiEdge,V,mapEdge[0],2);
            vR = interpolate(phiEdge,V,mapEdge[1],2);
			
			x = xEdge[0]*phiEdge[0]+xEdge[1]*phiEdge[1];
			y = yEdge[0]*phiEdge[0]+yEdge[1]*phiEdge[1];
			h = hEdge[0]*phiEdge[0]+hEdge[1]*phiEdge[1];
			
            unL = uL*nx+ vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;
            qe =  0.5*h*   ( (unL+unR) + sqrt(g/h)*( eL-eR ) );
            qu =  0.5*g*nx*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );
            qv =  0.5*g*ny*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );
			
            t = (4.0 * R*R + x*x + y*y) / (4.0*R*R);
			
            for (i=0; i < 2; i++) {
                BE[mapEdge[0][i]] -= qe*phiEdge[i]*t*jac*weight;
                BU[mapEdge[0][i]] -= qu*phiEdge[i]*t*jac*weight;
                BV[mapEdge[0][i]] -= qv*phiEdge[i]*t*jac*weight;
                BE[mapEdge[1][i]] += qe*phiEdge[i]*t*jac*weight;
                BU[mapEdge[1][i]] += qu*phiEdge[i]*t*jac*weight;
                BV[mapEdge[1][i]] += qv*phiEdge[i]*t*jac*weight;
			}
		}
	}
	
	return NULL;
}

void femShallowMultiplyInverseMatrix(femShallowProblem *myProblem)
{
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    femMesh *theMesh = myProblem->mesh;
    femDiscrete *theSpace = myProblem->space;
    femSolver *theSolver = myProblem->solver;
    femIntegration *theRule = myProblem->rule2d;
    
    int n = theSpace->n;
    double Xloc[n],Yloc[n],phi[n],dphidxsi[n],dphideta[n],Aloc[n*n],Uloc[n],jac;
    int iElem,iInteg,i,j,mapElem[n],mapE[n],mapU[n],mapV[n];
    
    
    for (i = 0; i < n; i++)   {
        Uloc[i] = 0;
        mapE[i] = i;
        mapU[i] = i + n;
        mapV[i] = i + 2*n; }
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femSolverInit(theSolver);
        for (i = 0; i < n*n; i++)  Aloc[i] = 0;
        int *mapCoord = &(myProblem->mesh->elem[iElem*n]);
        for (j=0; j < n; ++j) {
            mapElem[j] = iElem*n + j;
        	Xloc[j] = myProblem->mesh->X[mapCoord[j]];
        	Yloc[j] = myProblem->mesh->Y[mapCoord[j]]; }
		
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
            for (i = 0; i < n; i++) {
                dxdxsi += Xloc[i]*dphidxsi[i];
                dxdeta += Xloc[i]*dphideta[i];
                dydxsi += Yloc[i]*dphidxsi[i];
                dydeta += Yloc[i]*dphideta[i]; }
            jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                    Aloc[i*(theSpace->n)+j] += phi[i] * phi[j] * jac * weight; }}}
		
        femSolverAssemble(theSolver,Aloc,&BE[mapElem[0]],Uloc,mapE,theSpace->n);
        femSolverAssemble(theSolver,Aloc,&BU[mapElem[0]],Uloc,mapU,theSpace->n);
        femSolverAssemble(theSolver,Aloc,&BV[mapElem[0]],Uloc,mapV,theSpace->n);
        double *soluce = femSolverEliminate(theSolver);
     	for (i = 0; i < n; i++) {
        	BE[mapElem[i]] = soluce[i];
            BU[mapElem[i]] = soluce[i+n];
            BV[mapElem[i]] = soluce[i+2*n]; }}
	
}

/*void femShallowComputeEuler(femShallowProblem *theProblem)
 {
 
 int  size = theProblem->size, i;
 double *FE = theProblem->FE;
 double *FU = theProblem->FU;
 double *FV = theProblem->FV;
 double *E = theProblem->E;
 double *U = theProblem->U;
 double *V = theProblem->V;
 double theTimeStep = theProblem->timeStep;
 
 
 for (i=0; i < size; i++) {
 FE[i] = 0.0;
 FU[i] = 0.0;
 FV[i] = 0.0; }
 femShallowAddIntegralsElements(theProblem);
 femShallowAddIntegralsEdges(theProblem);
 femShallowMultiplyInverseMatrix(theProblem);
 for (i=0; i < size; i++) {
 E[i] += theTimeStep * FE[i];
 U[i] += theTimeStep * FU[i];
 V[i] += theTimeStep * FV[i]; }
 
 }*/

void femShallowComputeHeun(femShallowProblem *theProblem)
{
	
    int  size = theProblem->size, i;
    double *FE = theProblem->FE;
    double *FU = theProblem->FU;
    double *FV = theProblem->FV;
    double *E = theProblem->E;
    double *U = theProblem->U;
    double *V = theProblem->V;
    double theTimeStep = theProblem->timeStep;
	
	// For threads
	pthread_t threadIDs[NR_THREADS];
	
	double KE1[size],KU1[size],KV1[size],KE2[size],KU2[size],KV2[size];
    for (i=0; i < size; i++) {
        FE[i] = 0.0;
        FU[i] = 0.0;
        FV[i] = 0.0; }
	
	//------------------------------------------------------------------------------------
	// 1 - Creating the threads for Integration over elements
	
	for(i = 0; i < NR_THREADS; i++){
		// Fuck les free de ce truc, ça prend rien comme place !
		threadInfo* tI = malloc(sizeof(threadInfo));
		tI->theProblem = theProblem;
		tI->threadNr = i;
		if(pthread_create(&(threadIDs[i]), NULL, &femComputeThread, tI))
			printf("error in thread creation\n");
	}
	for(i = 0; i < NR_THREADS; i++){
		pthread_join(threadIDs[i], NULL);
	}
	
	
	// 3 - Inversion of matrix (not threaded)
    femShallowMultiplyInverseMatrix(theProblem);
	
    for (i=0; i < size; i++) {
		
		KE1[i] = FE[i];
		E[i] += theTimeStep * FE[i];
		FE[i] = 0;
        KU1[i] = FU[i];
		U[i] += theTimeStep * FU[i];
		FU[i] = 0;
        KV1[i] = FV[i];
		V[i] += theTimeStep * FV[i];
		FV[i] = 0;
	}
	
	
	//------------------------------------------------------------------------------------
	// 1' - Creating the threads for Integration over elements
	
	for(i = 0; i < NR_THREADS; i++){
		// Fuck les free de ce truc, ça prend rien comme place !
		threadInfo* tI = malloc(sizeof(threadInfo));
		tI->theProblem = theProblem;
		tI->threadNr = i;
		if(pthread_create(&(threadIDs[i]), NULL, &femComputeThread, tI))
			printf("error in thread creation\n");
	}
	for(i = 0; i < NR_THREADS; i++){
		pthread_join(threadIDs[i], NULL);
	}
	
	// 3' - Inversion of matrix (not threaded)
    femShallowMultiplyInverseMatrix(theProblem);
	
	for (i=0; i < size; i++) {
		
		KE2[i] = FE[i];
		E[i] += -theTimeStep * KE1[i] + theTimeStep/2.0 * (KE1[i] + KE2[i]);
        KU2[i] = FU[i];
		U[i] += -theTimeStep * KU1[i] + theTimeStep/2.0 * (KU1[i] + KU2[i]);
        KV2[i] = FV[i];
		V[i] += -theTimeStep * KV1[i] + theTimeStep/2.0 * (KV1[i] + KV2[i]);
	}
	
}

void* femComputeThread(void* args)
{
	femShallowAddIntegralsElements(args);
	femShallowAddIntegralsEdges(args);
	return NULL;
}


/*
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 FREE
 
 /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
 \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 */

void femShallowFree(femShallowProblem *myProblem)
{
    
    free(myProblem->FE);
    free(myProblem->FU);
    free(myProblem->FV);
    free(myProblem->E);
    free(myProblem->U);
    free(myProblem->V);
    femSolverFree(myProblem->solver);
    femIntegrationFree(myProblem->rule1d);
    femIntegrationFree(myProblem->rule2d);
    femDiscreteFree(myProblem->space);
    femEdgesFree(myProblem->edges);
    femMeshFree(myProblem->mesh);
    free(myProblem);
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
void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}

void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A);
    free(myBandSystem);
}

void femIterativeSolverFree(femIterativeSolver *mySolver)
{
    free(mySolver->R);
    free(mySolver);
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femEdgesFree(femEdges *theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->elem);
    free(theMesh);
}
