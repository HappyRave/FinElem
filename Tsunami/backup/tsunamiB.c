#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <GL/glfw.h>

/**
*   Projet - Simulation d'un tsunami
*   Mai 2013
*   Groupe 14
*   Alexis Godfrin & Benjamin Spitaels
*/

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Structures et noms des fonctions
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

#define Error(a)   femError(a,__LINE__,__FILE__)
#define Warning(a) femWarning(a,  __LINE__, __FILE__)

typedef enum {FEM_TRIANGLE,FEM_QUAD,FEM_EDGE} femElementType;
typedef enum {FEM_FULL,FEM_BAND,FEM_ITER} femSolverType;
//typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;

typedef struct { // Struture : grille des noeuds
    int *elem; // Sommets de chaque élément
    double *X; // Abscisses des noeuds
    double *Y; // Ordonnées des noeuds
    double *bath; // Bathymétrie
    int nElem; // Nombre total d'éléments
    int nNode; // Nombre total de noeuds
    int nLocalNode; // Nombre de noeuds par élément
} femMesh;

typedef struct { // Structure : segment
    int elem[2]; // Eléments contenant ce segment
    int node[2]; // noeuds du segment
} femEdge;

typedef struct { // Structure : tableau des segments
    femMesh *mesh; // grille
    femEdge *edges; // tableau des segments
    int nEdge; // Nombre total de segment
    int nBoundary; // Nombre de segments frontières
} femEdges;

typedef struct { // Structure : fonctions
    int n; // Nombre de noeuds de calcul
    int order; // Ordre des fonctions
    void (*x2)(double *xsi, double *eta); // Noeuds initiaux (2D)
    void (*phi2)(double xsi, double eta, double *phi); // Fonctions PHI (2D)
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta); // Dérivées des fonctions PHI (2D)
    void (*x1)(double *xsi); // Noeuds initiaux (1D)
    void (*phi1)(double xsi, double *phi); // Fonctions PHI (1D)
    void (*dphi1dx)(double xsi, double *dphidxsi); // Dérivées des fonctions PHI (1D)
} femDiscrete;

typedef struct { // Structure : intégration
    int n; // Nombre de points d'intégration
    const double *xsi; // Abscisses des points d'intégration
    const double *eta; // Ordonnées des points d'intégration
    const double *weight; // Poids d'intégration
} femIntegration;

typedef struct { // Structure : solver
    femSolverType type; // Type de solver
    void *solver; // Lien vers le solver
} femSolver;

typedef struct { // Structure : solveur plein
    double *B; // Vecteur B
    double **A; // Matrice A
    int size; // Taille
} femFullSystem;

typedef struct { // Structure : solveur bande
    double *B; // Vecteur B
    double **A; // Matrice A
    int size; // Taille
    int band; // Taille bande
} femBandSystem;

typedef struct { // Structure : solveur itératif
    double *R; // Vecteur R
    double *R0; // Vecteur R0
    double *D; // Vecteur D
    double *D0; // Vecteur D0
    double error; // Valeur de l'erreur
    int size; // Taille
    int iter; // Nombre d'itération
} femIterativeSolver;

typedef struct { // Structure : Tsunami
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
    int *elem; // extension pour l'ordre 2
    double *X; // extension pour l'ordre 2
    double *Y; // extension pour l'ordre 2
    double *bath; // extension pour l'ordre 2
    int nLocal; // Extension ordre 2
    int nNode; // Extension ordre 2
    double gravity;
    double gamma;
    double R;
    double Omega;
    double timeStep;
} femTsunamiProblem;

femIntegration*      femIntegrationCreate(int n, femElementType type);
void                 femIntegrationFree(femIntegration *theRule);

femMesh*             femMeshRead(const char *filename);
void                 femMeshFree(femMesh *theMesh);

femEdges*            femEdgesCreate(femMesh *theMesh);
void                 femEdgesFree(femEdges *theEdges);
void                 femEdgesPrint(femEdges *theEdges);
int                  femEdgesCompare(const void *edgeOne, const void *edgeTwo);
void                 femEdgesMap(femEdges *theEdges, int index, int map[2][2]);

femDiscrete*         femDiscreteCreate(int n, femElementType type);
void                 femDiscreteFree(femDiscrete* mySpace);
void                 femDiscretePrint(femDiscrete* mySpace);
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);
void                 femDiscreteXsi1(femDiscrete* mySpace, double *xsi);
void                 femDiscretePhi1(femDiscrete* mySpace, double xsi, double *phi);
void                 femDiscreteDphi1(femDiscrete* mySpace, double xsi, double *dphidxsi);
double               femDiscreteInterpolate(double *phi, double *U, int *map, int n);

femSolver*           femSolverFullCreate(int size);
femSolver*           femSolverBandCreate(int size,int band);
femSolver*           femSolverIterativeCreate(int size);
void                 femSolverFree(femSolver* mySolver);
void                 femSolverInit(femSolver* mySolver);
void                 femSolverPrint(femSolver* mySolver);
void                 femSolverPrintInfos(femSolver* mySolver);
double*              femSolverEliminate(femSolver* mySolver);
void                 femSolverConstrain(femSolver* mySolver, int myNode, double value);
void                 femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double               femSolverGet(femSolver* mySolver, int i, int j);
int                  femSolverConverged(femSolver *mySolver);

femFullSystem*       femFullSystemCreate(int size);
void                 femFullSystemFree(femFullSystem* mySystem);
void                 femFullSystemInit(femFullSystem* mySystem);
void                 femFullSystemPrint(femFullSystem* mySystem);
void                 femFullSystemPrintInfos(femFullSystem* mySystem);
double*              femFullSystemEliminate(femFullSystem* mySystem);
void                 femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);
void                 femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femFullSystemGet(femFullSystem* mySystem, int i, int j);

femBandSystem*       femBandSystemCreate(int size, int band);
void                 femBandSystemFree(femBandSystem* myBandSystem);
void                 femBandSystemInit(femBandSystem *myBand);
void                 femBandSystemPrint(femBandSystem *myBand);
void                 femBandSystemPrintInfos(femBandSystem *myBand);
double*              femBandSystemEliminate(femBandSystem *myBand);
void                 femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue);
void                 femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femBandSystemGet(femBandSystem* myBandSystem, int i, int j);

femIterativeSolver*  femIterativeSolverCreate(int size);
void                 femIterativeSolverFree(femIterativeSolver* mySolver);
void                 femIterativeSolverInit(femIterativeSolver* mySolver);
void                 femIterativeSolverPrint(femIterativeSolver* mySolver);
void                 femIterativeSolverPrintInfos(femIterativeSolver* mySolver);
double*              femIterativeSolverEliminate(femIterativeSolver* mySolver);
void                 femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue);
void                 femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double               femIterativeSolverGet(femIterativeSolver* mySolver, int i, int j);
int                  femIterativeSolverConverged(femIterativeSolver *mySolver);

femTsunamiProblem   *femTsunamiCreate(const char *meshFileName, int order);
void                 femTsunamiFree(femTsunamiProblem *myProblem);
void                 femTsunamiMeshOrder(femTsunamiProblem *myProblem);
void                 femTsunamiTriangleMap(femTsunamiProblem *myProblem, int index, int *map);
void                 femTsunamiEdgeMap(femTsunamiProblem *myProblem, int index, int map[2][2]);
void                 femTsunamiAddIntegralsElements(femTsunamiProblem *myProblem);
void                 femTsunamiAddIntegralsEdges(femTsunamiProblem *myProblem);
void                 femTsunamiMultiplyInverseMatrix(femTsunamiProblem *myProblem);
void                 femTsunamiCompute(femTsunamiProblem *myProblem);
double               tsunamiInitialConditionOkada(double x, double y);
void                 tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub);
int                  tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem);

double               femMin(double *x, int n);
double               femMax(double *x, int n);
void                 femError(char *text, int line, char *file);
void                 femWarning(char *text, int line, char *file);

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions femIntegration
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

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
    if (type == FEM_QUAD && n == 4)
    {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight;
    }
    else if (type == FEM_QUAD && n == 9)
    {
        theRule->n      = 9;
        theRule->xsi    = _gaussQuad9Xsi;
        theRule->eta    = _gaussQuad9Eta;
        theRule->weight = _gaussQuad9Weight;
    }
    else if (type == FEM_TRIANGLE && n == 3)
    {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight;
    }
    else if (type == FEM_TRIANGLE && n == 12)
    {
        theRule->n      = 12;
        theRule->xsi    = _gaussTri12Xsi;
        theRule->eta    = _gaussTri12Eta;
        theRule->weight = _gaussTri12Weight;
    }
    else if (type == FEM_EDGE && n == 2)
    {
        theRule->n      = 2;
        theRule->xsi    = _gaussEdge2Xsi;
        theRule->eta    = NULL;
        theRule->weight = _gaussEdge2Weight;
    }
    else Error("Cannot create such an integration rule !");
    return theRule;
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonction femMesh
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

femMesh *femMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash,*elem;

    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    fscanf(file, "Number of nodes %d \n", &theMesh->nNode);
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    theMesh->bath = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i)
    {
        fscanf(file,"%d : %le %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i],&theMesh->bath[i]);
    }

    char str[256]; fgets(str, sizeof(str), file);
    if (!strncmp(str,"Number of triangles",19))
    {
        sscanf(str,"Number of triangles %d \n", &theMesh->nElem);
        theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
        theMesh->nLocalNode = 3;
        for (i = 0; i < theMesh->nElem; ++i)
        {
            elem = &(theMesh->elem[i*3]);
            fscanf(file,"%d : %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2]);
        }
    }
    else if (!strncmp(str,"Number of quads",15))
    {
        sscanf(str,"Number of quads %d \n", &theMesh->nElem);
        theMesh->elem = malloc(sizeof(int)*4*theMesh->nElem);
        theMesh->nLocalNode = 4;
        for (i = 0; i < theMesh->nElem; ++i)
        {
            elem = &(theMesh->elem[i*4]);
            fscanf(file,"%d : %d %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2],&elem[3]);
        }
    }

    fclose(file);
    return theMesh;
}

void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->bath);
    free(theMesh->elem);
    free(theMesh);
}

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions femEdges
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

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

    for (i = 0; i < theMesh->nElem; i++)
    {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++)
        {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc];
        }
    }

    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), femEdgesCompare);

    int index = 0;
    int nBoundary = 0;

    for (i=0; i < theEdges->nEdge; i++)
    {
        if (i == theEdges->nEdge - 1 || femEdgesCompare(&edges[i],&edges[i+1]) != 0)
        {
            edges[index] = edges[i];
            nBoundary++;
        }
        else
        {
            edges[index] = edges[i];
            edges[index].elem[1] = edges[i+1].elem[0];
            i = i+1;
        }
        index++;
    }

    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
    return theEdges;
}

void femEdgesFree(femEdges *theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

void femEdgesPrint(femEdges *theEdges)
{
    int i;
    for (i = 0; i < theEdges->nEdge; ++i)
    {
        printf("%6d : %4d %4d : %4d %4d \n",i,
               theEdges->edges[i].node[0],theEdges->edges[i].node[1],
               theEdges->edges[i].elem[0],theEdges->edges[i].elem[1]);
    }
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

void femEdgesMap(femEdges *theEdges, int index, int map[2][2])
{
    int i,j,k;
    int n = theEdges->mesh->nLocalNode;

    for (j=0; j < 2; ++j)
    {
        int node = theEdges->edges[index].node[j];
        for (k=0; k < 2; k++)
        {
            int elem = theEdges->edges[index].elem[k];
            map[k][j] = (theEdges->mesh->nElem) * n;
            if (elem >= 0)
            {
                for (i=0; i < n; i++)
                {
                    if (theEdges->mesh->elem[elem*n + i] == node)
                    {
                        map[k][j] = elem*n + i;
                    }
                }
            }
        }
    }
}

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions femDiscrete
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

void _1c0_x(double *xsi) // EDGE - ORDER 1
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

void _q1c0_x(double *xsi, double *eta) // QUAD - ORDER 1
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

void _q2c0_x(double *xsi, double *eta) // QUAD - ORDER 2
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

void _p1c0_x(double *xsi, double *eta) //TRI - ORDER 1
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

void _p2c0_x(double *xsi, double *eta) // TRI - ORDER 2
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

femDiscrete *femDiscreteCreate(int n, femElementType type) // MODIFICATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_QUAD && n == 4)
    {
        theSpace->n       = 4;
        theSpace->order   = 1;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx;
        theSpace->x1      = _1c0_x;
        theSpace->phi1    = _1c0_phi;
        theSpace->dphi1dx = _1c0_dphidx;
    }
    else if (type == FEM_QUAD && n == 9)
    {
        theSpace->n       = 9;
        theSpace->order   = 2;
        theSpace->x2      = _q2c0_x;
        theSpace->phi2    = _q2c0_phi;
        theSpace->dphi2dx = _q2c0_dphidx;
    }
    else if (type == FEM_TRIANGLE && n == 3)
    {
        theSpace->n       = 3;
        theSpace->order   = 1;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx;
        theSpace->x1      = _1c0_x;
        theSpace->phi1    = _1c0_phi;
        theSpace->dphi1dx = _1c0_dphidx;
    }
    else if (type == FEM_TRIANGLE && n == 6)
    {
        theSpace->n       = 6;
        theSpace->order   = 2;
        theSpace->x2      = _p2c0_x;
        theSpace->phi2    = _p2c0_phi;
        theSpace->dphi2dx = _p2c0_dphidx;
    }
    else Error("Cannot create such a discrete space !");
    return theSpace;
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[16], eta[16], phi[16], dphidxsi[16], dphideta[16];

    femDiscreteXsi2(mySpace,xsi,eta);
    for (i=0; i < n; i++)
    {

        femDiscretePhi2(mySpace,xsi[i],eta[i],phi);
        femDiscreteDphi2(mySpace,xsi[i],eta[i],dphidxsi,dphideta);

        for (j=0; j < n; j++)
        {
            printf("(xsi=%+.1f,eta=%+.1f) : ",xsi[i],eta[i]);
            printf(" phi(%d)=%+.1f",j,phi[j]);
            printf("   dphidxsi(%d)=%+.1f",j,dphidxsi[j]);
            printf("   dphideta(%d)=%+.1f \n",j,dphideta[j]);
        }
        printf(" \n");
    }
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

double femDiscreteInterpolate(double *phi, double *U, int *map, int n)
{
    double u = 0.0; int i;
    for (i=0; i <n; i++)
    {
        u += phi[i]*U[map[i]];
    }
    return u;
}

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions femSolver
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

femSolver *femSolverFullCreate(int size) // FULL
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_FULL;
    mySolver->solver = (femSolver *)femFullSystemCreate(size);
    return(mySolver);
}

femSolver *femSolverBandCreate(int size, int band) // BANDE
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_BAND;
    mySolver->solver = (femSolver *)femBandSystemCreate(size,band);
    return(mySolver);
}

femSolver *femSolverIterativeCreate(int size) // ITER
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_ITER;
    mySolver->solver = (femSolver *)femIterativeSolverCreate(size);
    return(mySolver);
}

void femSolverFree(femSolver *mySolver) // TOUS
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemFree((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemFree((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverFree((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    free(mySolver);
}

void femSolverInit(femSolver *mySolver) // TOUS
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemInit((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemInit((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverInit((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

double femSolverGet(femSolver *mySolver,int i,int j) // TOUS
{
    double value = 0;
    switch (mySolver->type) {
        case FEM_FULL : value = femFullSystemGet((femFullSystem *)mySolver->solver,i,j); break;
        case FEM_BAND : value = femBandSystemGet((femBandSystem *)mySolver->solver,i,j); break;
        case FEM_ITER : value = (i==j); break;
        default : Error("Unexpected solver type"); }
    return(value);
}

void femSolverPrint(femSolver *mySolver) // TOUS
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrint((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrint((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrint((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

void femSolverPrintInfos(femSolver *mySolver) // TOUS
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrintInfos((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrintInfos((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrintInfos((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

void femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc) // TOUS
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemAssemble((femFullSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_BAND : femBandSystemAssemble((femBandSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        case FEM_ITER : femIterativeSolverAssemble((femIterativeSolver *)mySolver->solver,Aloc,Bloc,Uloc,map,nLoc); break;
        default : Error("Unexpected solver type"); }
}

void femSolverConstrain(femSolver *mySolver, int myNode, double myValue) // TOUS
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemConstrain((femFullSystem *)mySolver->solver,myNode,myValue); break;
        case FEM_BAND : femBandSystemConstrain((femBandSystem *)mySolver->solver,myNode,myValue); break;
        case FEM_ITER : femIterativeSolverConstrain((femIterativeSolver *)mySolver->solver,myNode,myValue); break;
        default : Error("Unexpected solver type"); }
}

double *femSolverEliminate(femSolver *mySolver) // TOUS
{
    double *soluce;
    switch (mySolver->type) {
        case FEM_FULL : soluce = femFullSystemEliminate((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : soluce = femBandSystemEliminate((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : soluce = femIterativeSolverEliminate((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(soluce);
}

int femSolverConverged(femSolver *mySolver) // TOUS
{
    int  testConvergence;
    switch (mySolver->type) {
        case FEM_FULL : testConvergence = 1; break;
        case FEM_BAND : testConvergence = 1; break;
        case FEM_ITER : testConvergence = femIterativeSolverConverged((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(testConvergence);
}

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions femFullSystem
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

femFullSystem *femFullSystemCreate(int size) // FULL
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

void femFullSystemFree(femFullSystem *theSystem) // FULL
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}

void femFullSystemInit(femFullSystem *mySystem) // FULL
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++)
        mySystem->B[i] = 0;}

void femFullSystemPrint(femFullSystem *mySystem) // FULL
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

void femFullSystemPrintInfos(femFullSystem *mySystem) // FULL
{
    int  size = mySystem->size;
    printf(" \n");
    printf("    Full Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(size+1));
}

double femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol) // FULL
{
    return(myFullSystem->A[myRow][myCol]);
}

void  femFullSystemConstrain(femFullSystem *mySystem, int myNode, double myValue) // FULL
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

void femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc) // FULL
{
    int i,j;
    for (i = 0; i < nLoc; i++) {
        for(j = 0; j < nLoc; j++) {
            mySystem->A[map[i]][map[j]] += Aloc[i*nLoc+j]; }
    mySystem->B[map[i]] += Bloc[i]; }
}

double* femFullSystemEliminate(femFullSystem *mySystem) // FULL
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

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions femBandSystem
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

femBandSystem *femBandSystemCreate(int size, int band) // BANDE
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

void femBandSystemFree(femBandSystem *myBandSystem) // BANDE
{
    free(myBandSystem->B);
    free(myBandSystem->A);
    free(myBandSystem);
}

void femBandSystemInit(femBandSystem *myBandSystem) // BANDE
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++)
        myBandSystem->B[i] = 0;
}

void femBandSystemPrint(femBandSystem *myBand) // BANDE
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

void femBandSystemPrintInfos(femBandSystem *myBand) // BANDE
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

double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol) // BANDE
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol];
    return(value);
}

void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc) // BANDE
{
    int i,j;
    for (i = 0; i < nLoc; i++) {
        int myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            if (myCol >= myRow)  myBandSystem->A[myRow][myCol] += Aloc[i*nLoc+j]; }
        myBandSystem->B[myRow] += Bloc[i]; }
}

void femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue) // BANDE
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

double  *femBandSystemEliminate(femBandSystem *myBand) // BANDE
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

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions femIterativeSolver
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

femIterativeSolver *femIterativeSolverCreate(int size) // ITER
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

void femIterativeSolverFree(femIterativeSolver *mySolver) // ITER
{
    free(mySolver->R);
    free(mySolver);
}

void femIterativeSolverInit(femIterativeSolver *mySolver) // ITER
{
    int i;
    mySolver->iter = 0;
    mySolver->error = 10.0e+12;
    for (i=0 ; i < mySolver->size*4 ; i++)
        mySolver->R[i] = 0;
}

void femIterativeSolverPrint(femIterativeSolver *mySolver) // ITER
{
    double  *R;
    int     i, size;
    R    = mySolver->R;
    size = mySolver->size;

    for (i=0; i < size; i++) {
        printf("%d :  %+.1e \n",i,R[i]); }
}

void femIterativeSolverPrintInfos(femIterativeSolver *mySolver) // ITER
{
    if (mySolver->iter == 1)     printf("\n    Iterative solver \n");
    printf("    Iteration %4d : %14.7e\n",mySolver->iter,mySolver->error);
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc) // ITER
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

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue) // ITER
{
    mySolver->R[myNode] = myValue;
    mySolver->D0[myNode] = myValue;
    mySolver->D[myNode] = myValue;
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver) // ITER
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

int femIterativeSolverConverged(femIterativeSolver *mySolver) // ITER
{
    int  testConvergence = 0;
    if (mySolver->iter  > 3000)     testConvergence = -1;
    if (mySolver->error < 10.0e-6)  testConvergence = 1;
    return(testConvergence);
}

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions femTsunami
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

femTsunamiProblem *femTsunamiCreate(const char *meshFileName, int order)
{
    femTsunamiProblem *myProblem = malloc(sizeof(femTsunamiProblem));

    int i;
    int Rayon = 6371220.0;

    myProblem->gravity = 9.81;
    myProblem->gamma = 1e-7;
    myProblem->Omega = (2.0 * M_PI)/86400.0;
    myProblem->R = Rayon;

    myProblem->mesh = femMeshRead(meshFileName);

    myProblem->edges = femEdgesCreate(myProblem->mesh);

    if ((myProblem->mesh->nLocalNode == 4) && (order == 1))
    {
        myProblem->space = femDiscreteCreate(4,FEM_QUAD);
        myProblem->rule2d = femIntegrationCreate(9,FEM_QUAD);
        myProblem->nLocal = myProblem->space->n;
        myProblem->nNode = myProblem->mesh->nNode; // Ordre 2 : nNode + nEdge + nElem;

        myProblem->rule1d = femIntegrationCreate(2,FEM_EDGE);
    }
    else if ((myProblem->mesh->nLocalNode == 3) && (order == 1))
    {
        myProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        myProblem->rule2d = femIntegrationCreate(3,FEM_TRIANGLE);
        myProblem->nLocal = myProblem->space->n;
        myProblem->nNode = myProblem->mesh->nNode; // Ordre 2 : nNode + nEdge;

        myProblem->rule1d = femIntegrationCreate(2,FEM_EDGE);
    }
    int sizeLoc = myProblem->nLocal;
    int nElem = myProblem->mesh->nElem;
    int sizeGlo = nElem * sizeLoc + 1;
    //myProblem->solver = femSolverFullCreate(3*sizeLoc);
    myProblem->solver = femSolverBandCreate(3*sizeLoc,sizeLoc+1);// femBandSystemCreate

    myProblem->size = sizeGlo;
    myProblem->E = malloc(sizeof(double)*sizeGlo);
    myProblem->U = malloc(sizeof(double)*sizeGlo);
    myProblem->V = malloc(sizeof(double)*sizeGlo);
    myProblem->FE = malloc(sizeof(double)*sizeGlo);
    myProblem->FU = malloc(sizeof(double)*sizeGlo);
    myProblem->FV = malloc(sizeof(double)*sizeGlo);

    myProblem->elem = malloc(sizeof(double) * sizeLoc * nElem);
    myProblem->X = malloc(sizeof(double)*myProblem->nNode);
    myProblem->Y = malloc(sizeof(double)*myProblem->nNode);
    myProblem->bath = malloc(sizeof(double)*myProblem->nNode);

    femTsunamiMeshOrder(myProblem);

    for (i=0; i < myProblem->size; i++)
    {
        if ( i == ((myProblem->size)-1) )
        {
            myProblem->E[i] = 0.0;
            myProblem->U[i] = 0.0;
            myProblem->V[i] = 0.0;
        }
        else
        {
            int node = myProblem->elem[i];// ATTENTION AUX ORDRES SUPERIEURS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            myProblem->E[i] = tsunamiInitialConditionOkada(myProblem->X[node], myProblem->Y[node]);
            myProblem->U[i] = 0.0;
            myProblem->V[i] = 0.0;
        }
    }
    return myProblem;
}

void femTsunamiFree(femTsunamiProblem *myProblem)
{
    free(myProblem->bath);
    free(myProblem->Y);
    free(myProblem->X);
    free(myProblem->elem);
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

void femTsunamiMeshOrder(femTsunamiProblem *myProblem)
{
    int i, j;
    int order = myProblem->space->order;
    int nLocal = myProblem->nLocal;
    int nLocalNode = myProblem->mesh->nLocalNode;
    int nNodeMesh = myProblem->mesh->nNode;
    int nNodeTot = myProblem->nNode;
    int nElem = myProblem->mesh->nElem;

    //for (i=0; i < nNodeMesh; i++)
    //{
    //    myProblem->X[i] = myProblem->mesh->X[i];
    //    myProblem->Y[i] = myProblem->mesh->Y[i];
    //    myProblem->bath[i] = myProblem->mesh->bath[i];
    //}
    //for (i=0; i < nElem; i++)
    //{
    //    for (j=0; j < nLocalNode; j++)
    //    {
    //        myProblem->elem[i*nLocal+j] = myProblem->mesh->elem[i*nLocalNode+j];
    //    }
    //}
    int nmax;
    if (nNodeMesh>=nElem) { nmax = nNodeMesh; }
    else { nmax = nElem; }

    for (i=0; i < nmax; i++)
    {
        if (i<nNodeMesh)
        {
            myProblem->X[i] = myProblem->mesh->X[i];
            myProblem->Y[i] = myProblem->mesh->Y[i];
            myProblem->bath[i] = myProblem->mesh->bath[i];
        }
        if (i<nElem)
        {
            for (j=0; j < nLocalNode; j++)
            {
                myProblem->elem[i*nLocal+j] = myProblem->mesh->elem[i*nLocalNode+j];
            }
        }

    }
    //Ajout ordre 2 !

}

void femTsunamiTriangleMap(femTsunamiProblem *myProblem, int index, int *map)
{
    int nLocal = myProblem->nLocal;

    int j;
    for (j=0; j < nLocal; ++j)
    {
        map[j] = index * nLocal + j;
    }
}

void femTsunamiEdgeMap(femTsunamiProblem *myProblem, int index, int map[2][2])
{
    int nLN = myProblem->mesh->nLocalNode; // Ordre 2 modifier!!!!!!!!!!!!!!!

    int *nodeElem = myProblem->elem;
    femEdge edge = myProblem->edges->edges[index];
    int *node = edge.node;
    int *elem = edge.elem;

    int i;
    for( i=0 ; i<nLN ; i++ )
    {
        int j = i+1;
        if ( i==(nLN-1) ) { j=0; }
        int n1E1 = nodeElem[elem[0]*nLN+i];
        int n2E1 = nodeElem[elem[0]*nLN+j];
        if ( (n1E1==node[0]) && (n2E1==node[1]) )
        {
            map[0][0] = elem[0]*nLN+i;
            map[0][1] = elem[0]*nLN+j;
        }
        if (elem[1] != -1)
        {
            int n1E2 = nodeElem[elem[1]*nLN+i];
            int n2E2 = nodeElem[elem[1]*nLN+j];
            if ( (n1E2==node[1]) && (n2E2==node[0]) )
            {
                map[1][0] = elem[1]*nLN+j;
                map[1][1] = elem[1]*nLN+i;
            }
        }
    }

    if (elem[1]==(-1))
    {
        map[1][1] = myProblem->size-1;
        map[1][0] = myProblem->size-1;
    }
}

void femTsunamiAddIntegralsElements(femTsunamiProblem *myProblem)
{
    /* Infos Tsunami */
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
    int *elem = myProblem->elem;
    double *X = myProblem->X;
    double *Y = myProblem->Y;
    double *Bath = myProblem->bath;
    int nLocal = myProblem->nLocal;
    double  g     = myProblem->gravity;
    double  gamma = myProblem->gamma;
    double  R     = myProblem->R;
    double  Omega = myProblem->Omega;

    /* Infos Mesh */
    femMesh *mesh = myProblem->mesh;
    int nElem = mesh->nElem;
    int nLN = mesh->nLocalNode;

    /* Infos Space */
    int nSpace = myProblem->space->n;

    /* Infos Rule2D */
    int nRule = myProblem->rule2d->n;
    const double *xsi = myProblem->rule2d->xsi;
    const double *eta = myProblem->rule2d->eta;
    const double *w = myProblem->rule2d->weight;

    /* Allocations */
    double *phi= malloc(nSpace * sizeof(double));
    double *dphidxsi = malloc(nSpace * sizeof(double));
    double *dphideta = malloc(nSpace * sizeof(double));
    double *dphidx = malloc(nSpace * sizeof(double));
    double *dphidy = malloc(nSpace * sizeof(double));
    int *map = malloc(nLocal * sizeof(int));

    int iElem, i, iRule, k;
    /* Parcours des éléments */
    for( iElem=0 ; iElem<nElem ; iElem ++ )
    {
        femTsunamiTriangleMap(myProblem, iElem, map);

        /* Boucle d'intégration */
        for( iRule=0 ; iRule<nRule ; iRule++ )
        {
            myProblem->space->dphi2dx(xsi[iRule], eta[iRule], dphidxsi, dphideta);
            myProblem->space->phi2(xsi[iRule], eta[iRule], phi);

            double dXdXsi = 0.0, dXdEta = 0.0, dYdXsi = 0.0, dYdEta = 0.0;
            double e = 0.0, u = 0.0, v = 0.0, x = 0.0, y = 0.0, h = 0.0;
            /* Calcul du vecteur dXdXsi et des interpolations e, u, v, h, x, y */
            for( k=0 ; k<nSpace ; k++ )
            {
                dXdXsi += X[elem[map[k]]] * dphidxsi[k];
                dXdEta += X[elem[map[k]]] * dphideta[k];
                dYdXsi += Y[elem[map[k]]] * dphidxsi[k];
                dYdEta += Y[elem[map[k]]] * dphideta[k];

                e += E[map[k]] * phi[k];
                u += U[map[k]] * phi[k];
                v += V[map[k]] * phi[k];
                h += Bath[elem[map[k]]] * phi[k];
                x += X[elem[map[k]]] * phi[k];
                y += Y[elem[map[k]]] * phi[k];
            }

            double jac = fabs(dXdXsi * dYdEta - dXdEta * dYdXsi);
            double dXsidX = dYdEta/jac;
            double dXsidY = -dXdEta/jac;
            double dEtadX = -dYdXsi/jac;
            double dEtadY = dXdXsi/jac;
            /* Calcul des dphidx et dphidy */
            for( k=0 ; k<nSpace ; k++ )
            {
                dphidx[k] = dphidxsi[k] * dXsidX + dphideta[k] * dEtadX;
                dphidy[k] = dphidxsi[k] * dXsidY + dphideta[k] * dEtadY;
            }

            /* Calcul du coefficient */
            double coef = (4.0 * R * R + x * x + y * y)/(4.0 * R * R);
            /* Calcul du facteur f du terme de Coriolis */
            double z3d = R*(4.0*R*R - x*x - y*y) / (4.0*R*R + x*x + y*y);
            double lat = asin(z3d/R);
            double f = 2.0 * Omega * sin(lat);

            /* Parcours des fonctions PHI */
            for( i=0 ; i<nLocal ; i++ )
            {
                /* Calculs des valeurs des degrées de liberté (ajout de l'intégrale de surface) */
                BE[map[i]] += (dphidx[i] * u + dphidy[i] * v) * h * jac * w[iRule] * coef;
                BE[map[i]] += phi[i] * ((h * (x * u + y * v))/(R * R)) * jac * w[iRule];

                BU[map[i]] += (phi[i] * (f * v - gamma * u) + dphidx[i] * g * e * coef) * jac * w[iRule];
                BU[map[i]] += phi[i] * ((g * x * e)/(2.0 * R * R)) * jac * w[iRule];

                BV[map[i]] += (phi[i] * (-f * u - gamma * v) + dphidy[i] * g * e * coef) * jac * w[iRule];
                BV[map[i]] += phi[i] * ((g * y * e)/(2.0 * R * R)) * jac * w[iRule];
            }
        }
    }
    /* Désallocations */
    free(phi);
    free(dphidxsi);
    free(dphideta);
    free(dphidx);
    free(dphidy);
    free(map);
}

void femTsunamiAddIntegralsEdges(femTsunamiProblem *myProblem)
{
    /* Infos Tsunami */
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
    double *X = myProblem->X;
    double *Y = myProblem->Y;
    double *Bath = myProblem->bath;
    double  g = myProblem->gravity;
    double  R = myProblem->R;
    int nNodeEdge = myProblem->space->order + 1;

    /* Infos Rule1D */
    int nRule = myProblem->rule1d->n;
    const double *xsi = myProblem->rule1d->xsi;
    const double *w = myProblem->rule1d->weight;

    /* Infos Edges */
    femEdges *edgesStruct = myProblem->edges;
    int nEdge = edgesStruct->nEdge;
    int mapEdge[2][2] = {{0,0},{0,0}};

    /* Allocation */
    double *phiEdge = malloc(nNodeEdge * sizeof(double));

    int iEdge, iRule, i;
    /* Parcours des segments */
    for ( iEdge=0 ; iEdge<nEdge ; iEdge++ )
    {
        femTsunamiEdgeMap(myProblem, iEdge, mapEdge);
        femEdge myEdge = edgesStruct->edges[iEdge];

        /* Information segment courant */
        double coordX[2] = { X[myEdge.node[0]] , X[myEdge.node[1]] }; // A MODIFIER POUR L'ORDRE 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        double coordY[2] = { Y[myEdge.node[0]] , Y[myEdge.node[1]] };
        double bathy[2] = { Bath[myEdge.node[0]] , Bath[myEdge.node[1]] };

        /* Calcul de la norme et du jacobien */
        double dx = coordX[1] - coordX[0];
        double dy = coordY[1] - coordY[0];
        double norme = sqrt(dx*dx+dy*dy);
        double Jac = norme/2.0;


        /* Boucle d'intégration */
        for ( iRule=0 ; iRule<nRule; iRule++)
        {
            myProblem->space->phi1(xsi[iRule], phiEdge);

            /* Calcul norme */
            double nx = dy/norme;
            double ny = -dx/norme;

            /* Calcul des interpolations */
            double eta_L = 0.0, u_L = 0.0,
                   eta_R = 0.0, u_R = 0.0,
                   h = 0.0, x = 0.0, y = 0.0;


            h = bathy[0] * phiEdge[0] + bathy[1] * phiEdge[1]; // A MODIFIER POUR L'ORDRE 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            x = coordX[0] * phiEdge[0] + coordX[1] * phiEdge[1];
            y = coordY[0] * phiEdge[0] + coordY[1] * phiEdge[1];

            u_L = nx * femDiscreteInterpolate(phiEdge, U, mapEdge[0], 2) + ny * femDiscreteInterpolate(phiEdge, V, mapEdge[0], 2); // A MODIFIER POUR L'ORDRE 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            eta_L = femDiscreteInterpolate(phiEdge, E, mapEdge[0], 2);
            if (myEdge.elem[1]==-1)
            {
                eta_R = eta_L;
                u_R = - u_L;
            }
            else
            {
                u_R = nx * femDiscreteInterpolate(phiEdge, U, mapEdge[1], 2) + ny * femDiscreteInterpolate(phiEdge, V, mapEdge[1], 2); // A MODIFIER POUR L'ORDRE 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                eta_R = femDiscreteInterpolate(phiEdge, E, mapEdge[1], 2);
            }

            double eta_etoile = (eta_L + eta_R)/2.0 + sqrt(h/g) * ((u_L - u_R)/2.0);
            double u_etoile = (u_L + u_R)/2.0 + sqrt(g/h) * ((eta_L - eta_R)/2.0);

            /* Coefficient */
            double coef = (4.0 * R * R + x * x + y * y)/(4.0 * R * R);

            /* Parcours des fonctions PHI */
            for( i=0 ; i<nNodeEdge ; i++ )
            {
                /* Calculs des valeurs des degrées de liberté (ajout de l'intégrale de ligne) */
            	BE[mapEdge[0][i]] -= phiEdge[i] * h * u_etoile * coef * w[iRule] * Jac;
            	BE[mapEdge[1][i]] += phiEdge[i] * h * u_etoile * coef * w[iRule] * Jac;

                BU[mapEdge[0][i]] -= phiEdge[i] * nx * g * eta_etoile * coef * w[iRule] * Jac;
                BU[mapEdge[1][i]] += phiEdge[i] * nx * g * eta_etoile * coef * w[iRule] * Jac;

                BV[mapEdge[0][i]] -= phiEdge[i] * ny * g * eta_etoile * coef * w[iRule] * Jac;
                BV[mapEdge[1][i]] += phiEdge[i] * ny * g * eta_etoile * coef * w[iRule] * Jac;
            }
        }
    }

    /* Désallocation */
    free(phiEdge);
}

void femTsunamiMultiplyInverseMatrix(femTsunamiProblem *myProblem)
{
    /* Infos Tsunami */
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *X = myProblem->X;
    double *Y = myProblem->Y;
    int *elem = myProblem->elem;
    int nLocal = myProblem->nLocal;

    /* Infos Solver */
    femSolver *theSolver = myProblem->solver;

    /* Infos Mesh */
    femMesh *theMesh = myProblem->mesh;
    int nElem = theMesh->nElem;
    int nLN = theMesh->nLocalNode;

    /* Infos Space */
    femDiscrete *theSpace = myProblem->space;
    int nSpace = theSpace->n;

    /* Infos Rule2D */
    femIntegration *theRule = myProblem->rule2d;
    int nRule = theRule->n;
    const double *xsi = theRule->xsi;
    const double *eta = theRule->eta;
    const double *w = theRule->weight;

    /* Allocations */
    double *phi= malloc(nSpace * sizeof(double));
    double *dphidxsi = malloc(nSpace * sizeof(double));
    double *dphideta = malloc(nSpace* sizeof(double));
    int *map = malloc(nLocal * sizeof(int));
    int *mapE = malloc( nLocal * sizeof(int));
    int *mapU = malloc( nLocal * sizeof(int));
    int *mapV = malloc( nLocal * sizeof(int));
    double *Aloc = malloc( nLocal * nLocal * sizeof(double));

    int iElem, i, j, iRule;
    /* Calcul des vecteurs mapE, mapU et mapV */
    for( i=0 ; i<nLocal ; i++ )
    {
        mapE[i] = i;
        mapU[i] = i + nLN;
        mapV[i] = i + 2 * nLN;
    }
    /* Parcours des éléments */
    for( iElem=0 ; iElem<nElem ; iElem++ )
    {
        femTsunamiTriangleMap(myProblem, iElem, map);

        /* Mise à zéro de Aloc */
        for( i=0 ; i<(nLocal*nLocal) ; i++ ) { Aloc[i] = 0.0; }
        /* Boucle d'intégration pour le calcul de Aloc */
        for( iRule=0 ; iRule<nRule ; iRule++ )
        {
            theSpace->dphi2dx(xsi[iRule], eta[iRule], dphidxsi, dphideta);
            theSpace->phi2(xsi[iRule], eta[iRule], phi);

            /* Calcul du vecteur dXdXsi*/
            double dXdXsi = 0.0, dXdEta = 0.0, dYdXsi = 0.0, dYdEta = 0.0;
            for( i=0 ; i<nSpace ; i++ )
            {
                dXdXsi += X[elem[map[i]]] * dphidxsi[i];
                dXdEta += X[elem[map[i]]] * dphideta[i];
                dYdXsi += Y[elem[map[i]]] * dphidxsi[i];
                dYdEta += Y[elem[map[i]]] * dphideta[i];
            }

            double jac = fabs(dXdXsi * dYdEta - dXdEta * dYdXsi);

            /* Double boucle pour donner ses valeurs à Aloc */
            for( i=0 ; i<nSpace ; i++ )
            {
                for( j=i ; j<nSpace ; j++ )
                {
                    Aloc[i*nLocal+j] += jac * w[iRule] * phi[i] * phi[j];
                    if (i!=j) { Aloc[j*nLocal+i] = Aloc[i*nLocal+j]; }
                }
            }
        }

        /* Résolution du sytème avec le solver */
        femSolverInit(theSolver);
        femSolverAssemble(theSolver, Aloc, &BE[map[0]], NULL, mapE, nSpace);
        femSolverAssemble(theSolver, Aloc, &BU[map[0]], NULL, mapU, nSpace);
        femSolverAssemble(theSolver, Aloc, &BV[map[0]], NULL, mapV, nSpace);
        double *soluce = femSolverEliminate(theSolver);

        /* Placement de la solution dans les vecteurs BE, BU et BV */
        for (i=0; i<nLocal; i++)
        {
            BE[map[i]] = soluce[mapE[i]];
            BU[map[i]] = soluce[mapU[i]];
            BV[map[i]] = soluce[mapV[i]];
        }
    }

    /* Désallocations */
    free(map);
    free(mapE);
    free(mapU);
    free(mapV);
    free(Aloc);
    free(dphidxsi);
    free(dphideta);
    free(phi);
}

void femTsunamiCompute(femTsunamiProblem *theProblem)
{
    /* Infos Tsunami */
    int  size = theProblem->size;
    double *FE = theProblem->FE;
    double *FU = theProblem->FU;
    double *FV = theProblem->FV;
    double *E = theProblem->E;
    double *U = theProblem->U;
    double *V = theProblem->V;
    double theTimeStep = theProblem->timeStep;

    /* Allocations de vecteurs de stockage temporaire*/
    double *Etemp = malloc( size * sizeof(double));
    double *Utemp = malloc( size * sizeof(double));
    double *Vtemp =  malloc( size * sizeof(double));

    int i;
    /* Utemp = Un */
    for (i=0; i < size; i++)
    {
        Etemp[i] = E[i];
        Utemp[i] = U[i];
        Vtemp[i] = V[i];
        FE[i] = 0.0;
        FU[i] = 0.0;
        FV[i] = 0.0;
    }
    femTsunamiAddIntegralsElements(theProblem);
    femTsunamiAddIntegralsEdges(theProblem);
    femTsunamiMultiplyInverseMatrix(theProblem);
    /* Utemp = Un + dt/2 * FU */
    /* U = Un + dt * FU */
    for (i=0; i < size; i++)
    {
        Etemp[i] += (theTimeStep/2.0) * FE[i];
        Utemp[i] += (theTimeStep/2.0) * FU[i];
        Vtemp[i] += (theTimeStep/2.0) * FV[i];

        E[i] += theTimeStep * FE[i];
        U[i] += theTimeStep * FU[i];
        V[i] += theTimeStep * FV[i];

        FE[i] = 0.0;
        FU[i] = 0.0;
        FV[i] = 0.0;
    }
    femTsunamiAddIntegralsElements(theProblem);
    femTsunamiAddIntegralsEdges(theProblem);
    femTsunamiMultiplyInverseMatrix(theProblem);
    /* Un+1 = Un + dt/2 * FU + dt/2 * F(Un + dt * FU) */
    for (i=0; i < size; i++)
    {
        E[i] = Etemp[i] + (theTimeStep/2.0) * FE[i];
        U[i] = Utemp[i] + (theTimeStep/2.0) * FU[i];
        V[i] = Vtemp[i] + (theTimeStep/2.0) * FV[i];
    }

    /* Désallocation */
    free(Etemp);
    free(Utemp);
    free(Vtemp);
}

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions d'erreur, ...
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

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

/*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*
* Fonctions principales
*
* &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

void tsunamiCompute(double dt, int nmax, int sub, int order, const char *meshFileName, const char *baseResultName)
{
    femTsunamiProblem *theProblem = femTsunamiCreate(meshFileName, order);
    theProblem->timeStep = dt;
    double theDiscreteTime = 0.0;
    double theStop = nmax * dt;
    int theIteration = 0;


    int nElem = theProblem->mesh->nElem;
    int nLN = theProblem->mesh->nLocalNode;

    while (theStop>=theDiscreteTime)
    {
        if (theIteration % sub == 0)
        {
            double *U = theProblem->U;
            double *V = theProblem->V;
            double *E = theProblem->E;
            printf("Iteration : %d\n",theIteration);
            tsunamiWriteFile(baseResultName, theIteration, U, V, E, nElem, nLN);
        }
        femTsunamiCompute(theProblem);
        theIteration += 1;
        theDiscreteTime += dt;
    }

    femTsunamiFree(theProblem);
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
    	frame = (int) ((t-t0) * 2); // PROBLEME!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(frame*nout<=nmax) //MODIFICATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {
            if (frame0 != frame)
            {

                glfwGetMousePos( &mouse, NULL );  mouse = 389;

                //char filename[256];
                //const char *basename = "tsunami-%08d.txt"; // "tsunami-%08d.txt";
                //sprintf(filename, basename, frame * nout);
                const char *basename = "%s-%08d.txt";
                char filename[256];
                sprintf(filename,basename,baseResultName,frame * nout);
        // A decommenter lorsque vous avez des fichiers de resultats
                //if (!access(filename, F_OK))
                //{
                //    glfwTerminate();
                //    exit( EXIT_SUCCESS );
                //}

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

                for (i=0; i < nElem; ++i)
                {
                    for (j=0; j < 3; ++j)
                    {
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
                        coord[j*3+2] = (4*R*R - x*x - y*y)*R / Factor;
                    }

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
                    glDisableClientState(GL_VERTEX_ARRAY);
                }
                glfwSwapBuffers();

            }
        }
        else //MODIFICATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        {
            glfwTerminate();
            exit( EXIT_SUCCESS );
        }
    }
    while( glfwGetKey( GLFW_KEY_ESC ) != GLFW_PRESS && glfwGetWindowParam( GLFW_OPENED ) );

    glfwTerminate();
    exit( EXIT_SUCCESS );
}

double tsunamiInitialConditionOkada(double x, double y)
{
    double R = 6371220;
    double x3d = 4*R*R*x / (4*R*R + x*x + y*y);
    double y3d = 4*R*R*y / (4*R*R + x*x + y*y);
    double z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
    double lat = asin(z3d/R)*180/M_PI;
    double lon = atan2(y3d,x3d)*180/M_PI;
    double lonMin = 142;
    double lonMax = 143.75;
    double latMin = 35.9;
    double latMax = 39.5;
    double olon = (lonMin+lonMax)/2;
    double olat = (latMin+latMax)/2;
    double angle = -12.95*M_PI/180;
    double lon2 = olon + (lon-olon)*cos(angle) + (lat-olat)*sin(angle);
    double lat2 = olat - (lon-olon)*sin(angle) + (lat-olat)*cos(angle);
    if ( lon2 <= lonMax && lon2 >= lonMin &&
         lat2 >= latMin && lat2 <= latMax )
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

void tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub)
{
    int i,j;
    const char *basename = "%s-%08d.txt";
    char filename[256];
    sprintf(filename,basename,baseResultName,iter);
    FILE* file = fopen(filename,"w");
    fprintf(file, "Number of elem %d \n", nelem);
    fprintf(file, "Number of local values per element %d \n", nsub);
    for (i = 0; i < nelem; ++i)
    {
    	for (j = 0; j < nsub; ++j)
        {
        	int index = i*nsub+j;
        	fprintf(file,"%d;%d;%le;%le;%le;\n",i,j,U[index],V[index],E[index]);
        }
    }
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
