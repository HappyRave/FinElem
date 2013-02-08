/*
 *  fem.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2012 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

static const double _gaussQuad4Xsi[4]     = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]     = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4]  = { 1.0, 1.0, 1.0, 1.0};


femIntegration *femIntegrationCreate(int n, femElementType type)
{
	femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
  		theRule->n      = 4;
  		theRule->xsi    = _gaussQuad4Xsi;
  		theRule->eta    = _gaussQuad4Eta;
  		theRule->weight = _gaussQuad4Weight; }
    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
	free(theRule);
}

femMesh *femMeshQuadRead(const char *filename)
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
    
    fscanf(file, "Number of quads %d \n", &theMesh->nElem);  
    theMesh->elem = malloc(sizeof(int)*4*theMesh->nElem);
    theMesh->nLocalNode = 4;
    for (i = 0; i < theMesh->nElem; ++i) {
        elem = &(theMesh->elem[i*4]);
        fscanf(file,"%d : %d %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2],&elem[3]);   }
    
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

void femMeshQuadWrite(const femMesh *theMesh, const char *filename)
{
    int i,*elem;
    
    FILE* file = fopen(filename,"w");
    
    fprintf(file, "Number of nodes %d \n", theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fprintf(file,"%6d : %14.7e %14.7e \n",i,theMesh->X[i],theMesh->Y[i]); }
    
    fprintf(file, "Number of quads %d \n", theMesh->nElem);  
    for (i = 0; i < theMesh->nElem; ++i) {
        elem = &(theMesh->elem[i*4]);
        fprintf(file,"%6d : %6d %6d %6d %6d \n", i,elem[0],elem[1],elem[2],elem[3]);   }
    
    fclose(file);
}

femGrid *femGridRead(const char *filename)
{
    femGrid *theGrid = malloc(sizeof(femGrid));

    int i,j,nx,ny;
    float missingData;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No grid file !");
    
    fscanf(file, "NCOLS        %d \n",&theGrid->nx); nx = theGrid->nx;
    fscanf(file, "NROWS        %d \n",&theGrid->ny); ny = theGrid->ny;
    fscanf(file, "XLLCENTER    %f \n",&theGrid->Ox);
    fscanf(file, "YLLCENTER    %f \n",&theGrid->Oy);
    fscanf(file, "CELLSIZE     %f \n",&theGrid->dx);
    fscanf(file, "NODATA_VALUE %f \n",&missingData);
    
    theGrid->elem = malloc(sizeof(int*) * nx);
    theGrid->elem[0] = malloc(sizeof(int) * nx * ny);
    for (i=1 ; i < nx ; i++) 
        theGrid->elem[i] = theGrid->elem[i-1] + ny;
    
    for (j = 0; j < ny; j ++)
        for (i = 0; i < nx; i++)
            fscanf(file,"%f",&theGrid->elem[i][ny-j-1]);
    
    //
    // -!- La première valeur dans les données correspond au max lat et min lon.
    //      
        
    fclose(file);
    return theGrid;
}

void femGridFree(femGrid *theGrid)
{
	free(theGrid->elem[0]);
    free(theGrid->elem);
    free(theGrid);
}

double femMin(double *x, int n) {
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n) {
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
