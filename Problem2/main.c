/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2012 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#ifdef graphic

// ======== Version graphique (pour une compilation sur votre ordinateur)

#include "glfem.h"

double earthRadius();
double earthGridInterpolate(femGrid *theGrid, double x, double y);
double earthJacobian(double x[4], double y[4], double xsi, double eta, double R);
double earthIntegrateBathymetry(femMesh *theMesh, femGrid *theGrid, femIntegration *theRule, double R);

int main(void)
{   
    femMesh *theMesh = femMeshQuadRead("mesh.txt");
    femGrid *theGrid = femGridRead("bath.txt");
    femIntegration *theRule = femIntegrationCreate(4,FEM_QUAD);
    double  R = earthRadius();
    double  I = earthIntegrateBathymetry(theMesh,theGrid,theRule,R); 
    printf("Sea volume = %10.3e m^3",I);

    
    // ======= Graphics 
     
    char theMessage[256];
    double *theField = malloc(sizeof(double)*theMesh->nNode); int i;
    for (i=0; i < theMesh->nNode; ++i)
        theField[i]  = -fmin(0.0,earthGridInterpolate(theGrid,theMesh->X[i],theMesh->Y[i]));          
    sprintf(theMessage,"Sea volume = %10.3e m^3 \n",I);
    fflush(stdout);
        
    glfemInit("MECA1120 : homework2 ");
    do
    {
        int w,h;
        glfwGetWindowSize(&w,&h);
    	glfemReshapeWindows(theMesh,w,h);
    	glfemPlotField(theMesh,theField); 
        glfemMessage(theMessage);              
    	glfwSwapBuffers();
    } 
    while( glfwGetKey( GLFW_KEY_ESC ) != GLFW_PRESS &&
		   glfwGetWindowParam( GLFW_OPENED ) );
           
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femIntegrationFree(theRule);
    femMeshFree(theMesh);
    femGridFree(theGrid);   
    exit(EXIT_SUCCESS);
    

}

#else
 
// ======== Version non graphique (test du serveur)

#include "fem.h"

double earthRadius();
double earthGridInterpolate(femGrid *theGrid, double x, double y);
double earthJacobian(double x[4], double y[4], double xsi, double eta, double R);
double earthIntegrateBathymetry(femMesh *theMesh, femGrid *theGrid, femIntegration *theRule, double R);

int main(void)
{   
    femMesh *theMesh = femMeshQuadRead("mesh.txt");
    femGrid *theGrid = femGridRead("bath.txt");
    femIntegration *theRule = femIntegrationCreate(4,FEM_QUAD);
    double  R = earthRadius();
    double  I = earthIntegrateBathymetry(theMesh,theGrid,theRule,R); 
    printf("Sea volume = %10.3e m^3\n",I);
    
    double *theField = malloc(sizeof(double) * theMesh->nNode); int i;
    for (i = 0; i < theMesh->nNode; ++i)
        theField[i]  = -fmin(0.0,earthGridInterpolate(theGrid,theMesh->X[i],theMesh->Y[i])); 
    int index=45;
    printf("Bathymetry at node %d is : %10.3e m\n",index,theField[index]);  
    
    femIntegrationFree(theRule);
    femMeshFree(theMesh);
    femGridFree(theGrid);   
    exit(0);
    
}

#endif
 
