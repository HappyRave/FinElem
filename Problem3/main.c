/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2012 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"
#include "glfem.h"



int         edgeCompare(const void* e0, const void *e1);
void        edgeExpand(femEdges *theEdges);
void        edgeSort(femEdges *theEdges);
void        edgeShrink(femEdges *theEdges);
double      edgeBoundaryLength(femEdges *theEdges);

int main(void)
{   
    femMesh  *theMesh  = femMeshRead("triangles.txt");
    femEdges *theEdges = femEdgesCreate(theMesh);    
    
    edgeExpand(theEdges);                   //femEdgesPrint(theEdges);
    edgeSort(theEdges);                     //femEdgesPrint(theEdges);
    edgeShrink(theEdges);                   //femEdgesPrint(theEdges);
    printf("Boundary edges  : %i \n", theEdges->nBoundary);
    printf("Boundary length : %14.7e \n", edgeBoundaryLength(theEdges));

    char theMessage[256];
    sprintf(theMessage, "Boundary edges : %i", theEdges->nBoundary);
      
//
//  On superpose le maillage (en bleu), 
//  tous les segments frontieres (en noir),
//  et la frontiere (en rouge)
//
//  Au depart de votre travail, vous devriez obtenir un maillage bleu....
//  et a la fin de l'exercice un maillage noir avec bord rouge :-)
//

    glfemInit("MECA1120 : homework3 ");
    do
    {
        int w,h;
        glfwGetWindowSize(&w,&h);
        glfemReshapeWindows(theMesh,w,h);
        glColor3f(0.4,0.4,1.0); glfemPlotMesh(theMesh);
        glColor3f(0.0,0.0,0.0); glfemPlotEdges(theEdges);  
        glColor3f(1.0,0.0,0.0); glfemPlotBnd(theEdges);          
        glfemMessage(theMessage);              
        glfwSwapBuffers();
    } 
    while( glfwGetKey(GLFW_KEY_ESC) != GLFW_PRESS &&
           glfwGetWindowParam(GLFW_OPENED) );
           
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femEdgesFree(theEdges);
    femMeshFree(theMesh);
    exit(EXIT_SUCCESS);

    
      
    

}

