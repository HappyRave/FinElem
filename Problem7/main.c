/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2012 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"

static femAdvProblem *theProblem;

static double draw(int ielem, double xsi, double eta, double x, double y)
{
    int m = theProblem->space->n;
 	double phi[m];
    femDiscretePhi2(theProblem->space,xsi,eta,phi);
	double c = 0.0; int i;
	for (i=0; i<m; i++)
		c += theProblem->C[ielem*m+i]*phi[i];
    return c;
}

int main(void)
{   
    theProblem = femAdvCreate("square.txt");
    femAdvInitialCondition(theProblem);
    
    int  size = theProblem->size, i;
    double *F = theProblem->F;  
    double *C = theProblem->C;    
    double theStartingTime = 0.0;
    double theDiscreteTime = 0.0;
    int    theIteration = 0;
    double theTimeStep = 0.04;
   
    
    glfemInit("MECA1120 : homework 7 ");
    do 
    {
        char theMessage[256];
   	    int w,h;
        glfwGetWindowSize(&w,&h); 
        double theTime = (glfwGetTime() - theStartingTime) * 5;   // 5
   
//
//  Pour figer/ne pas figer le resultat a un temps, commenter/decommenter les deux lignes ci-dessous
//  
        double theStop = 12000*theTimeStep;
        if (theTime >= theStop) theTime = theStop;
   
   		while (theTime >= theDiscreteTime) 
        {
        	theIteration += 1;
        	theDiscreteTime += theTimeStep;      	 
            for (i=0; i < size; i++) 
        		F[i] = 0.0;
            femAdvAddIntegralsTriangles(theProblem);
    		femAdvAddIntegralsEdges(theProblem);
    		femAdvMultiplyInverseMatrix(theProblem);
            for (i=0; i < size; i++) 
        		C[i] += theTimeStep * F[i];
        }
   		sprintf(theMessage,"Time = %.3f Iter = %d Step = %.3f ",theTime,theIteration-1,theTimeStep);

       	glfemReshapeWindows(theProblem->mesh,w,h);
 		glfemPlotFieldRecursive(theProblem->mesh,draw,0.0,1.0,1);	
		glColor3f(0.0, 0.0, 0.0);
    	glfemPlotMesh(theProblem->mesh);
     	glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);      
  	   	glfwSwapBuffers();

        if (glfwGetKey('R') == GLFW_PRESS) {
        	theStartingTime = glfwGetTime(); 
    		theDiscreteTime = 0;
    		theIteration = 0;
    		femAdvInitialCondition(theProblem); }
          
              
    } 
    while( glfwGetKey(GLFW_KEY_ESC) != GLFW_PRESS &&
           glfwGetWindowParam(GLFW_OPENED) );
  
    glfwTerminate(); 
    femAdvFree(theProblem);    
    exit(EXIT_SUCCESS);

}

