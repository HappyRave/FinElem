/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2013 UCL-IMMC : Vincent Legat    
 *  All rights reserved.
 *
 */


#include "glfem.h"
#include <time.h>

femDiffusionProblem* theProblem;

double f(int iElem, double xsi, double eta, double x, double y)
{
    double phi[20],Uloc[20],Xloc[20],Yloc[20];
    int map[20];
    femDiscretePhi2(theProblem->space,xsi,eta,phi);
    femDiffusionMeshLocal(theProblem,iElem,map,Xloc,Yloc,Uloc);  

    double u = 0; 
    int i;
    for (i = 0; i < theProblem->space->n; i++)    
        u += Uloc[i]*phi[i];       
    return u;   
}


int main(void)
{  
    printf("\n\n    V : Plot results \n");
    printf("    S : Spy matrix \n");
    printf("    N : Show all nodes \n");
    printf("    F-B-I : Full solver - Band solver - Conjugate Gradients solver \n");
    printf("    X-Y-N : Renumbering along x - along y - No renumbering \n");
    printf("    1-2-3 : Order of local polynomial approximations \n");
    
    femSolverType solverType = FEM_BAND;
    femRenumType  renumType  = FEM_YNUM;
    int           orderType = 1;
    char meshFileName[] = "rect_quad_1601.txt";  
    theProblem = femDiffusionCreate(meshFileName,solverType,renumType,orderType);
    femDiffusionPrintInfos(theProblem);

    clock_t tic = clock();
    femDiffusionCompute(theProblem);  
    femSolverPrintInfos(theProblem->solver);
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    fflush(stdout);
    
    glfemInit("MECA1120 : homework 6 ");

    int option = 1;    
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;
    int           newOrderType  = orderType;

    do 
    {
        int testConvergence,w,h;
        char theMessage[256];
        double theMaxValue = femMax(theProblem->soluce,theProblem->size);
        sprintf(theMessage, "Max : %.4f ",theMaxValue);
        glfwGetWindowSize(&w,&h);
        
        if (option == 2) {
            glfemReshapeWindows(theProblem->mesh,w,h);
            glfemPlotFieldRecursive(theProblem->mesh,f,0.0,theMaxValue,orderType);  
            glColor3f(0.0,0.0,0.0);
            int node = theProblem->mesh->nNode;
            glfemDrawNodes(theProblem->X,theProblem->Y,node); 
            glColor3f(1.0,0.0,0.0);
            glfemDrawNodes(theProblem->X+node,theProblem->Y+node,theProblem->size-node); }
        else if (option == 1) {
            glfemReshapeWindows(theProblem->mesh,w,h);
            glfemPlotFieldRecursive(theProblem->mesh,f,0.0,theMaxValue,orderType);  }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h); }
            glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
        glfwSwapBuffers();
   
        if (solverType != newSolverType || renumType != newRenumType || orderType != newOrderType) { 
            solverType = newSolverType;
            renumType = newRenumType;
            orderType = newOrderType;
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType,orderType);
            femDiffusionPrintInfos(theProblem);
            clock_t tic = clock();
            do {
                femDiffusionCompute(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            switch (renumType) {
                case FEM_XNUM : printf("    Renumbering along the x-direction\n"); break;
                case FEM_YNUM : printf("    Renumbering along the y-direction\n"); break;
                default : break; }
                
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout); }
            
        if (glfwGetKey('N') == GLFW_PRESS)   option = 2;
        if (glfwGetKey('V') == GLFW_PRESS)   option = 1;
        if (glfwGetKey('S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey('F') == GLFW_PRESS)   newSolverType = FEM_FULL; 
        if (glfwGetKey('B') == GLFW_PRESS)   newSolverType = FEM_BAND; 
        if (glfwGetKey('I') == GLFW_PRESS)   newSolverType = FEM_ITER; 
        if (glfwGetKey('X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey('Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey('N') == GLFW_PRESS)   newRenumType  = FEM_NO; 
        if (glfwGetKey('1') == GLFW_PRESS)   newOrderType  = 1; 
        if (glfwGetKey('2') == GLFW_PRESS)   newOrderType  = 2; 
        if (glfwGetKey('3') == GLFW_PRESS)   newOrderType  = 3; 
        
    } 
    while( glfwGetKey(GLFW_KEY_ESC) != GLFW_PRESS &&
           glfwGetWindowParam(GLFW_OPENED) );
           
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);
    
    

}

