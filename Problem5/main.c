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



int main(void)
{  
    printf("\n\n    V : Plot results \n");
    printf("    S : Spy matrix \n");
    printf("    F-B-I : Full solver - Band solver - Ierative solver \n");
    printf("    X-Y-N : Renumbering along x - along y - No renumbering \n");
    
    femSolverType solverType = FEM_BAND;
    femRenumType  renumType  = FEM_YNUM;
    char meshFileName[] = "triangles_6315.txt";  
    femDiffusionProblem* theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
    clock_t tic = clock();
    femDiffusionCompute(theProblem);  
    femSolverPrintInfos(theProblem->solver);
    printf("    CPU time : %.9f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    fflush(stdout);
    
    glfemInit("MECA1120 : homework 5 ");

    int option = 1;    
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;

    do 
    {
        int testConvergence,w,h;
        char theMessage[256];
        sprintf(theMessage, "Max : %.4f ",femMax(theProblem->soluce,theProblem->size));
        glfwGetWindowSize(&w,&h);
        
        if (option == 1) {
            glfemReshapeWindows(theProblem->mesh,w,h);
            glfemPlotField(theProblem->mesh,theProblem->soluce);   }
        else {
            glColor3f(1.0,0.0,0.0);
            glfemPlotSolver(theProblem->solver,theProblem->size,w,h); }
        glColor3f(0.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
        glfwSwapBuffers();
   
        if (solverType != newSolverType || renumType != newRenumType) { 
            solverType = newSolverType;
            renumType = newRenumType;
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
            clock_t tic = clock();
            do {
                femDiffusionCompute(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.9f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            switch (renumType) {
                case FEM_XNUM : printf("    Renumbering along the x-direction\n"); break;
                case FEM_YNUM : printf("    Renumbering along the y-direction\n"); break;
                default : break; }
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout); }
            
        if (glfwGetKey('V') == GLFW_PRESS)   option = 1;
        if (glfwGetKey('S') == GLFW_PRESS)   option = 0;
        if (glfwGetKey('F') == GLFW_PRESS)   newSolverType = FEM_FULL; 
        if (glfwGetKey('B') == GLFW_PRESS)   newSolverType = FEM_BAND; 
        if (glfwGetKey('I') == GLFW_PRESS)   newSolverType = FEM_ITER; 
        if (glfwGetKey('X') == GLFW_PRESS)   newRenumType  = FEM_XNUM; 
        if (glfwGetKey('Y') == GLFW_PRESS)   newRenumType  = FEM_YNUM; 
        if (glfwGetKey('N') == GLFW_PRESS)   newRenumType  = FEM_NO; 
        
  
    } 
    while( glfwGetKey(GLFW_KEY_ESC) != GLFW_PRESS &&
           glfwGetWindowParam(GLFW_OPENED) );
           
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);
    
    

}

