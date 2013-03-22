/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2013 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */


#include "glfem.h"

femDiffusionProblem *diffusionCreate(const char *filename);
void     diffusionFree(femDiffusionProblem *theProblem);
void     diffusionMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y);
void     diffusionSolve(femDiffusionProblem *theProblem);

int main(void)
{   
    femDiffusionProblem* theProblem = diffusionCreate("triangles_6315.txt");
    diffusionSolve(theProblem);   
 
    printf("Number of elements : %4d\n", theProblem->mesh->nElem);
    printf("Number of unknowns : %4d\n", theProblem->system->size);
    printf("Maximum value : %.4f\n", femMax(theProblem->system->B,theProblem->system->size));
    fflush(stdout);
    
    char theMessage[256];
    sprintf(theMessage, "Max : %.4f", femMax(theProblem->system->B,theProblem->system->size));
    
    glfemInit("MECA1120 : homework4 ");
    do
    {
        int w,h;
        glfwGetWindowSize(&w,&h);
        glfemReshapeWindows(theProblem->mesh,w,h);
        glfemPlotField(theProblem->mesh,theProblem->system->B);            
        glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);              
        glfwSwapBuffers();
    } 
    while( glfwGetKey(GLFW_KEY_ESC) != GLFW_PRESS &&
           glfwGetWindowParam(GLFW_OPENED) );
           
    // Check if the ESC key was pressed or the window was closed
               
    glfwTerminate(); 
    diffusionFree(theProblem);
    exit(EXIT_SUCCESS);
    
    

}

