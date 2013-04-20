/*
 *  main.c
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2013 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

#include "glfem.h"
# define NORMAL_MODE 0
# define NOPLOT_MODE 1
# define PLOT_MODE   2

static femShallowProblem *theProblem;
static double *thePlotField;
int theProgramMode = NORMAL_MODE;

static double draw(int ielem, double xsi, double eta, double x, double y)
{
    int m = theProblem->space->n;
    double phi[m];
    femDiscretePhi2(theProblem->space,xsi,eta,phi);
    double c = 0.0; int i;
    for (i=0; i<m; i++)
        c += fabs(thePlotField[ielem*m+i])*phi[i];
    return c;
}



int main(int argc, char *argv[])
{   

//
//  Options possibles
//    -noplot : version non graphique avec creation de fichiers de resultats
//    -plot : lecture des fichiers et animation
//    (par défaut) : version graphique interactive
//

    while(argc--) {
        if (strcmp( *argv, "-noplot" ) == 0)    theProgramMode = NOPLOT_MODE;
        if (strcmp( *argv++, "-plot" ) == 0)    theProgramMode = PLOT_MODE; }
    
    theProblem = femShallowCreate("triangles_4254.txt");
	//theProblem = femShallowCreate("quads_313.txt");
  
    double theStartingTime = 0.0;
    double theDiscreteTime = 0.0;
    double theTime;
    int    theIteration = 0;
    
//
//  Pas de temps plus ou moins optimal pour le maillage triangles_4254.txt
//  Pour un maillage plus fin, il faut le diminuer....
//  Pour un maillage plus grossier, on peut l'augmenter.
//  La sélection automatique du pas de temps fera partie du projet final :-)
//    
    theProblem->timeStep = 0.0004;  // 0.0008 rigolo
    double theTimeStep = theProblem->timeStep; 
    double theStop = 2000*theTimeStep;
    char   filename[256], theMessage[256];
    
    if (theProgramMode == NOPLOT_MODE) {
        while (theStop >= theDiscreteTime) {
                theIteration += 1;
                theDiscreteTime += theTimeStep;          
                femShallowCompute(theProblem);
                if (theIteration % 100 == 0) {
                    printf("Iteration = %6d \n",theIteration);
                    sprintf(filename,"data/result-%d.txt",theIteration);
                    femShallowWrite(theProblem,filename);
                    fflush(stdout); }}}
    else {
        thePlotField = theProblem->U;
        glfemInit("MECA1120 : homework 8 ");
        do 
        {
            int w,h;
            glfwGetWindowSize(&w,&h); 
            if (theProgramMode != PLOT_MODE) 
                 theTime = (glfwGetTime() - theStartingTime) / 100;   // 5
            else theTime = (glfwGetTime() - theStartingTime) / 5;
                
       
    //
    //  Pour figer/ne pas figer le resultat a un temps, commenter/decommenter la ligne ci-dessous
    //  
            if (theTime >= theStop) theTime = theStop;
       
            while (theTime >= theDiscreteTime) 
            {
                if (theProgramMode != PLOT_MODE) {
                    theIteration += 1;
                    theDiscreteTime += theTimeStep;          
                    femShallowCompute(theProblem); }
                else { 
                    theIteration += 100;
                    theDiscreteTime += theTimeStep*100;
                    sprintf(filename,"data/result-%d.txt",theIteration);
                    femShallowRead(theProblem,filename);
                    fflush(stdout); }   
                
            }
            sprintf(theMessage,"Time = %.3f Iter = %d Step = %.4f ",theTime,theIteration-1,theTimeStep);

            glfemReshapeWindows(theProblem->mesh,w,h);
            double maxScale = femMax(thePlotField,theProblem->size);
            double minScale = femMin(thePlotField,theProblem->size);
            maxScale = fmax(fabs(maxScale),fabs(minScale));

            glfemPlotFieldRecursive(theProblem->mesh,draw,0.0,maxScale,1);  
            glColor3f(0.0, 0.0, 0.0);
            glfemPlotMesh(theProblem->mesh);
            glColor3f(1.0,0.0,0.0); glfemDrawMessage(20,460,theMessage);      
            glfwSwapBuffers();

            if (glfwGetKey('U') == GLFW_PRESS)   thePlotField = theProblem->U; 
            if (glfwGetKey('V') == GLFW_PRESS)   thePlotField = theProblem->V; 
            if (glfwGetKey('E') == GLFW_PRESS)   thePlotField = theProblem->E; 
            if (glfwGetKey('R') == GLFW_PRESS) {
                theStartingTime = glfwGetTime(); 
                theDiscreteTime = 0;
                theIteration = 0;
                }
              
                  
        } 
        while( glfwGetKey(GLFW_KEY_ESC) != GLFW_PRESS &&
               glfwGetWindowParam(GLFW_OPENED) );
      
        glfwTerminate(); 
    }
    
    femShallowFree(theProblem);    
    exit(EXIT_SUCCESS);

}

