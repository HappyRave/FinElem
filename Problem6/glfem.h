/*
 *  glfem.h
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2013 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 *  Pour GLFW (version utilisée 2.7.7)
 *  Pour l'installation de la librairie, voir http://www.glfw.org/
 *
 */

#ifndef _GLFEM_H_
#define _GLFEM_H_

#include <GL/glfw.h>
#include "fem.h"

void glfemDrawColorElement(double *x, double *y, double *u, int n);
void glfemDrawColorRecursiveElement(double *x, double *y, 
                                     double *xsi, double *eta, 
                                     double (*f)(int,double,double,double,double),
                                     double fMin, double fMax,
                                     int iElem,
                                     int level, int nLocalNode);
void glfemDrawElement(double *x, double *y, int n);
void glfemDrawNodes(double* x, double* y,int n);

void glfemReshapeWindows(femMesh *theMesh, int width, int heigh);
void glfemPlotField(femMesh *theMesh, double *u);
void glfemPlotFieldRecursive(femMesh *theMesh, double (*f)(int,double,double,double,double),double fMin, double fMax,int level);
void glfemPlotMesh(femMesh *theMesh);
void glfemPlotEdges(femEdges *theEdges);
void glfemPlotBnd(femEdges *theEdges);

void glfemMatrix(double **A, int size, int width, int heigh);
void glfemPlotSolver(femSolver *theSolver, int size, int width, int heigh);

void glfemMessage(char *message);
void glfemDrawMessage(int h, int v, char *message);
void glfemSetRasterSize(int width, int height);
void glfemInit(char *windowName);

#endif