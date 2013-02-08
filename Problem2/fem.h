
/*
 *  fem.h
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2012 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

# define Error(a)   femError(a,  __LINE__, __FILE__)

typedef enum {FEM_TRIANGLE,FEM_QUAD} femElementType;

typedef struct {
    int *elem;
    double *X;
    double *Y;
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;

typedef struct {
    float **elem;
    float Ox;
    float Oy;
    float dx;
    int nx;
    int ny;
} femGrid;

typedef struct {
  	int n;
  	const double *xsi;
  	const double *eta;
  	const double *weight;
} femIntegration;

femIntegration 	    *femIntegrationCreate(int n, femElementType type);
femMesh             *femMeshQuadRead(const char *filename);
femGrid			    *femGridRead(const char *filename);
void				 femIntegrationFree(femIntegration *theRule);
void				 femMeshFree(femMesh *theMesh);
void				 femGridFree(femGrid *theGrid);

void        		 femMeshQuadWrite(const femMesh* myMesh, const char *filename);
double      		 femMin(double *x, int n);
double      		 femMax(double *x, int n);
void        	     femError(char *text, int line, char *file);

#endif