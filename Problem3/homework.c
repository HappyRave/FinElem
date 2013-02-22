//
//  Created by Laurent Debersaques and Maxime De Mol on 14/02/13.
//  Copyright (c) 2013 DaKot. All rights reserved.
//

#include "fem.h"

int			edgeCompare(const void* e0, const void *e1);
void		edgeExpand(femEdges *theEdges);
void		edgeSort(femEdges *theEdges);
void		edgeShrink(femEdges *theEdges);
double		edgeBoundaryLength(femEdges *theEdges);
void 		removeAndCompress(femEdge* allEdges, int i, int* n);

void edgeExpand(femEdges *theEdges)
{
	femMesh *mesh=theEdges->mesh;
	int nLocalNode=mesh->nLocalNode;
	int nElem=mesh->nElem;
	int *elem=mesh->elem;
	femEdge* edges=theEdges->edges;
	
	int quadNr;
	int edgeNr;
	for (quadNr=0; quadNr<nElem; quadNr++) {
		
		femEdge *segmentID = &edges[quadNr*nLocalNode];
		int *verticesID = &elem[quadNr*nLocalNode];
		
		for (edgeNr=0; edgeNr<nLocalNode; edgeNr++) {
			segmentID[edgeNr]=(femEdge){{quadNr,-1}, {verticesID[edgeNr], verticesID[(edgeNr+1)%nLocalNode]}};
		}
	}
}

void edgeSort(femEdges *theEdges)
{
	qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), edgeCompare);
}

int edgeCompare(const void* e0, const void *e1)
{
	int node0min;
	int node0max;
	int node1min;
	int node1max;
	
	if (((femEdge*) e0)->node[0]<((femEdge*) e0)->node[1]) {
		node0min=((femEdge*) e0)->node[0];
		node0max=((femEdge*) e0)->node[1];
	} else {
		node0min=((femEdge*) e0)->node[1];
		node0max=((femEdge*) e0)->node[0];
	}
	if (((femEdge*) e1)->node[0]<((femEdge*) e1)->node[1]) {
		node1min=((femEdge*) e1)->node[0];
		node1max=((femEdge*) e1)->node[1];
	} else {
		node1min=((femEdge*) e1)->node[1];
		node1max=((femEdge*) e1)->node[0];
	}
	
	if (node0min>node1min) {
		return(-1);
	} else if (node0min<node1min) {
		return(1);
	} else {
		if (node0max>node1max) {
			return(-1);
		} else if (node0max<node1max) {
			return(1);
		} else {
			return 0;
		}
	}
}

void edgeShrink(femEdges *theEdges)
{
	
	/* --------- Legatus ---------  */
	
    // Ici commence votre contribution a faire :-)
	
    int n = 0;          // Nouveau nombre total de segments : A MODIFIER
    int nBoundary = 0;  // Nombre de segments frontieres : A MODIFIER
	
	/* ------- End Legatus -------  */
	
	
	// n = 0 ? What a joke !
	n = theEdges->nEdge;
	
	// Shorcut variable
	femEdge* allEdges = theEdges->edges;
	
	// Iteration variables
	int i = 0;
	femEdge* currentFemEdge = NULL;
	femEdge* nextFemEdge = NULL;
	
	while(i+1 < n){
		
		// Current edge in the array
		currentFemEdge = &(allEdges[i]);
		// Getting the next edge for comparison with the current one
		nextFemEdge = &(allEdges[i+1]);
		
		// Checking if the two edges are actually equal
		if(edgeCompare(currentFemEdge,nextFemEdge) == 0){
			currentFemEdge->elem[1] = nextFemEdge->elem[0];
			removeAndCompress(allEdges,i+1,&n);
		} else {
			nBoundary++;
		}
		
		i++;
	}
	
	// Treating the last edge of the array, not considered in the previous loop
	nBoundary++;
	
	
	/* --------- Legatus ---------  */
	
    // Ici, finit votre contribution
    
    // Reallocation du tableau des edges
    
    theEdges->edges = realloc(theEdges->edges, n * sizeof(femEdge));
    theEdges->nEdge = n;
    theEdges->nBoundary = nBoundary;
	
	/* ------- End Legatus -------  */
	
}

double edgeBoundaryLength(femEdges *theEdges)
{
	
	// Init data
	double R = 6371000;
	double toReturn = 0;
	
	// Shorcut variable
	femEdge* allEdges = theEdges->edges;
	double* X = theEdges->mesh->X;
	double* Y = theEdges->mesh->Y;
	int n = theEdges->nEdge;
	
	// Loop on all edges
	int i = 0;
	femEdge currentFemEdge;
	for(i = 0; i < n; i++){
		
		currentFemEdge = allEdges[i];
		
		if(currentFemEdge.elem[1] == -1){
			
			double Long1 = X[currentFemEdge.node[0]]*M_PI/180;
			double Lat1 = Y[currentFemEdge.node[0]]*M_PI/180;
			double Long2 = X[currentFemEdge.node[1]]*M_PI/180;
			double Lat2 = Y[currentFemEdge.node[1]]*M_PI/180;
			
			double X1 = R*cos(Long1)*cos(Lat1);
			double Y1 = R*sin(Long1)*cos(Lat1);
			double X2 = R*cos(Long2)*cos(Lat2);
			double Y2 = R*sin(Long2)*cos(Lat2);
			double Z1 = R*sin(Lat1);
			double Z2 = R*sin(Lat2);
			
			toReturn += sqrt((Y2-Y1)*(Y2-Y1) + (X2-X1)*(X2-X1) + (Z2-Z1)*(Z2-Z1));
		}
		
	}
	
    return toReturn;
}

/**
 * removeAndCompress(femEdge** allEdges, int i, int* n)
 * Removes entry in the array (*allEdges) at position i by
 * decrementing by 1 all the entries at position > i.
 * The last element is then set to a femEdge considered empty
 */
void removeAndCompress(femEdge* allEdges, int i, int* n){
	
	int k;
	for (k=i; k<(*n)-1; k++){
		allEdges[k] = allEdges[k+1];
	}
	allEdges[k]=(femEdge){{-1,-1}, {-1, -1}};
	(*n)--;
}