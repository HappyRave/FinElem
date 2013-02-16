#include "fem.h"

int         edgeCompare(const void* e0, const void *e1);
void        edgeExpand(femEdges *theEdges);
void        edgeSort(femEdges *theEdges);
void        edgeShrink(femEdges *theEdges);
double      edgeBoundaryLength(femEdges *theEdges);


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
		
		femEdge *segmentID=&edges[quadNr*nLocalNode];
		int *verticesID=&elem[quadNr*nLocalNode];
		
		for (edgeNr=0; edgeNr<nLocalNode; edgeNr++) {
			if (edgeNr==nLocalNode-1) {
				segmentID[edgeNr]=(femEdge){{quadNr,-1}, {verticesID[edgeNr], verticesID[0]}};
			} else {
				segmentID[edgeNr]=(femEdge){{quadNr,-1}, {verticesID[edgeNr], verticesID[edgeNr+1]}};
			}
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
    // Ici commence votre contribution a faire :-)

    int n = 0;          // Nouveau nombre total de segments : A MODIFIER
    int nBoundary = 0;  // Nombre de segments frontieres : A MODIFIER
    
    // Ici, finit votre contribution
    
    // Reallocation du tableau des edges
    
    theEdges->edges = realloc(theEdges->edges, n * sizeof(femEdge));
    theEdges->nEdge = n;
    theEdges->nBoundary = nBoundary;
}

double edgeBoundaryLength(femEdges *theEdges)
{
    double L = 777;
    return L;
}

