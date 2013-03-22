Identification réussie : bienvenue à Maxime :-)
Inscription au cours Meca1120 enregistrée :-)
Votre numéro de groupe pour le cours Meca1120 est : 10 :-)
Devoir 5 du groupe 10 (soumise par Laurent Debersaques )


Voici la soumission 1 du devoir 5 du groupe 10 (soumise par Laurent Debersaques )

#include"fem.h"

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef struct {
	int num;
	double x;
	double y;
} node;

int compareX( const void* e0, const void* e1);
int compareY( const void* e0, const void* e1);

void femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType)
{
    int i;
	//int* elem = theProblem->mesh->elem;
	double* X = theProblem->mesh->X;
	double* Y = theProblem->mesh->Y;
	node tmp;
	node* nodes = (node*)malloc(sizeof(node)*theProblem->mesh->nNode);
	
	for (i = 0; i < theProblem->mesh->nNode; i++) {
		tmp.num = i; tmp.x = X[i]; tmp.y = Y[i];
		nodes[i] = tmp;
	}
	
	
    switch (renumType) {
        case FEM_NO :
      		for (i = 0; i < theProblem->mesh->nNode; i++)
        		theProblem->number[i] = i;
            break;
			//
			// A modifier :-)
			// debut
			//
        case FEM_XNUM :
			qsort(nodes,theProblem->mesh->nNode,sizeof(node),compareX);
			for (i = 0; i < theProblem->mesh->nNode; i++)
        		theProblem->number[nodes[i].num] = i;
			break;
        case FEM_YNUM :
            qsort(nodes,theProblem->mesh->nNode,sizeof(node),compareY);
			for (i = 0; i < theProblem->mesh->nNode; i++)
        		theProblem->number[nodes[i].num] = i;
            break;
			//
			// end
			//
			
        default	: Error("Unexpected renumbering option"); }
	free(nodes);
}

int compareX( const void* e0, const void* e1) {
	node* n0 = ((node*)e0);
	node* n1 = ((node*)e1);
	
	double x0 = n0->x;
	double x1 = n1->x;
	
	if (x0 == x1) {
		return 0;
	} else if (x0 > x1) {
		return(1);
	} else {
		return(-1);
	}
}

int compareY( const void* e0, const void* e1) {
	node* n0 = ((node*)e0);
	node* n1 = ((node*)e1);
	
	double x0 = n0->y;
	double x1 = n1->y;
	
	if (x0 == x1) {
		return 0;
	} else if (x0 > x1) {
		return(1);
	} else {
		return(-1);
	}
}

int femDiffusionComputeBand(femDiffusionProblem *theProblem)
{
	int* elem = theProblem->mesh->elem;
	int nLocalNode = theProblem->mesh->nLocalNode;
    int nElem = theProblem->mesh->nElem;
	int i = 0;
	int localMax = 0;
	int globMax = 0;
	for (i=0; i<nElem; i++) {
		int iMax;
		int jMax;
		int* seg = &elem[i*nLocalNode];
		if (nLocalNode == 3) {
			iMax = MIN(theProblem->number[seg[0]], MIN(theProblem->number[seg[1]], theProblem->number[seg[2]]));
			jMax = MAX(theProblem->number[seg[0]], MAX(theProblem->number[seg[1]], theProblem->number[seg[2]]));
		} else {
			iMax = MIN(MIN(theProblem->number[seg[0]], theProblem->number[seg[1]]), MIN(theProblem->number[seg[2]], theProblem->number[seg[3]]));
			jMax = MAX(MAX(theProblem->number[seg[0]], theProblem->number[seg[1]]), MAX(theProblem->number[seg[2]], theProblem->number[seg[3]]));
		}
		localMax = abs(jMax - iMax)+1;
		if (localMax > globMax) globMax = localMax;
	}
    return(globMax);
}

// Question 3

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
	
	// d^0
	if(mySolver->iter == 0){
		int i,j;
		for (i = 0; i < nLoc; i++) {
			
			int myRow = map[i];
			
			// Computing R0 (i.e. : computing r^0)
			for(j = 0; j < nLoc; j++) {
			    mySolver->R0[myRow] += Aloc[i*nLoc+j]*Uloc[j];
			}
			mySolver->R0[myRow] -= Bloc[i];
			
			// Updating other data structures
			mySolver->R[myRow] = mySolver->R0[myRow];
			mySolver->D[myRow] = -mySolver->R0[myRow];
			mySolver->D0[myRow] = -mySolver->R0[myRow];
		}
		
		for (i = 0; i < nLoc; i++) {
			
			int myRow = map[i];
			
			// Computing S0 (i.e. : computing Ad^0)
			for(j = 0; j < nLoc; j++) {
			    mySolver->S0[myRow] += Aloc[i*nLoc+j]*mySolver->D0[map[j]];
			}
			
			// Updating other data structures
			mySolver->S[myRow] = mySolver->S0[myRow];
			
		}
		
		
		// d^k
	} else {
		
		int i,j;
		
		// Computing R (i.e., r^k) and S (i.e., Ad^k)
		for (i = 0; i < nLoc; i++) {
			int myRow = map[i];
			for(j = 0; j < nLoc; j++) {
			    mySolver->R[myRow] += Aloc[i*nLoc+j]*Uloc[j];
			}
			mySolver->R[myRow] -= Bloc[i];
		}
		
		// Computing D (i.e. : computing d^k)
		double Adk_Product[nLoc];
		double num = 0, deno = 0;
		for (i = 0; i < nLoc; i++) {
			int myRow = map[i];
			Adk_Product[i] = mySolver->S0[myRow];
			num += mySolver->R[myRow]*Adk_Product[i];
			deno += mySolver->R0[myRow]*Adk_Product[i];
		}
		
		double frac = num/deno;
		for(i = 0; i < nLoc; i++){
			int myRow = map[i];
			mySolver->D[myRow] += -mySolver->R[myRow] + frac*mySolver->D0[myRow];
		}
		
		for(i = 0; i < nLoc; i++){
			int myRow = map[i];
			for(j = 0; j < nLoc; j++) {
				mySolver->S[myRow] += Aloc[i*nLoc+j]*mySolver->D[map[j]];
			}
		}
		
	}
	
}

void femIterativeSolverConstrain(femIterativeSolver* mySolver, int myNode, double myValue)
{
	
	mySolver->D[myNode] = myValue;
	
}

double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
	
	// Descente de gradients conjuguÃ©s
	mySolver->iter++;
	double error = 0.0; int i;
	
	// Computing alpha_k
	double num = 0, deno = 0;
	for(i = 0; i < mySolver->size; i++){
		num -= mySolver->R[i] * mySolver->D[i];
		deno += mySolver->D[i] * mySolver->S[i];
	}
	double alpha_k = num/deno;
	
	// Computing error, multiply direction vector by alpha_k and updating datastructures
	for (i=0; i < mySolver->size; i++) {
	    error += (mySolver->R[i])*(mySolver->R[i]);
	    mySolver->D0[i] = mySolver->D[i]*alpha_k;
		mySolver->D[i] = 0;
		mySolver->S0[i] = mySolver->S[i]; mySolver->S[i] = 0;
		mySolver->R0[i] = mySolver->R[i]; mySolver->R[i] = 0;
	}
	mySolver->error = sqrt(error);
	return mySolver->D0;
	
}
