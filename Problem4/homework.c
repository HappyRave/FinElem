#include"fem.h"

// Utils (for program)
void computeAe(double** Ae, double* x, double* y, femDiscrete* space, femIntegration* rule);
void gradientProducts(double** gP, double xsi, double eta, double* x, double* y, femDiscrete* space);
double jacobian(double* x, double* y, double xsi, double eta, femDiscrete* space);
double** fForGauss(double** gP, double xsi, double eta, double* x, double* y, femDiscrete* space);

// Utils (for matrix handling)
double** mallocMatrix(int rows, int coll);
void freeMatrix(double** A, int rA);
void fillMatrixWithZeros(double** M, int rows, int cols);
void matrixTranspose(double** result, double** A, int rA, int cA);
void matrixProduct(double** result, double** A, int rA, int cA, double** B, int cB);
void multiplyAllInMatrix(double** A,int rA, int cA, double v);
void sumMatrices(double** result, double** A, double** B, int rows, int cols);

// Compute Be
void computeBLocal(femDiffusionProblem* theProblem, double* B,double* x,double* y);

femDiffusionProblem *diffusionCreate(const char *filename)
{
    femDiffusionProblem *theProblem = malloc(sizeof(femDiffusionProblem));
    theProblem->mesh  = femMeshRead(filename);
    theProblem->edges = femEdgesCreate(theProblem->mesh);
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}

void diffusionFree(femDiffusionProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}


void diffusionMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y) {
	int nbrNodes = theMesh->nLocalNode;
	int *elem = &(theMesh->elem[nbrNodes*i]);
	int j;
	for (j=0; j<nbrNodes; j++) {
		*(map+j)=*(elem+j);
		*(x+j)=theMesh->X[*(map+j)];
		*(y+j)=theMesh->Y[*(map+j)];
	}
}



void diffusionSolve(femDiffusionProblem *theProblem) {
	
	// Shortcut variables
	femMesh* mesh = theProblem->mesh;
	femDiscrete* space = theProblem->space;
    femIntegration* rule = theProblem->rule;
    femFullSystem* system = theProblem->system;
	double** A = system->A;
	double* B = system->B;
	int n = space->n;
	
	// Data structures
	int* map = malloc(sizeof(int)*n);
	double* x = malloc(sizeof(double)*n);
	double* y = malloc(sizeof(double)*n);
		
	int i = 0;
	for(i=0; i < mesh->nElem; i++){
		
		diffusionMeshLocal(mesh,i,map,x,y);
		
		// Computing Ae
		double** Ae = mallocMatrix(n,n);
		computeAe(Ae,x,y,space,rule);
		int j = 0;
		for(j = 0; j<n; j++){
			int k = 0;
			for(k = 0; k<n; k++){
				A[map[j]][map[k]]+=Ae[j][k];
			}
		}
		freeMatrix(Ae,n);
		
		// Computing Be
		double *Be=(double*)malloc(sizeof(double)*n);
		computeBLocal(theProblem, Be, x,y);
		j=0;
		for (j=0; j<n; j++) {
			B[map[j]]+=Be[j];
		}
	}
	
	//
	//femFullSystemPrint(system);
	
	free(map);
	free(x);
	free(y);
	
	femEdges* theEdges = theProblem->edges;
	femEdge* edge = theEdges->edges;
	i=0;
	for (i=0; i<theEdges->nEdge; i++) {
		if (edge[i].elem[1]==-1) {
		femFullSystemConstrain(system,edge[i].node[0],0.0);
		femFullSystemConstrain(system,edge[i].node[1],0.0);
		}
	}
	//
	//femFullSystemPrint(system);
	
	femFullSystemEliminate(system);
	
	//
	//femFullSystemPrint(system);
	
}
/* UTILS (for program) */

void computeAe(double** Ae, double* x, double* y, femDiscrete* space, femIntegration* rule){
	
	int n = space->n;
	
	if(n==3){
		gradientProducts(Ae,0,0,x,y,space);
		multiplyAllInMatrix(Ae,n,n,jacobian(x,y,0,0,space)/2);
	}
	else {
		
		fillMatrixWithZeros(Ae,n,n);
		
		// Integration rule
	  	int nRule = rule->n;
	  	const double *xsik = rule->xsi;
	  	const double *etak = rule->eta;
	  	const double *weightk = rule->weight;
		
		double** gP = mallocMatrix(n,n);
		
		// Gauss integration
		int i = 0; double xsi = 0; double eta = 0; double weight = 1;
		for(i = 0; i< nRule ; i++){
			xsi = xsik[i]; eta = etak[i]; weight = weightk[i];
			gradientProducts(gP,xsi,eta,x,y,space);
			multiplyAllInMatrix(gP,n,n,weight*jacobian(x, y, xsi, eta, space));
			sumMatrices(Ae, Ae, gP, n, n);
		}
		
		freeMatrix(gP,n);
		
	}
	
}

void gradientProducts(double** gP, double xsi, double eta, double* x, double* y, femDiscrete* space){
	
	// Computing the matrix dphi_i/dXsi
	double** dphidxiM = mallocMatrix(space->n, 2);
	double *dphidxsi = malloc(sizeof(double)* space->n);
	double *dphideta = malloc(sizeof(double)* space->n);
	space->dphi2dx(xsi, eta, dphidxsi, dphideta);
	int i = 0;
	for(i = 0; i < space->n; i++){
		dphidxiM[i][0] = dphidxsi[i];
		dphidxiM[i][1] = dphideta[i];
	}
	
	// Computing the matrix dXsi/dX
	double** dXsidXM = mallocMatrix(2,2);
	dXsidXM[0][0] = 0; dXsidXM[1][0] = 0;
	dXsidXM[0][1] = 0; dXsidXM[1][1] = 0;
	for(i = 0; i < space->n; i++){
		dXsidXM[0][0] += y[i]*dphideta[i];
		dXsidXM[0][1] -= x[i]*dphideta[i];
		dXsidXM[1][0] -= y[i]*dphidxsi[i];
		dXsidXM[1][1] += x[i]*dphidxsi[i];
	}
	double Je = (dXsidXM[0][0]*dXsidXM[1][1] - dXsidXM[0][1]*dXsidXM[1][0]);
	multiplyAllInMatrix(dXsidXM,2,2,1/Je);
	
	// Computing dphidX
	double** dphidXM = mallocMatrix(space->n,2);
	matrixProduct(dphidXM, dphidxiM, space->n, 2, dXsidXM, 2);
	
	// Free useless matrices
	freeMatrix(dphidxiM,space->n);
	freeMatrix(dXsidXM,2);
	
	// Compute the transpose of dphidX
	double** transposedMatrix = mallocMatrix(2,space->n);
	matrixTranspose(transposedMatrix,dphidXM,space->n,2);
	
	// Finally computing the gradient product
	matrixProduct(gP, dphidXM, space->n, 2, transposedMatrix, space->n);
	
	// Free useless matrices
	freeMatrix(dphidXM,space->n);
	freeMatrix(transposedMatrix,2);
	
}

// Computes the jacobian for any space type
double jacobian(double* x, double* y, double xsi, double eta, femDiscrete* space){
	
	double jacobian;
	double j[2][2] = {{0,0},{0,0}};
	
	double *dphidxsi = malloc(sizeof(double)* space->n);
	double *dphideta = malloc(sizeof(double)* space->n);
	space->dphi2dx(eta,xsi,dphidxsi, dphideta);
	int i = 0;
	for(i = 0; i < space->n; i++){
		j[0][0] += x[i]*dphidxsi[i];
		j[0][1] += x[i]*dphideta[i];
		j[1][0] += y[i]*dphidxsi[i];
		j[1][1] += y[i]*dphideta[i];
	}
	jacobian = (j[0][0]*j[1][1] - j[0][1]*j[1][0]);
	
	return fabs(jacobian);
}

/* UTILS (for matrices handling) */

// Dynamic allocation of a matrix and returns pointer
double** mallocMatrix(int rows, int coll){
	double** A = malloc(sizeof(double*)*rows);
	int i = 0;
	for(i=0; i<rows; i++){
		A[i] = malloc(sizeof(double)*coll);
	}
	return A;
}

// Free memory of a dynamically allocated matrix
void freeMatrix(double** A, int rA){
	int i=0;
	for(i=0; i<rA; i++){
		free(A[i]);
	}
	free(A);
}

// Fills the matix M with zeros
void fillMatrixWithZeros(double** M, int rows, int cols){
	int i = 0; int j = 0;
	for(i = 0; i < rows ; i++){
		for(j = 0; j < cols ; j++){
			M[i][j] = 0;
		}
	}
}

// Dot product between matrices
void matrixProduct(double** result, double** A, int rA, int cA, double** B, int cB){
	int i = 0;
	int j = 0;
	int k = 0;
	
	for(i = 0; i < rA ; i++){
		for(j = 0; j < cB ; j++){
			result[i][j] = 0;
			for(k = 0; k < cA; k++){
				result[i][j] += A[i][k]*B[k][j];
			}
		}
	}
}

// Computes the transpose of a matrix
void matrixTranspose(double** result, double** A, int rA, int cA){
	int i = 0;
	int j = 0;
	for(i = 0; i<rA ; i++){
		for(j=0; j<cA; j++){
			result[j][i] = A[i][j];
		}
	}
}

// Multiplies all elements in the matrix A by the constant v
void multiplyAllInMatrix(double** A,int rA, int cA, double v){
	int i = 0;
	int j = 0;
	for(i = 0; i < rA ; i++){
		for(j = 0; j < cA ; j++){
			A[i][j] *= v;
		}
	}
}


// Sums matrices A and B and stores result in matrix result
void sumMatrices(double** result, double** A, double** B, int rows, int cols){
	int i = 0; int j = 0;
	for(i = 0; i < rows ; i++){
		for(j = 0; j < cols ; j++){
			result[i][j] = A[i][j]+B[i][j];
		}
	}
}


void computeBLocal(femDiffusionProblem* theProblem, double* B,double* x,double* y)
{
	femIntegration* theRule = theProblem->rule;
	femDiscrete* theSpace = theProblem->space;
	femMesh *theMesh = theProblem->mesh;
	int n = theMesh->nLocalNode;
	int i;
	int j;
	
	double* phi = (double*)malloc(sizeof(double)*n);
	
	double* xsi = (double*)theRule->xsi;
	double* eta = (double*)theRule->eta;
	double* w = (double*)theRule->weight;
	int nPoints = theRule->n;
	
	for (i=0; i<nPoints; i++) {
		femDiscretePhi2(theSpace,xsi[i],eta[i],phi);
		for (j=0; j<n; j++) {
			B[j]+=w[i]*phi[j]*jacobian(x, y, xsi[i], eta[i], theSpace);
		}
	}
	free(phi);
}