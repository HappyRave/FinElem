# include "fem.h"

double interpolate(double *phi, double *Vector, int *map, int n);

void femAdvInitialCondition(femAdvProblem *myProblem)
{
    double  u,v,elevation;
    int     j,elem,*node;    
    double *X = myProblem->mesh->X;
    double *Y = myProblem->mesh->Y;
    double *C = myProblem->C;
    double *U = myProblem->U;
    double *V = myProblem->V;
 
    for (elem=0; elem < myProblem->mesh->nElem; elem++) {
        node = &(myProblem->mesh->elem[elem*3]);
        for (j=0; j < 3; ++j) {
        	double x = X[node[j]];
        	double y = Y[node[j]]; 
            double c = exp(- pow((x - 0.8), 2) / 0.01) * exp(- pow((y - 0.5), 2) /0.01);
            femStommel(x, y, &u, &v, &elevation);
			C[elem*3+j] = c; 
            U[elem*3+j] = u; 
            V[elem*3+j] = v; }}
    int last = myProblem->size - 1;
    C[last] = 0.0; U[last] = 0.0; V[last] = 0.0;
}

# ifndef NOINTEGRALSTRIANGLES

void femAdvAddIntegralsTriangles(femAdvProblem *myProblem)
{
	double phi[3];
	double dphix[3];
	double dphiy[3];
	
	int i;
	for (i=0; i<myProblem->mesh->nElem; i++) {
		int map[3];
		femAdvTriangleMap(myProblem,i,map);
		int *vertices = &(myProblem->mesh->elem[i*3]);
		double x[3] = {myProblem->mesh->X[vertices[0]],myProblem->mesh->X[vertices[1]],myProblem->mesh->X[vertices[2]]};
		double y[3] = {myProblem->mesh->Y[vertices[0]],myProblem->mesh->Y[vertices[1]],myProblem->mesh->Y[vertices[2]]};
		
		dphix[0]=(y[1]-y[2]);
		dphix[1]=(y[2]-y[0]);
		dphix[2]=(y[0]-y[1]);
		
		dphiy[0]=(x[2]-x[1]);
		dphiy[1]=(x[0]-x[2]);
		dphiy[2]=(x[1]-x[0]);
		
		int k;
		for(k=0;k<3;k++) {
			
			int j;
			for (j=0; j<myProblem->rule2d->n; j++) {
				
				double xsi = myProblem->rule2d->xsi[j];
				double eta = myProblem->rule2d->eta[j];
				double weight = myProblem->rule2d->weight[j];
				
				myProblem->space->phi2(xsi, eta, phi);
				
				double u=interpolate(phi,myProblem->U, map, 3);
				double v=interpolate(phi,myProblem->V, map, 3);
				double c=interpolate(phi,myProblem->C, map, 3);
			
				myProblem->F[map[k]]+=weight*(u*dphix[k]+v*dphiy[k])*c;
			}
		}
	}
}

# endif
# ifndef NOINTEGRALSEDGES

void femAdvAddIntegralsEdges(femAdvProblem *myProblem)
{
	int n = myProblem->rule1d->n;
	double *xsi = (double*)myProblem->rule1d->xsi;
	double *weight = (double*)myProblem->rule1d->weight;
	int map[2][2];

	int i;
	for (i=0; i<myProblem->edges->nEdge; i++) {
		
		femAdvEdgeMap(myProblem,i,map);
		
		int nodeOne = myProblem->edges->edges[i].node[0];
		int nodeTwo = myProblem->edges->edges[i].node[1];

		double x1 = myProblem->mesh->X[nodeOne];
		double x2 = myProblem->mesh->X[nodeTwo];
		double y1 = myProblem->mesh->Y[nodeOne];
		double y2 = myProblem->mesh->Y[nodeTwo];
		double len = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
		
		double nx = -(y2-y1)/len;
		double ny = (x2-x1)/len;
	
		int j;
		for (j=0; j<n; j++) {
			double phi[n];
			myProblem->space->phi1(xsi[j], phi);
			double u = (interpolate(phi,myProblem->U,map[0],n)+interpolate(phi,myProblem->U,map[1],n))/2.0;
			double v = (interpolate(phi,myProblem->V,map[0],n)+interpolate(phi,myProblem->V,map[1],n))/2.0;
			double beta = u*nx+v*ny;
			
			double c;
			if (beta <=0){
				c = interpolate(phi,myProblem->C,map[0],n);
			} else {
				c = interpolate(phi,myProblem->C,map[1],n);
			}
			
			myProblem->F[map[0][0]]+=weight[j]*phi[0]*beta*c*len/2.0;
			myProblem->F[map[0][1]]+=weight[j]*phi[1]*beta*c*len/2.0;
			myProblem->F[map[1][0]]-=weight[j]*phi[0]*beta*c*len/2.0;
			myProblem->F[map[1][1]]-=weight[j]*phi[1]*beta*c*len/2.0;
		}
    }
}

# endif
# ifndef NOINVERSE

void femAdvMultiplyInverseMatrix(femAdvProblem *myProblem)
{

	double invA[3][3] = {{18.0,-6.0,-6.0},{-6.0,18.0,-6.0},{-6.0,-6.0,18.0}};
	int map[3];
	double Je;
	
	int el;
	for(el=0; el<myProblem->mesh->nElem; el++){
		
		int *vertices = &(myProblem->mesh->elem[el*3]);
		double x[3] = {myProblem->mesh->X[vertices[0]],myProblem->mesh->X[vertices[1]],myProblem->mesh->X[vertices[2]]};
		double y[3] = {myProblem->mesh->Y[vertices[0]],myProblem->mesh->Y[vertices[1]],myProblem->mesh->Y[vertices[2]]};

		Je=(x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
		Je = fabs(Je);

		femAdvTriangleMap(myProblem,el,map);

		double Bej[3] = {myProblem->F[map[0]], myProblem->F[map[1]], myProblem->F[map[2]]};

		int i;
		for(i=0; i<3; i++){
			myProblem->F[map[i]] = 0.0;

			int j;
			for(j=0; j<3; j++){
				myProblem->F[map[i]] += invA[i][j]*Bej[j]/Je;
			}	

		}
	}
}

# endif
# ifndef TRIANGLEMAP

void femAdvTriangleMap(femAdvProblem *myProblem, int index, int map[3])
{

    // Fourni gracieusement en l'honneur de l'election papale de Francois 1er :-)
 
    int j;
    for (j=0; j < 3; ++j) 
        map[j] = index*3 + j; 
}

# endif
# ifndef EDGEMAP

void femAdvEdgeMap(femAdvProblem *myProblem, int index, int map[2][2])
{
 
	int triLeft = myProblem->edges->edges[index].elem[0];
	int triRight = myProblem->edges->edges[index].elem[1];
	int nodeOne = myProblem->edges->edges[index].node[0];
	int nodeTwo = myProblem->edges->edges[index].node[1];
	
	int ext = 0;
	if (triRight < 0) {
		ext++;
		triRight = 0; //dummy value
	}
	
	int mapLeft[3];
	int mapRight[3];
	femAdvTriangleMap(myProblem,triLeft,mapLeft);
	femAdvTriangleMap(myProblem,triRight,mapRight);
	
	int leftNodes[3];
	int *verticesLeft = &(myProblem->mesh->elem[triLeft*3]);
	int rightNodes[3];
	int *verticesRight = &(myProblem->mesh->elem[triRight*3]);
	int i;
	for (i=0; i<3; i++) {
		leftNodes[i] = verticesLeft[i];
		rightNodes[i] = verticesRight[i];
	}
	
	int j;
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			if (leftNodes[i] == nodeOne && rightNodes[j] == leftNodes[i]) {
				map[0][0] = mapLeft[i];
				map[1][0] = mapRight[j];
			}
			if (leftNodes[i] == nodeTwo && rightNodes[j] == leftNodes[i]) {
				map[0][1] = mapLeft[i];
				map[1][1] = mapRight[j];
			}
			if (ext == 1 && leftNodes[i] == nodeTwo) {
				map[0][1] = mapLeft[i];
			}
		}
	}
	
	if (ext == 1) {
		map[1][0] = myProblem->mesh->nElem*3;
		map[1][1] = myProblem->mesh->nElem*3;
	}
}

# endif

double interpolate(double *phi, double *Vector, int *map, int n) {
	double u = 0.0;
	int i;
	for (i=0; i < n ; i++) {
		u += phi[i]*Vector[map[i]];
	}
	return u;
}






