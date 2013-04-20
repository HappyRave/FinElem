# include "fem.h"

void femShallowMap(int index, int* map, int nLocal);    
void femShallowEdgeMap(femShallowProblem *myProblem, int index, int map[2][2]);

    
void femShallowMap(int index, int* map, int nLocal){
    int j;
    for (j=0; j < nLocal; ++j) 
        map[j] = index*nLocal + j; 
}

void femShallowEdgeMap(femShallowProblem *myProblem, int index, int map[2][2]){

    int i,j,k;
    
    for (j=0; j < 2; ++j) {
        int node = myProblem->edges->edges[index].node[j];
        for (k=0; k < 2; k++) {
            int elem = myProblem->edges->edges[index].elem[k];
            map[k][j] = (myProblem->mesh->nElem)*3;
            if (elem >= 0) {
                for (i=0; i < 3; i++) {
                    if (myProblem->mesh->elem[elem*3 + i] == node) {
                        map[k][j] = elem*3 + i;
					}
				}
			}
		}
	}

}


# ifndef NOINTEGRALSTRIANGLES

void femShallowAddIntegralsElements(femShallowProblem *myProblem)
{

	printf("IntEl\n");    
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
    double *X = myProblem->mesh->X;
    double *Y = myProblem->mesh->Y;
    
    double  h     = myProblem->height;
    double  g     = myProblem->gravity;
    double  L     = myProblem->L;
    double  gamma = myProblem->gamma;
    double  tau0  = myProblem->tau0;
    double  rho   = myProblem->rho;
    double  f0    = myProblem->f0;
    double  beta  = myProblem->beta;

    //
    // A completer
    //

    femMesh *mesh = myProblem->mesh;
	femIntegration *theRule = myProblem->rule2d;
	femDiscrete *theSpace = myProblem->space;

	// Number of nodes by element (3 for triangles, 4 for quads)
	int nLocal = myProblem->space->n; 
	// Getting the map of vertices composing the element i
    int map[nLocal];

	int el=0; int i = 0;
    for(el=0; el < mesh->nElem; el++){

        femShallowMap(el,map,nLocal);

		// Getting coordinates & values stored for each vertex
		double e[nLocal], u[nLocal], v[nLocal], x[nLocal], y[nLocal];

		int j = 0;
        for(j=0; j<nLocal; j++){
            x[j] = X[mesh->elem[map[j]]];
            y[j] = Y[mesh->elem[map[j]]];
            e[j] = E[map[j]];
            u[j] = U[map[j]];
            v[j] = V[map[j]];
        }

		// What is required for our integral
		double dphidx[nLocal];
		double dphidy[nLocal];
		double jacobian = 0;
		
		// In case of a triangle, nothing depends on (xsi,eta) and we can compute outside
		if(nLocal == 3){			
			if(el==0) printf("This is a triangular mesh ");
			jacobian = (x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0]);
		    dphidx[0] = (y[1] - y[2])/jacobian;
		    dphidx[1] = (y[2] - y[0])/jacobian;
		    dphidx[2] = (y[0] - y[1])/jacobian;
		    dphidy[0] = (x[2] - x[1])/jacobian;
		    dphidy[1] = (x[0] - x[2])/jacobian;
		    dphidy[2] = (x[1] - x[0])/jacobian;  
		}

		// Shall contain the values of the phis for a given (xsi,eta)
		double phi[nLocal];
		// Values given by the integration rule
		double xsi = 0, eta = 0, weight = 0;
		int k = 0;
		for (k=0; k < theRule->n; k++) {

            xsi = theRule->xsi[k];
            eta = theRule->eta[k];
            weight = theRule->weight[k];     

			
			if(nLocal == 4){
				
				// Getting values of dphi/dxsi
				double dphidxsi[4]; double dphideta[4];
				theSpace->dphi2dx(xsi, eta, dphidxsi, dphideta);

				// Computing elements of the matrix dx/dxsi + the jacobian
				double el1 = 0, el2 = 0, el3 = 0, el4 = 0;

				for(i = 0; i < nLocal; i++){
					el1 += x[i]*dphidxsi[i];
					el2 += x[i]*dphideta[i];
					el3 += y[i]*dphidxsi[i];
					el4 += y[i]*dphideta[i];
				}
				jacobian = el1*el4-el2*el3;

				// Computing dphidx = dphi/dxsi * dxsi/dx
				for(i=0; i< nLocal; i++){
					dphidx[i] = (dphidxsi[i]*el4 + dphideta[i]*-el3)/jacobian;
					dphidy[i] = (dphidxsi[i]*-el2 + dphideta[i]*el1)/jacobian;
				}

			}
			
			// Getting values of phis for this pair (xsi,eta)
            femDiscretePhi2(theSpace,xsi,eta,phi);

			// Interpolation of the values at (xsi,eta)
            double zeU = femDiscreteInterpolate(phi,U,map,nLocal);
            double zeV = femDiscreteInterpolate(phi,V,map,nLocal);
            //double zeE = femDiscreteInterpolate(phi,E,map,nLocal);

			// Integral computation
            for (i=0; i < nLocal; i++) {

				// Shortcut variables for integration
				double f = f0 + beta * (y[i] - L/2);
				double tau = tau0*sin( M_PI*(y[i]-L/2)/L );
				double eeh = 0; double ueh = 0; double veh = 0;
				int theBigVariableToIncrementMan = 0;
				for(theBigVariableToIncrementMan = 0; theBigVariableToIncrementMan < nLocal; theBigVariableToIncrementMan++){
					eeh += e[theBigVariableToIncrementMan]*phi[theBigVariableToIncrementMan];
					ueh += u[theBigVariableToIncrementMan]*phi[theBigVariableToIncrementMan];
					veh += v[theBigVariableToIncrementMan]*phi[theBigVariableToIncrementMan];
				}

				// Integration
                BE[map[i]] += (zeU*h*dphidx[i] + zeV*h*dphidy[i])*jacobian*weight; 
				BU[map[i]] += (phi[i]*(f*veh + tau/(rho*h) - gamma*ueh) + dphidx[i]*g*eeh)*jacobian*weight;
				BV[map[i]] += (phi[i]*(-f*ueh-gamma*veh) + dphidy[i]*g*eeh)*jacobian*weight;
			}
		}

	}

}

# endif
# ifndef NOINTEGRALSEDGES

void femShallowAddIntegralsEdges(femShallowProblem *myProblem)
{

	printf("IntEdges\n");
	myProblem = myProblem;
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
     
    double  h = myProblem->height;
    double  g = myProblem->gravity;

    //
    // A completer
    //

	femIntegration *theRule = myProblem->rule1d;
	femDiscrete *theSpace = myProblem->space;

	double  xEdge[2],yEdge[2],phiEdge[2];
    double  xsi,weight,jac;
    int     i,j,k,edge,mapEdge[2][2];

    for (edge=0; edge < myProblem->edges->nEdge; edge++) {
        
		// Getting the map of edges
		femShallowEdgeMap(myProblem,edge,mapEdge);

		// Getting coordinates of the edge's endpoints
        for (j=0; j < 2; ++j) {
            int node = myProblem->edges->edges[edge].node[j];
            xEdge[j] = myProblem->mesh->X[node];
            yEdge[j] = myProblem->mesh->Y[node]; 
		}

        double dxdxsi = xEdge[1] - xEdge[0];
        double dydxsi = yEdge[1] - yEdge[0];
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double normal[2] = {dydxsi/norm, -dxdxsi/norm};
        jac = norm / 2.0;

        for (k=0; k < theRule->n; k++) {

            xsi = theRule->xsi[k];
            weight = theRule->weight[k];     

            femDiscretePhi1(theSpace,xsi,phiEdge);

			double etaL = femDiscreteInterpolate(phiEdge,E,mapEdge[0],2);
			double etaR = femDiscreteInterpolate(phiEdge,E,mapEdge[1],2);

			double UL = femDiscreteInterpolate(phiEdge,U,mapEdge[0],2);
			double VL = femDiscreteInterpolate(phiEdge,V,mapEdge[0],2);
			double UnL = UL*normal[0]+ VL*normal[1];

			double UR = femDiscreteInterpolate(phiEdge,U,mapEdge[1],2);
			double VR = femDiscreteInterpolate(phiEdge,V,mapEdge[1],2);
			double UnR = UR*normal[0]+ VR*normal[1];

			double etaStar = (etaL+etaR)/2 + sqrt(h/g)*(UnL-UnR)/2;
			double UnStar = (UnL+UnR)/2 + sqrt(g/h)*(etaL-etaR)/2;

            for (i=0; i < 2; i++) {
                BE[mapEdge[0][i]] += (phiEdge[i]*h*UnStar)*jac*weight;
                BE[mapEdge[1][i]] += (phiEdge[i]*h*UnStar)*jac*weight;
                BU[mapEdge[0][i]] += (phiEdge[i]*normal[0]*g*etaStar)*jac*weight;
                BU[mapEdge[1][i]] += (phiEdge[i]*normal[0]*g*etaStar)*jac*weight;
                BV[mapEdge[0][i]] += (phiEdge[i]*normal[1]*g*etaStar)*jac*weight;
                BV[mapEdge[1][i]] += (phiEdge[i]*normal[1]*g*etaStar)*jac*weight;
			}

		}
	}
}

# endif
# ifndef NOMULTIPLYINVERSEMATRIX


void femShallowMultiplyInverseMatrix(femShallowProblem *myProblem){
	
	//printf("MultInv\n");
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    femMesh *theMesh = myProblem->mesh;
    femDiscrete *theSpace = myProblem->space;
    femSolver *theSolver = myProblem->solver;
    femIntegration *theRule = myProblem->rule2d;
    
    int n = theSpace->n;
    double Xloc[n],Yloc[n],Aloc[n*n];
    int iElem,i,j,mapElem[n],mapE[n],mapU[n],mapV[n];
    
    for (i = 0; i < n; i++){
        mapE[i] = i;
        mapU[i] = i + n;
        mapV[i] = i + 2*n;
	}
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femSolverInit(theSolver);
        for (i = 0; i < n*n; i++)  Aloc[i] = 0;
        int *mapCoord = &(myProblem->mesh->elem[iElem*n]);
        for (j=0; j < n; ++j) {
            mapElem[j] = iElem*n + j;
            Xloc[j] = myProblem->mesh->X[mapCoord[j]];
            Yloc[j] = myProblem->mesh->Y[mapCoord[j]]; }
		
		int iInteg;
		double phi[n],dphidxsi[n],dphideta[n];
		//double dphidx[n],dphidy[n];
		
		for (iInteg=0; iInteg < theRule->n; iInteg++) {
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
			int i;
			int j;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += Xloc[i]*dphidxsi[i];
                dxdeta += Xloc[i]*dphideta[i];
                dydxsi += Yloc[i]*dphidxsi[i];
                dydeta += Yloc[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
			/*
			 for (i = 0; i < theSpace->n; i++) {
			 dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
			 dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }
			 */
            for (i = 0; i < theSpace->n; i++) {
                for(j = 0; j < theSpace->n; j++) {
                    //Aloc[i+n*j] += (dphidx[i] * dphidx[j]+ dphidy[i] * dphidy[j]) * jac * weight;
					//printf("Aloc: %f\n",Aloc[i+n*j]);
					Aloc[i+n*j] += phi[i]*phi[j]*jac*weight;
				}
			}
		}
		
        femSolverAssemble(theSolver,Aloc,&BE[mapElem[0]],NULL,mapE,theSpace->n);
        femSolverAssemble(theSolver,Aloc,&BU[mapElem[0]],NULL,mapU,theSpace->n);
        femSolverAssemble(theSolver,Aloc,&BV[mapElem[0]],NULL,mapV,theSpace->n);
		//puts("Bitch");
		//femFullSystemPrint((femFullSystem*)theSolver);
		double *soluce = femSolverEliminate(theSolver);  // A decommenter
        for (i = 0; i < n; i++) {
			
			BE[mapElem[i]] = soluce[mapE[i]];
			BU[mapElem[i]] = soluce[mapU[i]];
			BV[mapElem[i]] = soluce[mapV[i]];
            
		}
	}
}


# endif


