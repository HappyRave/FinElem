# include "fem.h"


double interpolate(double *phi, double *U, int *map, int n)
{
    double u = 0.0; int i;
    for (i=0; i <n; i++)
        u += phi[i]*U[map[i]];
    return u;
}

void femShallowEdgeMap(femShallowProblem *myProblem, int index, int map[2][2], int n)
{
    int i,j,k;
    
    for (j=0; j < 2; ++j) {
        int node = myProblem->edges->edges[index].node[j];
        for (k=0; k < 2; k++) {
            int elem = myProblem->edges->edges[index].elem[k];
            map[k][j] = (myProblem->mesh->nElem)*n;
            if (elem >= 0) {
                for (i=0; i < n; i++) {
                    if (myProblem->mesh->elem[elem*n + i] == node) {
                        map[k][j] = elem*n + i;  }}}}}
}


# ifndef NOINTEGRALSTRIANGLES

void femShallowAddIntegralsElements(femShallowProblem *myProblem)
{
	
    
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
    double *Y = myProblem->mesh->Y;
    
    double  h     = myProblem->height;
    double  g     = myProblem->gravity;
    double  L     = myProblem->L;
    double  gamma = myProblem->gamma;
    double  tau0  = myProblem->tau0;
    double  rho   = myProblem->rho;
    double  f0    = myProblem->f0;
    double  beta  = myProblem->beta;
	
    femIntegration *theRule = myProblem->rule2d;
    femDiscrete *theSpace = myProblem->space;
	
	int n = theSpace->n;
    double  Xloc[n],Yloc[n];
    double  u,v,e,f,tau,y;
    int     j,iElem,mapElem[n];
    
    for (iElem=0; iElem < myProblem->mesh->nElem; iElem++) {
        int *mapCoord = &(myProblem->mesh->elem[iElem*n]);
        for (j=0; j < n; ++j) {
            mapElem[j] = iElem*n + j;
            Xloc[j] = myProblem->mesh->X[mapCoord[j]];
            Yloc[j] = myProblem->mesh->Y[mapCoord[j]]; }
		
        int iInteg;
		double phi[n],dphidxsi[n],dphideta[n],dphidx[n],dphidy[n];
		
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
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += Xloc[i]*dphidxsi[i];
                dxdeta += Xloc[i]*dphideta[i];
                dydxsi += Yloc[i]*dphidxsi[i];
                dydeta += Yloc[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }
			
            u = interpolate(phi,U,mapElem,n);
            v = interpolate(phi,V,mapElem,n);
			e = interpolate(phi,E,mapElem,n);
            y = interpolate(phi,Y,mapElem,n);
			f = f0 + beta*(Yloc[i] - L/2.0);
			tau = tau0*sin( M_PI*(Yloc[i]-L/2.0)/L );
            
            for (i=0; i < n; i++) {
				f = f0 + beta*(Yloc[i] - L/2.0);
				tau = tau0*sin( M_PI*(Yloc[i]-L/2.0)/L );
                BE[mapElem[i]] += (dphidx[i]*h*u + dphidy[i]*h*v)*jac*weight;
				BU[mapElem[i]] += (phi[i]*(f*v + tau/(rho*h) - gamma*u) + dphidx[i]*g*e)*jac*weight;
				BV[mapElem[i]] += (phi[i]*(-f*u -gamma*v) + dphidy[i]*g*e)*jac*weight;
			}
		}
	}
}

# endif
# ifndef NOINTEGRALSEDGES

void femShallowAddIntegralsEdges(femShallowProblem *myProblem)
{
    double *BE = myProblem->FE;
    double *BU = myProblem->FU;
    double *BV = myProblem->FV;
    double *E = myProblem->E;
    double *U = myProblem->U;
    double *V = myProblem->V;
	
    double  h = myProblem->height;
    double  g = myProblem->gravity;
	
    femIntegration *theRule = myProblem->rule1d;
    femDiscrete *theSpace = myProblem->space;
    
	int n = theSpace->n;
    double  xEdge[2],yEdge[2],phiEdge[2];
    double  xsi,weight,jac;
    double  ul,ur,vl,vr,el,er,uStar,eStar,unl,unr;
    int     i,j,k,edge,mapEdge[2][2];
	
    for (edge=0; edge < myProblem->edges->nEdge; edge++) {
        femShallowEdgeMap(myProblem,edge,mapEdge,n);
        for (j=0; j < 2; ++j) {
            int node = myProblem->edges->edges[edge].node[j];
            xEdge[j] = myProblem->mesh->X[node];
            yEdge[j] = myProblem->mesh->Y[node]; }
        double dxdxsi = xEdge[1] - xEdge[0];
        double dydxsi = yEdge[1] - yEdge[0];
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double normal[2] = {dydxsi/norm, -dxdxsi/norm};
        jac = norm / 2.0;
        for (k=0; k < theRule->n; k++) {
            xsi = theRule->xsi[k];
            weight = theRule->weight[k];
            femDiscretePhi1(theSpace,xsi,phiEdge);
            ul = interpolate(phiEdge,U,mapEdge[0],2);
			ur = interpolate(phiEdge,U,mapEdge[1],2);
            vl = interpolate(phiEdge,V,mapEdge[0],2);
			vr = interpolate(phiEdge,V,mapEdge[1],2);
			el = interpolate(phiEdge,E,mapEdge[0],2);
			er = interpolate(phiEdge,E,mapEdge[1],2);
			unl = ul*normal[0]+vl*normal[1];
			unr = ur*normal[0]+vr*normal[1];
			uStar = (unl+unr)/2.0 + sqrt(g/h)*(el-er)/2.0;
			eStar = (el+er)/2.0 + sqrt(h/g)*(unl-unr)/2.0;
			
            for (i=0; i < 2; i++) {
                BE[mapEdge[0][i]] -= (phiEdge[i]*h*uStar)*jac*weight;
                BE[mapEdge[1][i]] += (phiEdge[i]*h*uStar)*jac*weight;
                BU[mapEdge[0][i]] -= (phiEdge[i]*normal[0]*g*eStar)*jac*weight;
                BU[mapEdge[1][i]] += (phiEdge[i]*normal[0]*g*eStar)*jac*weight;
                BV[mapEdge[0][i]] -= (phiEdge[i]*normal[1]*g*eStar)*jac*weight;
                BV[mapEdge[1][i]] += (phiEdge[i]*normal[1]*g*eStar)*jac*weight;
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

