
#include"fem.h"


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


void diffusionMeshLocal(const femMesh *theMesh, const int iElem, int *map, double *x, double *y)
{
    int j,nLocal = theMesh->nLocalNode;
    
    for (j=0; j < nLocal; ++j) {
        map[j] = theMesh->elem[iElem*nLocal+j];
        x[j]   = theMesh->X[map[j]];
        y[j]   = theMesh->Y[map[j]]; }
}

void diffusionSolve(femDiffusionProblem *theProblem)
{
    femMesh *theMesh = theProblem->mesh;
    femEdges *theEdges = theProblem->edges;
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
	
    if (theSpace->n > 4) Error("Unexpected discrete space size !");
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,map[4];
	
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        diffusionMeshLocal(theMesh,iElem,map,x,y);
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
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i]*dphidxsi[i];
                dxdeta += x[i]*dphideta[i];
                dydxsi += y[i]*dphidxsi[i];
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }
            for (i = 0; i < theSpace->n; i++) {
                for(j = 0; j < theSpace->n; j++) {
                    theSystem->A[map[i]][map[j]] += (dphidx[i] * dphidx[j]
													 + dphidy[i] * dphidy[j]) * jac * weight; }}
            for (i = 0; i < theSpace->n; i++) {
                theSystem->B[map[i]] += phi[i] * jac *weight; }}}
	
	for (iEdge= 0; iEdge < theEdges->nEdge; iEdge++) {
        if (theEdges->edges[iEdge].elem[1] < 0) {
            femFullSystemConstrain(theSystem,theEdges->edges[iEdge].node[0],0.0);
            femFullSystemConstrain(theSystem,theEdges->edges[iEdge].node[1],0.0); }}
	
    femFullSystemEliminate(theSystem);
}