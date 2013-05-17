
#include <stdio.h>
#include <stdlib.h>

void tsunamiCompute(double dt, int nmax, int sub, int order, const char *meshFileName, const char *baseResultName);
void tsunamiAnimate(double dt, int nmax, int sub, int order, const char *meshFileName, const char *baseResultName);
void guiSandbox(double dt, int nmax, int sub, int order, const char *meshFileName, const char *baseResultName);


int main(void)
{
	//tsunamiCompute(1.0,24000,100,2,"PacificQuadFine.txt","data/tsunamiQuad2Fine");

	//tsunamiAnimate(1.0,1000,100,1,"PacificQuadSmall.txt","data/tsunamiSmall");
	tsunamiAnimate(1.0,24000,100,2,"PacificQuadFine.txt","data/tsunamiQuad2Fine");
	//guiSandbox(1.0,4000,50,2,"PacificQuadFine.txt","data/tsunamiSmall");
	
	exit(0);

    
}