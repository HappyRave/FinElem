
#include <stdio.h>
#include <stdlib.h>

void tsunamiAnimate(double dt, int nmax, int sub, int order, const char *meshFileName, const char *baseResultName);


int main(void)
{   
        tsunamiAnimate(1.0,100,25,1,"tsunamiSmall.txt","data/tsunamiSmall");
        exit(0);
       
}