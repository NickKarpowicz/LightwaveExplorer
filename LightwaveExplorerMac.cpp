#include <stdio.h>
#include "LightwaveExplorerCoreCPU.h"
int main(int argc, char* argv[]){
    printf("hi.\n");
    if(argc<2){
        printf("Running default\n");
        mainCPU(2, "DefaultValues.ini");
        return 1;
    }
    mainCPU(argc,argv[1]);
    return 0;
}