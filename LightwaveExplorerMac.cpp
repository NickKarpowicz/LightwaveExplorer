#include <stdio.h>
#include "LightwaveExplorerCoreCPU.h"
int main(int argc, char* argv[]){
    if(argc<2){
        printf("Running default\n");
        mainCPU(2, "DefaultValues.ini");
        return 1;
    }
    mainCPU(argc,argv[1]);
    return 0;
}