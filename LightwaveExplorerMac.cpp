#include <stdio.h>
#include "LightwaveExplorerCoreCPU.h"
int main(int argc, char* argv[]){
    if(argc<2){
        printf("Running default\n");
        char defaultPath[] = "DefaultValues.ini\0";
        mainCPU(2, defaultPath);
        return 1;
    }
    mainCPU(argc,argv[1]);
    return 0;
}