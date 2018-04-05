#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"serial_functions.h"

#define N 3600
#define GENERATION 500

int main(int argc,char* argv[]){
	int *world=NULL,*newWorld=NULL;
	int i,j;

	clock_t begin=clock(),end;

	world = malloc(N*N*sizeof(int));
	newWorld = malloc(N*N*sizeof(int));

	if(argc==2)
		getInputData(argv[1],world,N*N);
	else
		makeRandomWorld(world,N*N);

//	printWorld(world,N);

	for(i=0;i<GENERATION;i++){
		computeNewWorld(world,newWorld,N);
		if((i+1)%20==0)
			if(!world_changed(world,newWorld,N*N))
				break;
		swapArrays(&world,&newWorld);
//		putchar('\n');
//		printWorld(world,N);
	}

	end = clock();
	printf("time:%lf\n",(double)(end-begin)/CLOCKS_PER_SEC);

	free(world);
	free(newWorld);
	return 0;
}
