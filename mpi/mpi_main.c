#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mpi.h>
#include"mpi_functions.h"

#define N 360
#define GENERATION 200

int main(int argc, char* argv[]) {
	int sqrt_P, k; //sq_of_P=sqrt of P_size and k=N/sq_of_P (must be integers)
	int *inArray = NULL, *world = NULL, *newWorld = NULL;
	int *myWorld = NULL, *myNewWorld = NULL;
	int i = 0, j = 0, z = 0;
	int local_changed = 0, total_changed = 0;
	int P_size, P_rank;
	int up, down, right, left, up_left, up_right, down_left, down_right; //the 8 neighbors of every P
	int *recup = NULL, *recdown = NULL, *recright = NULL, *recleft = NULL; //pointers and
	int recup_left, recup_right, recdown_left, recdown_right; //ints to receive data from neighbors
	MPI_Datatype LineType;
	MPI_Datatype ColType;
  	MPI_Status fstatus;
  	MPI_Offset offset;
	MPI_Status status[8];
	MPI_Request *send_req, *receive_req; //send_requestX,receive_requestX
  	MPI_File filehandle;

  	double localstart,localfinish,localelapsed,elapsed;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &P_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &P_rank);

  	MPI_Barrier(MPI_COMM_WORLD);
  	localstart=MPI_Wtime();

	sqrt_P = sqrt(P_size);
	k = N / sqrt_P;

    //check if (number of processes) or (problem's dimension) are accepted numbers
  	if((sqrt_P*sqrt_P != P_size) || ((N % sqrt_P) != 0)){//otherwise program stops here
    	if(P_rank==0)
      		printf("Wrong argument(s)!\n");
    	return -1;
  	}

	send_req = malloc(8 * sizeof(MPI_Request));
	receive_req = malloc(8 * sizeof(MPI_Request));
	myWorld = malloc(k*k * sizeof(int));
	myNewWorld = malloc(k*k * sizeof(int));
	recup = malloc(k * sizeof(int));
	recdown = malloc(k * sizeof(int));
	recright = malloc(k * sizeof(int));
	recleft = malloc(k * sizeof(int));

  	offset = P_rank*k;
  	MPI_File_open(MPI_COMM_WORLD,"input.txt",MPI_MODE_CREATE|MPI_MODE_DELETE_ON_CLOSE|MPI_MODE_RDWR,MPI_INFO_NULL,&filehandle);
  	MPI_File_set_view(filehandle,0,MPI_INT,MPI_INT,"native",MPI_INFO_NULL);

	if (P_rank == 0) {
    	world = malloc(N*N * sizeof(int));
		// Uncomment for PRINT (if you comment out print below)
        //newWorld = malloc(N*N*sizeof(int));
		//create the world
		if (argc == 5) { //from text
			inArray = malloc(N*N * sizeof(int));
			getInputData(argv[4], inArray, N*N);
			makeArrayInBlocks(inArray, world, sqrt_P, k);
			free(inArray);
		}
		else
			makeRandomWorld(world, N*N);
		//split data to all processes
		//MPI_Scatter(world, k*k, MPI_INT, myWorld, k*k, MPI_INT, 0, MPI_COMM_WORLD);
        //P 0 writes data in shared file
    	MPI_File_write_shared(filehandle,world,N*N,MPI_INT,&fstatus);
	}
	//else
		//MPI_Scatter(world, k*k, MPI_INT, myWorld, k*k, MPI_INT, 0, MPI_COMM_WORLD);

  	MPI_Barrier(MPI_COMM_WORLD);//wait P 0 to complete writing and then read
  	MPI_File_read_at(filehandle,offset,myWorld,k,MPI_INT,&fstatus);

	findNeighbors(P_rank, sqrt_P, P_size, &up, &down, &right, &left, &up_left, &up_right, &down_left, &down_right);

	//create datatypes for line and column;
	MPI_Type_contiguous(k, MPI_INT, &LineType);
	MPI_Type_commit(&LineType);
	MPI_Type_vector(k, 1, k, MPI_INT, &ColType);
	MPI_Type_commit(&ColType);

	for (i = 0; i < GENERATION; i++) {

		//send * 8
    MPI_Isend(&myWorld[0], 1, MPI_INT, up_left, 0, MPI_COMM_WORLD, &send_req[0]);	//ul
    MPI_Isend(&myWorld[k - 1], 1, MPI_INT, up_right, 0, MPI_COMM_WORLD, &send_req[1]);	//ur
    MPI_Isend(&myWorld[k*(k - 1)], 1, MPI_INT, down_left, 0, MPI_COMM_WORLD, &send_req[2]);	//dl
    MPI_Isend(&myWorld[k*k - 1], 1, MPI_INT, down_right, 0, MPI_COMM_WORLD, &send_req[3]);	//dr
    MPI_Isend(&myWorld[0], 1, LineType, up, 0, MPI_COMM_WORLD, &send_req[4]);	//u
    MPI_Isend(&myWorld[k*(k - 1)], 1, LineType, down, 0, MPI_COMM_WORLD, &send_req[5]);	//d
    MPI_Isend(&myWorld[0], 1, ColType, left, 0, MPI_COMM_WORLD, &send_req[6]);	//l
    MPI_Isend(&myWorld[k-1], 1, ColType, right, 0, MPI_COMM_WORLD, &send_req[7]);	//r

		//recv*8
    MPI_Irecv(&recup_left, 1, MPI_INT, up_left, 0, MPI_COMM_WORLD, &receive_req[0]);	//ul
    MPI_Irecv(&recup_right, 1, MPI_INT, up_right, 0, MPI_COMM_WORLD, &receive_req[1]);	//ur
    MPI_Irecv(&recdown_left, 1, MPI_INT, down_left, 0, MPI_COMM_WORLD, &receive_req[2]);	//dl
    MPI_Irecv(&recdown_right, 1, MPI_INT, down_right, 0, MPI_COMM_WORLD, &receive_req[3]);	//dr
    MPI_Irecv(recup, 1, LineType, up, 0, MPI_COMM_WORLD, &receive_req[4]);	//u
    MPI_Irecv(recdown, 1, LineType, down, 0, MPI_COMM_WORLD, &receive_req[5]);	//d
    MPI_Irecv(recleft, k, MPI_INT, left, 0, MPI_COMM_WORLD, &receive_req[6]);	//l
    MPI_Irecv(recright, k, MPI_INT, right, 0, MPI_COMM_WORLD, &receive_req[7]);	//r


	checkInside(myWorld, myNewWorld, k); //save myWorld's changes in myNewWorld

	MPI_Waitall(8, send_req, status);
	MPI_Waitall(8, receive_req, status);

	checkPerimeter(myWorld, k, recdown_right, recdown_left, recup_right, recup_left, recdown, recup, recright, recleft, myNewWorld);//same here ^

	// check if grid changed or not...
	if( !((i + 1) % 20) ) {	//check if world didn't changed every 20 times
   		if( !(world_changed(myWorld, myNewWorld, k*k)) )
     		local_changed=1;

		MPI_Allreduce(&local_changed, &total_changed,1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	 	if(total_changed == P_size) {
   			if(P_rank==0)
     			printf("Grid didn't change in %d generation! Exiting..\n",i+1);
	 		break;
	  	}
 	}

		// Uncomment if you want to print the current two consecutive generations
 /*
 		if(P_rank == 0){//print current and new grid
			MPI_Gather(myWorld,k*k,MPI_INT,world,k*k,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Gather(myNewWorld,k*k,MPI_INT,newWorld,k*k,MPI_INT,0,MPI_COMM_WORLD);
		}
		else{
			MPI_Gather(myWorld,k*k,MPI_INT,NULL,k*k,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Gather(myNewWorld,k*k,MPI_INT,NULL,k*k,MPI_INT,0,MPI_COMM_WORLD);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		if(P_rank == 0){
			printWorld(world,sqrt_P,k);

            printWorld(newWorld,sqrt_P,k);
		}
*/

		swapArrays(&myWorld, &myNewWorld);
    local_changed = 0;
	}

  localfinish = MPI_Wtime();
  localelapsed = localfinish - localstart;
  MPI_Reduce(&localelapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if(P_rank == 0)
  	printf("Execution time: %f\n",elapsed);

	free(recup);
	free(recdown);
	free(recleft);
	free(recright);
	free(myWorld);
	free(myNewWorld);
	free(world);
  free(newWorld);
	free(receive_req);
	free(send_req);
  MPI_File_close(&filehandle);
	MPI_Type_free(&LineType);
	MPI_Type_free(&ColType);
	MPI_Finalize();

	return 0;
}
