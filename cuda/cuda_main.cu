#include<stdio.h>
#include<math.h>
#include<stdlib.h>

extern "C" {

#include"cuda_functions.h"

}

#define N 3600
#define GENERATION 500
#define MAX_THREAD_NUM 1024

// Kernel that calculates each creatures new state, depending on its neighbors
__global__ void calc_gen(int *in, int *neigh, int *out) {
	extern __shared__ int sha[];
	// Divide extern for readability
	int *state = sha;

	int tid = threadIdx.x;
	int creature_id = (blockIdx.x * blockDim.x) + threadIdx.x;
	// Calculate k (number of creatures in row)
	int k = sqrt((float) blockDim.x);

	// Save state of current creature
	state[tid] = in[creature_id];

	// Outside rows and columns
	// The total states are equal to the number of threads (creatures) in each
	// block (blockDim.x or k*k) so start the shared memory for the Outside
	// rows and columns after that
	int *up = &state[blockDim.x];
	int *down = &up[k];
	int *right = &down[k];
	int *left = &right[k];

	// Useful results for determining creature's position in each block
	int tid_div = tid / k;
	int tid_mod = tid % k;

	// Get outside rows and cols and save in shared
	// If its the upper row, then each upper thread (creature) stores the
	// neighboring creature that is up
	if (tid_div == 0) {
		up[tid_mod] = in[(neigh[blockIdx.x * 8] * blockDim.x) + blockDim.x - k + tid];
	}
	// If its the last row, then each thread (creature) stores the
	// neighboring creature that is down
	else if (tid_div == k - 1) {
		down[tid_mod] = in[(neigh[blockIdx.x * 8 + 1] * blockDim.x) + tid];
	}
	// If its the left column, then each thread (creature) stores the
	// neighboring creature that is left
	else if (tid_mod == 0) {
		left[tid_div] = in[(neigh[blockIdx.x * 8 + 2] * blockDim.x) + (tid_div*k + k - 1)];
	}
	// If its the right column, then each thread (creature) stores the
	// neighboring creature that is right
	else if (tid_mod == k - 1) {
		right[tid_div] = in[(neigh[blockIdx.x * 8 + 3] * blockDim.x) + (tid_div * k)];
	}

	int neighbor_count = 0;
	// Synchronise so that shared memory is all there
	__syncthreads();

	// Calculate next creature state
	if (tid_div == 0) {
		// For the upper left creature
		if (tid == 0) {
			int up_left = in[(neigh[blockIdx.x * 8 + 4] * blockDim.x) + k*k - 1];
			neighbor_count = up_left + up[0] + up[1] + left[0] + left[1] + state[1]
							+ state[k] + state[k + 1];
		}
		// For the upper right creature
		else if (tid == k - 1) {
			int up_right = in[(neigh[blockIdx.x * 8 + 5] * blockDim.x) + k * (k - 1)];
			neighbor_count = up_right + up[tid_mod] + up[tid_mod - 1] + right[0] + right[1]
							+ state[tid - 1] + state[tid + k - 1] + state[tid + k];
		}
		// For the upper row
		else {
			neighbor_count = up[tid_mod - 1] + up[tid_mod] + up[tid_mod + 1] + state[tid - 1]
							+ state[tid + 1] + state[tid + k - 1] + state[tid + k]
							+ state[tid + k + 1];
		}
	}
	else if (tid_div == k - 1) {
		// For the downmost left creature
		if (tid == k * (k - 1)) {
			int down_left = in[(neigh[blockIdx.x * 8 + 6] * blockDim.x) + k - 1];
			neighbor_count = down_left + left[tid_div - 1] + left[tid_div] + down[0] + down[1]
							+ state[tid + 1] + state[tid - k] + state[tid - k + 1];
		}
		// For the downmost right creature
		else if (tid == k*k - 1) {
			int down_right = in[(neigh[blockIdx.x * 8 + 7] * blockDim.x)];
			neighbor_count = down_right + right[tid_div - 1] + right[tid_div]
							+ down[k - 1] + down[k - 2] + state[tid - 1]
							+ state[tid - k] + state[tid - k - 1];
		}
		// For the downmost row
		else {
			neighbor_count = down[tid_mod - 1] + down[tid_mod] + down[tid_mod + 1] +
							+ state[tid - 1] + state[tid + 1]
							+ state[tid - k] + state[tid - k - 1] + state[tid - k + 1];
		}
	}
	// For the leftmost column
	else if (tid_mod == 0) {
		neighbor_count = left[tid_div] + left[tid_div - 1] + left[tid_div + 1]
						+ state[tid + 1] + state[tid - k] + state[tid - k + 1]
						+ state[tid + k] + state[tid + k + 1];
	}
	// For the rightmost column
	else if (tid_mod == k - 1) {
		neighbor_count = right[tid_div] + right[tid_div - 1] + right[tid_div + 1]
						+ state[tid - 1] + state[tid - k] + state[tid - k - 1]
						+ state[tid + k] + state[tid + k - 1];
	}
	// Normally inside the perimeter
	else {
		neighbor_count = state[tid - k - 1] + state[tid - k] + state[tid - k + 1]
						+ state[tid - 1] + state[tid + 1]
						+ state[tid + k - 1] + state[tid + k] + state[tid + k + 1];
	}

	// Finally calculate the new state of the creature
	if (state[tid] == 1) {
		if ((neighbor_count == 2) || (neighbor_count == 3)) {
			out[creature_id] = 1;
		}
		else {
			out[creature_id] = 0;
		}
	}
	// Else if state[tid] == 0
	else {
		if (neighbor_count == 2) {
			out[creature_id] = 1;
		}
		else {
			out[creature_id] = 0;
		}
	}
}

// Same kernel with added shared memory to check if there is no change between
// generations (called every 20 generations, heavier)
// Contains a form of reduction as a means to quickly determine if
// a change has been made in the world
// The changed data for each block is stored in *chagned
__global__ void calc_gen_wcheck(int *in, int *neigh, int *out, int *changed) {
	extern volatile __shared__ int sha_wcheck[];
	// Divide extern for readability
 	volatile int *state = sha_wcheck;

	int tid = threadIdx.x;
	int creature_id = (blockIdx.x * blockDim.x) + threadIdx.x;
	// Calculate k (number of creatures in row)
	int k = sqrt((float) blockDim.x);

	// Save state of current creature
	state[tid] = in[creature_id];

	// Outside rows and columns
	// The total states are equal to the number of threads (creatures) in each
	// block (blockDim.x or k*k) so start the shared memory for the Outside
	// rows and columns after that
	volatile int *up = &state[blockDim.x];
 	volatile int *down = &up[k];
 	volatile int *right = &down[k];
 	volatile int *left = &right[k];

	// Useful results for determining creature's position in each block
	int tid_div = tid / k;
	int tid_mod = tid % k;

	// Get outside rows and cols and save in shared
	// If its the upper row, then each upper thread (creature) stores the
	// neighboring creature that is up
	if (tid_div == 0) {
		up[tid_mod] = in[(neigh[blockIdx.x * 8] * blockDim.x) + blockDim.x - k + tid];
	}
	// If its the last row, then each thread (creature) stores the
	// neighboring creature that is down
	else if (tid_div == k - 1) {
		down[tid_mod] = in[(neigh[blockIdx.x * 8 + 1] * blockDim.x) + tid];
	}
	// If its the left column, then each thread (creature) stores the
	// neighboring creature that is left
	else if (tid_mod == 0) {
		left[tid_div] = in[(neigh[blockIdx.x * 8 + 2] * blockDim.x) + (tid_div*k + k - 1)];
	}
	// If its the right column, then each thread (creature) stores the
	// neighboring creature that is right
	else if (tid_mod == k - 1) {
		right[tid_div] = in[(neigh[blockIdx.x * 8 + 3] * blockDim.x) + (tid_div * k)];
	}

	int neighbor_count = 0;
	// Synchronise so that shared memory is all there
	__syncthreads();

	// Calculate next creature state
	if (tid_div == 0) {
		// For the upper left creature
		if (tid == 0) {
			int up_left = in[(neigh[blockIdx.x * 8 + 4] * blockDim.x) + k*k - 1];
			neighbor_count = up_left + up[0] + up[1] + left[0] + left[1] + state[1]
							+ state[k] + state[k + 1];
		}
		// For the upper right creature
		else if (tid == k - 1) {
			int up_right = in[(neigh[blockIdx.x * 8 + 5] * blockDim.x) + k * (k - 1)];
			neighbor_count = up_right + up[tid_mod] + up[tid_mod - 1] + right[0] + right[1]
							+ state[tid - 1] + state[tid + k - 1] + state[tid + k];
		}
		// For the upper row
		else {
			neighbor_count = up[tid_mod - 1] + up[tid_mod] + up[tid_mod + 1] + state[tid - 1]
							+ state[tid + 1] + state[tid + k - 1] + state[tid + k]
							+ state[tid + k + 1];
		}
	}
	else if (tid_div == k - 1) {
		// For the downmost left creature
		if (tid == k * (k - 1)) {
			int down_left = in[(neigh[blockIdx.x * 8 + 6] * blockDim.x) + k - 1];
			neighbor_count = down_left + left[tid_div - 1] + left[tid_div] + down[0] + down[1]
							+ state[tid + 1] + state[tid - k] + state[tid - k + 1];
		}
		// For the downmost right creature
		else if (tid == k*k - 1) {
			int down_right = in[(neigh[blockIdx.x * 8 + 7] * blockDim.x)];
			neighbor_count = down_right + right[tid_div - 1] + right[tid_div]
							+ down[k - 1] + down[k - 2] + state[tid - 1]
							+ state[tid - k] + state[tid - k - 1];
		}
		// For the downmost row
		else {
			neighbor_count = down[tid_mod - 1] + down[tid_mod] + down[tid_mod + 1] +
							+ state[tid - 1] + state[tid + 1]
							+ state[tid - k] + state[tid - k - 1] + state[tid - k + 1];
		}
	}
	// For the leftmost column
	else if (tid_mod == 0) {
		neighbor_count = left[tid_div] + left[tid_div - 1] + left[tid_div + 1]
						+ state[tid + 1] + state[tid - k] + state[tid - k + 1]
						+ state[tid + k] + state[tid + k + 1];
	}
	// For the rightmost column
	else if (tid_mod == k - 1) {
		neighbor_count = right[tid_div] + right[tid_div - 1] + right[tid_div + 1]
						+ state[tid - 1] + state[tid - k] + state[tid - k - 1]
						+ state[tid + k] + state[tid + k - 1];
	}
	// Normally inside the perimeter
	else {
		neighbor_count = state[tid - k - 1] + state[tid - k] + state[tid - k + 1]
						+ state[tid - 1] + state[tid + 1]
						+ state[tid + k - 1] + state[tid + k] + state[tid + k + 1];
	}

	// Finally calculate the new state of the creature
	int new_state = 0;
	int prev_state = state[tid];
	if (prev_state == 1) {
		if ((neighbor_count == 2) || (neighbor_count == 3)) {
			new_state = 1;
		}
		else {
			new_state = 0;
		}
	}
	// Prev_state == 0
	else {
		if (neighbor_count == 2) {
			new_state = 1;
		}
		else {
			new_state = 0;
		}
	}
	out[creature_id] = new_state;
	// Sync because the shared memory state will be repurposed
	__syncthreads();
	// Now use the state shared memory to save if the state of current creature has changed
	// Save 0 if the creature's state has not changed, 1 otherwise after calculating
	// the creature's new state
	if (new_state != prev_state) {
		state[tid] = 1;
	}
	else {
		state[tid] = 0;
	}
	// Begin custom reduction loop
	for(unsigned int s = blockDim.x/2; s > 32; s>>=1) {
		if (tid < s)
			state[tid] += state[tid + s];
		__syncthreads();
	}
	// Unroll last six iterations (same warp)
	if (tid < 32) {
		state[tid] += state[tid + 32];
		state[tid] += state[tid + 16];
		state[tid] += state[tid +  8];
		state[tid] += state[tid +  4];
		state[tid] += state[tid +  2];
		state[tid] += state[tid +  1];
	}
	// Write reduction result back to global memory
	if (tid == 0)
		changed[blockIdx.x] = state[0];
}

// Cuda kernel that calculates, for each block all its neighboring blocks (8)
// and outputs it in the out array
// (nearly the same code as the mpi and mpi_openmp implementations)
__global__ void calc_neighbors(int sq_of_P, int P_size , int *out) {
	int rank = (blockIdx.x * blockDim.x) + threadIdx.x;
	// Some threads will exceed the wanted rank outputs
	// (P_size % MAX_THREAD_NUM)
	if (rank < P_size) {

		int Pmod = rank % sq_of_P;
		int Pdiv = rank / sq_of_P;  //rank = div*sP+mod

		int D = ((Pdiv + 1) % sq_of_P)*sq_of_P + Pmod;
		int R  = Pdiv*sq_of_P + (Pmod + 1)%sq_of_P;
		int Dr = ((Pdiv + 1)%sq_of_P)*sq_of_P + (Pmod + 1)%sq_of_P;
		int U, L, Ul, Ur, Dl;

		if (Pdiv == 0) {
			U  = P_size - sq_of_P + Pmod;
			Ur = P_size - sq_of_P + (Pmod + 1)%sq_of_P;
			if (Pmod == 0)
				Ul = P_size - 1;
			else
				Ul = P_size - sq_of_P + Pmod - 1;
		}
		else {
			U  = (Pdiv - 1)*sq_of_P + Pmod;
			Ur = (Pdiv - 1)*sq_of_P + (Pmod + 1)%sq_of_P;
	    	if (Pmod == 0)
				Ul = Pdiv*sq_of_P - 1;
	        else
	        	Ul = (Pdiv - 1)*sq_of_P + Pmod - 1;
		}
		if (Pmod == 0) {
			L  = (Pdiv + 1)*sq_of_P - 1;
			Dl =((Pdiv + 1)%sq_of_P)*sq_of_P + sq_of_P - 1;
		}
		else {
			L  = Pdiv*sq_of_P + Pmod - 1;
			Dl = ((Pdiv + 1)%sq_of_P)*sq_of_P + Pmod - 1;
		}

		// Store results
		out[rank * 8] = U;
		out[rank * 8 + 1] = D;
		out[rank * 8 + 2] = L;
		out[rank * 8 + 3] = R;
		out[rank * 8 + 4] = Ul;
		out[rank * 8 + 5] = Ur;
		out[rank * 8 + 6] = Dl;
		out[rank * 8 + 7] = Dr;
	}
}

// Main
int main(int argc, char* argv[]) {
	// P_size is the number of blocks that the problem will be split into
	int P_size;
	if (argc >= 2)
		P_size = atoi(argv[1]);
	else
		return 1;

	int sqrt_P = sqrt(P_size);
	// sq_of_P=sqrt of P_size and k = N / sq_of_P (must be integers)
	int k  = N / sqrt_P;

	// k*k (creature number) should not be greater than the maximum thread number
	// of a block (depends of the gpu)
	if (k*k > MAX_THREAD_NUM) {
		printf("Please insert higher block number, so that each block has less than or equal to %d creatures \n", MAX_THREAD_NUM);
		return 1;
	}
	int world_datasize = N*N * sizeof(int);
	int *world = (int*)malloc(world_datasize);
	// Later used to check if world changed
	int *newWorld = (int*)malloc(world_datasize);

	// Create the world
	if (argc == 3) { //from text
		int *inArray = (int*)malloc(N*N * sizeof(int));
		getInputData(argv[4], inArray, N*N);
		makeArrayInBlocks(inArray, world, sqrt_P, k);
		free(inArray);
	}
	else
		makeRandomWorld(world, N*N);

	// Calculate neighbors
	int *d_neigh;
	cudaMalloc((void**) &d_neigh, 8 * P_size * sizeof(int));
	// Divide problem use as many threads as possible
	calc_neighbors<<<MAX_THREAD_NUM, (P_size / MAX_THREAD_NUM) + 1>>>(sqrt_P, P_size, d_neigh);

	// Initialize main loop cuda memory
	int *d_world, *d_newWorld;
	cudaMalloc((void**) &d_world, world_datasize);
	cudaMalloc((void**) &d_newWorld, world_datasize);

	// A value stores 1 if the block that it represents has
	// changed, 0 otherwise
	int *changed = (int*)malloc(P_size * sizeof(int));
	for (int i = 0; i < P_size; i++)
		changed[i] = 0;

	int *d_changed;
	cudaMalloc((void**) &d_changed, P_size * sizeof(int));
	cudaMemcpy(d_changed, changed, P_size * sizeof(int), cudaMemcpyHostToDevice);

	// Copy main world data
	cudaMemcpy(d_world, &world, world_datasize, cudaMemcpyHostToDevice);
	// Calculate shared memory size
	// sha_size = k*k (block creatures) + 4*k (4 neighboring rows and columns) * sizeof(int)
	int changed_flag = 0;
	for (int i = 0; i < GENERATION; i++) {
		// Check if world didn't change every 20 generations
		if( !((i+1) % 20) ) {
			// Mind the shared memory size
			calc_gen_wcheck <<<P_size, k*k, (k*k + 4*k) * sizeof(int)>>> (d_world, d_neigh, d_newWorld, d_changed);
			cudaMemcpy(changed, d_changed, P_size * sizeof(int), cudaMemcpyDeviceToHost);
			for (int ins = 0; ins < P_size; i++) {
				if (changed[ins] > 0) {
					changed_flag = 1;
					break;
				}
			}
			if(changed_flag == 0) {
				printf("Grid didn't change in %d generation! Exiting..\n", i + 1);
				break;
			}
			else {
				changed_flag = 0;
				// Reinitialise the d_changed array
				for (int i = 0; i < P_size; i++)
					changed[i] = 0;
				cudaMemcpy(d_changed, changed, P_size * sizeof(int), cudaMemcpyHostToDevice);
			}
		}
		// Else don't check or change global variables
		else {
			calc_gen <<<P_size, k*k, (k*k + 4*k) * sizeof(int)>>> (d_world, d_neigh, d_newWorld);
		}
		// Swap arrays for next generation iteration
		swapArrays(&d_world, &d_newWorld);
	}

	// Commented out, used for printing
	//cudaMemcpy(newWorld, &d_newWorld, world_datasize, cudaMemcpyDeviceToHost);
	// printworld(newWorld, sqrt_P, k);

	// Free memory
	free(world);
	free(newWorld);
	free(changed);
	// Free allocated cuda memory
	cudaFree(d_neigh);
	cudaFree(d_world);
	cudaFree(d_newWorld);
	cudaFree(d_changed);

	return 0;
}
